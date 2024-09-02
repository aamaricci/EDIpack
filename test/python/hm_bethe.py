import numpy as np
import scipy as sp
from edipy import global_env as ed
import mpi4py
from mpi4py import MPI
import os,sys

#INIT MPI 
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
print("I am process",rank,"of",comm.Get_size())
master = (rank==0)

def dens_bethe(x,d):
  root=(1-(x/d)**2)+0.j
  root=np.sqrt(root)
  dens_bethe=(2/(np.pi*d))*root
  return dens_bethe.real

#READ ED INPUT:
ed.read_input("inputED.conf")
if(ed.Nspin!=1 or ed.Norb!=1):
    print("This test code is for Nspin=1 + Norb=1.")
    ed.Nspin=1
    ed.Norb=1
Le      = 1000
wmixing = 0.3
wband   = 1.0

#BUILD Density of States:
Eband,de = np.linspace(-wband,wband,Le,retstep=True)
Dband = dens_bethe(Eband,wband)

#BUILD frequency array:
wm = np.pi/ed.beta*(2*np.arange(ed.Lmats)+1)
wr = np.linspace(ed.wini,ed.wfin,ed.Lreal)


#Allocate minimal required arrays:
#ALL functions must have shape [Nspin,Nspin,Norb,Norb(,L)]:
Smats=np.zeros((ed.Nspin,ed.Nspin,ed.Norb,ed.Norb,ed.Lmats),dtype='complex',order='F')
Gmats=np.zeros((ed.Nspin,ed.Nspin,ed.Norb,ed.Norb,ed.Lmats),dtype='complex',order='F')
Sreal=np.zeros((ed.Nspin,ed.Nspin,ed.Norb,ed.Norb,ed.Lreal),dtype='complex',order='F')
Greal=np.zeros((ed.Nspin,ed.Nspin,ed.Norb,ed.Norb,ed.Lreal),dtype='complex',order='F')
Delta=np.zeros((ed.Nspin,ed.Nspin,ed.Norb,ed.Norb,ed.Lmats),dtype='complex',order='F')
Hloc =np.zeros((ed.Nspin,ed.Nspin,ed.Norb,ed.Norb),dtype='float',order='F')



#SETUP SOLVER
Nb=ed.get_bath_dimension()
bath = np.zeros(Nb,dtype='float',order='F')
ed.init_solver(bath)
bath_prev = np.copy(bath)

#DMFT CYCLE
converged=False;iloop=0
while (not converged and iloop<ed.Nloop ):
    iloop=iloop+1
    print("DMFT-loop:",iloop,"/",ed.Nloop)

    #Solve the EFFECTIVE IMPURITY PROBLEM (first w/ a guess for the bath)
    ed.solve(bath,Hloc)
    

    #Get Self-energies
    ed.get_sigma_matsubara(Smats)
    ed.get_sigma_realaxis(Sreal)

    #Compute the local gf:
    for i in range(ed.Lmats):
        zeta             = wm[i]*1j - Smats[0,0,0,0,i]
        Gmats[0,0,0,0,i] = sum(Dband/(zeta - Eband))*de

    for i in range(ed.Lreal):
        zeta             = wr[i]+ed.eps*1j - Sreal[0,0,0,0,i]
        Greal[0,0,0,0,i] = sum(Dband/(zeta - Eband))*de
        
    if(rank==0):
        np.savetxt('Gloc_iw.dat', np.transpose([wm,Gmats[0,0,0,0,:].imag,Gmats[0,0,0,0,:].real]))
        np.savetxt('Gloc_realw.dat', np.transpose([wr,Greal[0,0,0,0,:].imag,Greal[0,0,0,0,:].real]))


    #Get the Delta function and FIT:    
    Delta[0,0,0,0,:] = 0.25*wband*Gmats[0,0,0,0,:]
    if(rank==0):
        np.savetxt('Delta_iw.dat', np.transpose([wm,Delta[0,0,0,0,:].imag,Delta[0,0,0,0,:].real]))
    ed.chi2_fitgf(Delta,bath,ispin=1,iorb=1)
    
    if(iloop>1):
        bath = wmixing*bath + (1.0-wmixing)*bath_prev
    bath_prev=np.copy(bath)
 
    err,converged=ed.check_convergence(Delta[0,0,0,0,:],ed.dmft_error,1,ed.Nloop)

print("Done...")

