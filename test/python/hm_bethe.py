import numpy as np
import scipy as sp
from edipy import *
import edipy as ed



def dens_bethe(x,d):
  root=(1-(x/d)**2)+0.j
  root=np.sqrt(root)
  dens_bethe=(2/(np.pi*d))*root
  return dens_bethe.real



wband=1
Le=1000

ed.read_input("inputED.conf")



#Build Bethe DOS:
Eband,de = np.linspace(-wband,wband,Le,retstep=True)
Dband = dens_bethe(Eband,wband)


#Allocate minimal required arrays:
Smats=np.zeros((ed.Nspin,ed.Nspin,ed.Norb,ed.Norb,ed.Lmats),dtype='complex',order='F')
Gmats=np.zeros((ed.Nspin,ed.Nspin,ed.Norb,ed.Norb,ed.Lmats),dtype='complex',order='F')
Sreal=np.zeros((ed.Nspin,ed.Nspin,ed.Norb,ed.Norb,ed.Lreal),dtype='complex',order='F')
Greal=np.zeros((ed.Nspin,ed.Nspin,ed.Norb,ed.Norb,ed.Lreal),dtype='complex',order='F')
Weiss=np.zeros((ed.Nspin,ed.Nspin,ed.Norb,ed.Norb,ed.Lmats),dtype='complex',order='F')
Hloc =np.zeros((ed.Nspin,ed.Nspin,ed.Norb,ed.Norb),dtype='float',order='F')



#Get bath dimension:
Nb = ed.get_bath_dimension()
bath = np.zeros(Nb,dtype='float',order='F')
bath_prev = np.zeros(Nb,dtype='float',order='F')

#INIT SOLVER
ed.init_solver(bath,Hloc)

wm = np.pi/ed.beta*(2*np.arange(ed.Lmats)+1)
wr = np.linspace(ed.wini,ed.wfin,ed.Lreal)

#START DMFT LOOP:
converged=False;iter=0
while (not converged and iter<ed.Nloop ):
    iter=iter+1
    print("DMFT-loop:",iter,"/",ed.Nloop)
    
    ed.solve(bath,Hloc)
    ed.get_sigma_matsubara(Smats)
    
    for i in range(ed.Lmats):
        Gmats[0,0,0,0,i] = sum(Dband/(wm[i]*1j - Smats[0,0,0,0,i] - Eband))*de
    np.savetxt('Gloc.dat', np.transpose([wm,Gmats[0,0,0,0,:].imag,Gmats[0,0,0,0,:].real]))

    Weiss[0,0,0,0,:] = 0.25*wband*Gmats[0,0,0,0,:] #1.0/(1.0/Gmats[0,0,0,0,:] + Smats[0,0,0,0,:])
    np.savetxt('Weiss.dat', np.transpose([wm,Weiss[0,0,0,0,:].imag,Weiss[0,0,0,0,:].real]))

    if(iter>1):
        bath = 0.5*bath + 0.5*bath_prev
    bath_prev=bath

    ed.chi2_fitgf(Weiss,bath,Hloc,ispin=1,iorb=1)
    err,converged=ed.check_convergence(Weiss[0,0,0,0,:],ed.dmft_error,1,ed.Nloop)

ed.get_sigma_realaxis(Sreal)
for i in range(ed.Lreal):
    Greal[0,0,0,0,i] = sum(Dband/(wr[i]+ed.eps*1j - Sreal[0,0,0,0,i] - Eband))*de
np.savetxt('Greal.dat', np.transpose([wr,Greal[0,0,0,0,:].imag,Greal[0,0,0,0,:].real]))



