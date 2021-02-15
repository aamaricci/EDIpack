import numpy as np
import ed2py                         #import the Fortran module
import lancpy                          #import the PYTHON library itself
from ed2py import *
from lancpy import *


global Norb         #Norb =# of impurity orbitals
global Nspin        #Nspin=# spin degeneracy (max 2)
global Nloop        #max dmft loop variables
global xmu          #chemical potential
global beta         #inverse temperature
global w0_ph        #w0_ph: phonon frequency (constant)
global eps          #broadening
global wini,wfin    #frequency range
global xmin,xmax    #x-range for the local lattice probability distribution function (phonons)
global Nsuccess     #Number of repeated success to fall below convergence threshold  
global dmft_error   #dmft convergence threshold
global sb_field     #symmetry breaking field  
global cg_scheme    #fit scheme: delta (default), weiss for G0
global bath_type    #flag to set bath type: normal (1bath/imp), hybrid(1bath)
global nread        #fixed density. if 0.d0 fixed chemical potential calculation.
global Lmats        #Number of Matsuabra Freq.
global Lreal        #Number of real-axis freq.
global Lpos         #


#FROM ED_INPUT_VARS
Norb       = ed2py.ed_input_vars.norb
Nspin      = ed2py.ed_input_vars.nspin
Nloop      = ed2py.ed_input_vars.nloop
xmu        = ed2py.ed_input_vars.xmu
beta       = ed2py.ed_input_vars.beta
w0_ph      = ed2py.ed_input_vars.w0_ph
eps        = ed2py.ed_input_vars.eps
wini       = ed2py.ed_input_vars.wini
wfin       = ed2py.ed_input_vars.wfin
xmin       = ed2py.ed_input_vars.xmin
xmax       = ed2py.ed_input_vars.xmax
Nsuccess   = ed2py.ed_input_vars.nsuccess
dmft_error = ed2py.ed_input_vars.dmft_error
sb_field   = ed2py.ed_input_vars.sb_field
cg_scheme  = ed2py.ed_input_vars.cg_scheme
bath_type  = ed2py.ed_input_vars.bath_type
nread      = ed2py.ed_input_vars.nread
Lmats      = ed2py.ed_input_vars.lmats
Lreal      = ed2py.ed_input_vars.lreal
Lpos       = ed2py.ed_input_vars.lpos




#FROM ED_HLOC_DECOMPOSITION
def set_hloc(hloc,lam=None):
    if lam is None:
        ln=len(np.shape(hloc))
        if(ln==4):
            ed2py.set_hloc_nn(hloc)
        elif(ln==2):
            ed2py.set_hloc_so(hloc)
        else:
            raise
        ValueError('shape(Hloc) /= [Nspin,Nspin,Norb,Norb] or [Nspin*Norb,Nspin*Norb]')
    else:
        ed2py.set_hloc_symm(hloc,lam)
    return;




#FROM ED_BATH:
def get_bath_dimension(hloc=None):    
    if not hloc:
        Nb = ed2py.get_bath_dimension_direct()
    else:
        assert isinstance(hloc, complex), 'type(hloc)/=complex'
        DimHloc=len(np.shape(hloc))
        if(ln==5):
            Nb = ed2py.get_bath_dimension_symmetries(hloc)
        else:
            raise ValueError('Dim(Hloc) /= [Nsym,Nspin,Nspin,Norb,Norb]')
    return Nb;

# missing interface in ed2py:
# get_bath_component_dimension
# get_bath_component
# set_bath_component
# copy_bath_component

def spin_symmetrize_bath(bath,save):
    assert isinstance(save, bool), 'type(save)/=bool'
    DimBath=len(np.shape(bath))
    if(DimBath==1):
        ed2py.spin_symmetrize_bath_site(bath,save)
    elif(DimBath==2):
        ed2py.spin_symmetrize_bath_ineq(bath,save)
    else:
        raise ValueError('Dim(Bath) /= [{Nlat},Nb]')
    return bath;

def orbs_symmetrize_bath(bath,save:bool):
    assert isinstance(save, bool), 'type(save)/=bool'
    DimBath=len(np.shape(bath))
    if(DimBath==1):
        ed2py.orbs_symmetrize_bath_site(bath,save)
    elif(DimBath==2):
        ed2py.orbs_symmetrize_bath_ineq(bath,save)
    else:
        raise ValueError('Dim(Bath) /= [{Nlat},Nb]')
    return bath;

def orb_equality_bath(bath,save:bool):
    assert isinstance(save, bool), 'type(save)/=bool'
    DimBath=len(np.shape(bath))
    if(DimBath==1):
        ed2py.orb_equality_bath_site(bath,save)
    elif(DimBath==2):
        ed2py.orb_equality_bath_ineq(bath,save)
    else:
        raise ValueError('Dim(Bath) /= [{Nlat},Nb]')
    return bath;

def ph_symmetrize_bath(bath,save:bool):
    assert isinstance(save, bool), 'type(save)/=bool'
    DimBath=len(np.shape(bath))
    if(DimBath==1):
        ed2py.ph_symmetrize_bath_site(bath,save)
    elif(DimBath==2):
        ed2py.ph_symmetrize_bath_ineq(bath,save)
    else:
        raise ValueError('Dim(Bath) /= [{Nlat},Nb]')
    return bath;

def ph_trans_bath(bath,save:bool):
    assert isinstance(save, bool), 'type(save)/=bool'
    DimBath=len(np.shape(bath))
    if(DimBath==1):
        ed2py.ph_trans_bath_site(bath,save)
    elif(DimBath==2):
        ed2py.ph_trans_bath_ineq(bath,save)
    else:
        raise ValueError('Dim(Bath) /= [{Nlat},Nb]')
    return bath;

def break_symmetry_bath(bath,field,sign:float,save:bool):
    assert isinstance(field, float), 'type(field)/=float'
    assert isinstance(sign, float), 'type(sign)/=float'
    assert isinstance(save, bool), 'type(save)/=bool'
    DimBath=len(np.shape(bath))
    if(DimBath==1):
        bath=ed2py.break_symmetry_bath_site(field,sign,save)
    elif(DimBath==2):
        bath=ed2py.break_symmetry_bath_ineq(field,sign,save)
    else:
        raise ValueError('Dim(Bath) /= [{Nlat},Nb]')
    return bath;




# #FROM ED_AUX_FUNX:
# missing interface in ed2py:
# ed_search_variable





#FROM ED_IO:
def get_sigma_matsubara(arg):
    DimArg=len(np.shape(arg))
    if(DimArg==5):
        ed2py.get_sigma_matsubara_site(arg)
    elif(DimArg==6):
        ed2py.get_sigma_matsubara_ineq(arg)
    else:
        raise ValueError('Dim(Sigma) /= [{Nlat},Nspin,Nspin,Norb,Norb,:]')
    return arg;

def get_sigma_realaxis(arg):
    DimArg=len(np.shape(arg))
    if(DimArg==5):
        ed2py.get_sigma_realaxis_site(arg)
    elif(DimArg==6):
        ed2py.get_sigma_realaxis_ineq(arg)
    else:
        raise ValueError('Dim(Sigma) /= [{Nlat},Nspin,Nspin,Norb,Norb,:]')
    return arg;

def get_gimp_matsubara(arg):
    DimArg=len(np.shape(arg))
    if(DimArg==5):
        ed2py.get_gimp_matsubara_site(arg)
    elif(DimArg==6):
        ed2py.get_gimp_matsubara_ineq(arg)
    else:
        raise ValueError('Dim(Gimp) /= [{Nlat},Nspin,Nspin,Norb,Norb,:]')
    return arg;

def get_gimp_realaxis(arg):
    DimArg=len(np.shape(arg))
    if(DimArg==5):
        ed2py.get_gimp_realaxis_site(arg)
    elif(DimArg==6):
        ed2py.get_gimp_realaxis_ineq(arg)
    else:
        raise ValueError('Dim(Gimp) /= [{Nlat},Nspin,Nspin,Norb,Norb,:]')
    return arg;

def get_dens(arg):
    DimArg=len(np.shape(arg))
    if(DimArg==1):
        ed2py.get_dens_site(arg)
    elif(DimArg==2):
        ed2py.get_dens_ineq(arg)
    else:
        raise ValueError('Dim(Dens) /= [{Nlat},Norb]')
    return arg;

def get_mag(arg):
    DimArg=len(np.shape(arg))
    if(DimArg==1):
        ed2py.get_mag_site(arg)
    elif(DimArg==2):
        ed2py.get_mag_ineq(arg)
    else:
        raise ValueError('Dim(Mag) /= [{Nlat},Norb]')
    return arg;

def get_docc(arg):
    DimArg=len(np.shape(arg))
    if(DimArg==1):
        ed2py.get_docc_site(arg)
    elif(DimArg==2):
        ed2py.get_docc_ineq(arg)
    else:
        raise ValueError('Dim(Docc) /= [{Nlat},Norb]')
    return arg;

def get_eimp(arg):
    DimArg=len(np.shape(arg))
    if(DimArg==1):
        ed2py.get_eimp_site(arg)
    elif(DimArg==2):
        ed2py.get_eimp_ineq(arg)
    else:
        raise ValueError('Dim(Eimp) /= [{Nlat},4]')
    return arg;

def get_double(arg):
    DimArg=len(np.shape(arg))
    if(DimArg==1):
        ed2py.get_doubles_site(arg)
    elif(DimArg==2):
        ed2py.get_doubles_ineq(arg)
    else:
        raise ValueError('Dim(Doubles) /= [{Nlat},4]')
    return arg;






#FROM ED_BATH_FIT
def chi2_fitgf(func,bath,hloc,ispin=1,iorb=None):
    DimArg=len(np.shape(bath))
    if(DimArg==1):                        #single site(isite,iorb)
        if iorb is None:
            ed2py.chi2_fitgf_site(func,bath,hloc,ispin)
        else:
            ed2py.chi2_fitgf_site(func,bath,hloc,ispin,iorb)        
    elif(DimArg==2):                      #ineq sites(hloc)
        if iorb is None:
            ed2py.chi2_fitgf_ineq(func,bath,hloc,ispin)
        else:
            ed2py.chi2_fitgf_ineq(func,bath,hloc,ispin,iorb)        
    else:
        raise ValueError('Dim(bath) /= [{Nlat},Nb] in Chi2_FitGf')
    return bath;




#FROM ED_MAIN:
def init_solver(bath,hloc):
    DimArg=len(np.shape(bath))
    if(DimArg==1):
        ed2py.init_solver_site(bath,hloc)
    elif(DimArg==2):
        ed2py.init_solver_ineq(bath,hloc)
    else:
        raise ValueError('Dim(Bath) /= [{Nlat},Nb]')
    return ;

def solve(bath,hloc,sflag=True):
    DimArg=len(np.shape(bath))
    if(DimArg==1):
        ed2py.solve_site(bath,hloc,sflag)
    elif(DimArg==2):
        ed2py.solve_ineq(bath,hloc)
    else:
        raise ValueError('Dim(Bath) /= [{Nlat},Nb]')
    return ;






