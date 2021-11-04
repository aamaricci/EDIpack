import numpy as np
import edi2py                         #import the Fortran module
import edipy                          #import the PYTHON library itself
from edi2py import *
from edipy import *


global Norb         #Norb =# of impurity orbitals
global Nspin        #Nspin=# spin degeneracy (max 2)
global Nloop        #max dmft loop variables
global Nph
global Uloc,Ust
global Jh,Jx,Jp
global xmu          #chemical potential
global beta         #inverse temperature
global g_ph
global w0_ph        #w0_ph: phonon frequency (constant)
global eps          #broadening
global wini,wfin    #frequency range
global xmin,xmax    #x-range for the local phonons probability distribution function
global Nsuccess     #Number of repeated success to fall below convergence threshold  
global dmft_error   #dmft convergence threshold
global sb_field     #symmetry breaking field  
global cg_scheme    #fit scheme: delta (default), weiss for G0
global nread        #fixed density. if 0.d0 fixed chemical potential calculation.
global Lmats        #Number of Matsuabra Freq.
global Lreal        #Number of real-axis freq.
global Lpos         #
global Hfile
global LOGfile
global bath_type    #flag to set bath type: normal (1bath/imp), hybrid(1bath)



#FROM ED_INPUT_VARS
Norb       = edi2py.ed_input_vars.norb
Nspin      = edi2py.ed_input_vars.nspin
Nloop      = edi2py.ed_input_vars.nloop
Nph        = edi2py.ed_input_vars.nph
Uloc       = edi2py.ed_input_vars.uloc
Ust        = edi2py.ed_input_vars.ust
Jh         = edi2py.ed_input_vars.jh
Jx         = edi2py.ed_input_vars.jx
Jp         = edi2py.ed_input_vars.jp
xmu        = edi2py.ed_input_vars.xmu
beta       = edi2py.ed_input_vars.beta
g_ph       = edi2py.ed_input_vars.g_ph
w0_ph      = edi2py.ed_input_vars.w0_ph
eps        = edi2py.ed_input_vars.eps
wini       = edi2py.ed_input_vars.wini
wfin       = edi2py.ed_input_vars.wfin
xmin       = edi2py.ed_input_vars.xmin
xmax       = edi2py.ed_input_vars.xmax
Nsuccess   = edi2py.ed_input_vars.nsuccess
dmft_error = edi2py.ed_input_vars.dmft_error
sb_field   = edi2py.ed_input_vars.sb_field
cg_scheme  = edi2py.ed_input_vars.cg_scheme
nread      = edi2py.ed_input_vars.nread
Lmats      = edi2py.ed_input_vars.lmats
Lreal      = edi2py.ed_input_vars.lreal
Lpos       = edi2py.ed_input_vars.lpos
Hfile      = edi2py.ed_input_vars.hfile
LOGfile    = edi2py.ed_input_vars.logfile
bath_type  = edi2py.ed_input_vars.bath_type




#FROM ED_BATH:
def get_bath_dimension():    
    Nb = edi2py.get_bath_dimension()
    return Nb;

def set_Hreplica(hloc=None,hvec=None,lambdavec=None):
    if hloc is not None:
        Dimhloc=len(np.shape(hloc))
        if(Dimhloc==2):                       #hloc[Nso,Nso]
            edi2py.init_Hreplica_direct_so(hloc)
        elif(Dimhloc==4):                     #hloc[Nspin,Nspin,Norb,Norb]
            edi2py.init_Hreplica_direct_nn(hloc)  
        else:
            raise ValueError('Shape(hloc) != [Nso,Nso] or [Nspin,Nspin,Norb,Norb] in set_Hreplica')
    elif None not in (hvec,lambdavec):
        Dimhvec=len(np.shape(hvec))
        if(Dimhvec is not 5):
            raise ValueError("Dim(hvec) != [Nspin,Nspin,Norb,Norb,Nsym] in set_Hreplica")
        DimLam=len(np.shape(lambdavec))
        if(DimLam==1):
            edi2py.init_Hreplica_symmetries_site(hvec,lambdavec)
        elif(DimLam==2):
            edi2py.init_Hreplica_symmetries_lattice(hvec,lambdavec)
        else:
             raise ValueError('Shape(lambdavec) != 1 or 2 [Nsym] or [Nineq,Nsym] in set_Hreplica')
    else:
        raise ValueError("set_Hreplica requires either hloc or (hvec,lambdavec) as arguments")
    return ;


def spin_symmetrize_bath(bath,save:bool=True):
    assert isinstance(save, bool), 'type(save)!=bool'
    DimBath=len(np.shape(bath))
    if(DimBath==1):
        edi2py.spin_symmetrize_bath_site(bath,save)
    elif(DimBath==2):
        edi2py.spin_symmetrize_bath_ineq(bath,save)
    else:
        raise ValueError('Shape(Bath) != [{Nlat},Nb]')
    return bath;

def orbs_symmetrize_bath(bath:float,save:bool=True):
    assert isinstance(bath, float), 'type(bath)!=float'
    assert isinstance(save, bool), 'type(save)!=bool'
    DimBath=len(np.shape(bath))
    if(DimBath==1):
        edi2py.orbs_symmetrize_bath_site(bath,save)
    elif(DimBath==2):
        edi2py.orbs_symmetrize_bath_ineq(bath,save)
    else:
        raise ValueError('Shape(Bath) != [{Nlat},Nb]')
    return bath;

def orb_equality_bath(bath,save:bool=True):
    assert isinstance(save, bool), 'type(save)!=bool'
    DimBath=len(np.shape(bath))
    if(DimBath==1):
        edi2py.orb_equality_bath_site(bath,save)
    elif(DimBath==2):
        edi2py.orb_equality_bath_ineq(bath,save)
    else:
        raise ValueError('Shape(Bath) != [{Nlat},Nb]')
    return bath;

def ph_symmetrize_bath(bath,save:bool=True):
    assert isinstance(save, bool), 'type(save)!=bool'
    DimBath=len(np.shape(bath))
    if(DimBath==1):
        edi2py.ph_symmetrize_bath_site(bath,save)
    elif(DimBath==2):
        edi2py.ph_symmetrize_bath_ineq(bath,save)
    else:
        raise ValueError('Shape(Bath) != [{Nlat},Nb]')
    return bath;

def ph_trans_bath(bath,save:bool=True):
    assert isinstance(save, bool), 'type(save)!=bool'
    DimBath=len(np.shape(bath))
    if(DimBath==1):
        edi2py.ph_trans_bath_site(bath,save)
    elif(DimBath==2):
        edi2py.ph_trans_bath_ineq(bath,save)
    else:
        raise ValueError('Shape(Bath) != [{Nlat},Nb]')
    return bath;

def break_symmetry_bath(bath,field:float,sign:float,save:bool=True):
    assert isinstance(field, float), 'type(field)!=float'
    assert isinstance(sign, float), 'type(sign)!=float'
    assert isinstance(save, bool), 'type(save)!=bool'
    DimBath=len(np.shape(bath))
    if(DimBath==1):
        bath=edi2py.break_symmetry_bath_site(field,sign,save)
    elif(DimBath==2):
        bath=edi2py.break_symmetry_bath_ineq(field,sign,save)
    else:
        raise ValueError('Shape(Bath) != [{Nlat},Nb]')
    return bath;

def get_bath_component_dimension(type:str,ndim):
    assert isinstance(type, str), 'type(type)!=str'
    DimNdim=len(np.shape(ndim))
    if not DimNdim==1:
        raise ValueError('Shape(ndim) != 1 in get_bath_component_dimension')
    edi2py.get_bath_component_dimension(type,ndim)
    return ndim;

def get_bath_component(array,bath,type:str):
    assert isinstance(type, str), 'type(type)!=str'
    DimArray=len(np.shape(array))
    DimBath=len(np.shape(bath))
    if not DimArray==3:
        raise ValueError('Shape(array) != 3 in get_bath_component')
    if not DimBath==1:
        raise ValueError('Shape(bath) != 1 in get_bath_component')
    edi2py.get_bath_component(array,bath,type)
    return ;

def set_bath_component(array,bath,type:str):
    assert isinstance(type, str), 'type(type)!=str'
    DimArray=len(np.shape(array))
    DimBath=len(np.shape(bath))
    if not DimArray==3:
        raise ValueError('Shape(array) != 3 in set_bath_component')
    if not DimBath==1:
        raise ValueError('Shape(bath) != 1 in set_bath_component')
    edi2py.set_bath_component(array,bath,type)
    return ;

def copy_bath_component(bathIN,bathOUT,type:str):
    assert isinstance(type, str), 'type(type)!=str'
    DimArray=len(np.shape(array))
    DimBathIn=len(np.shape(bathIN))
    DimBathOut=len(np.shape(bathOUT))
    if not DimBathIn==1:
        raise ValueError('Shape(bathIN) != 1 in copy_bath_component')
    if not DimBathOut==1:
        raise ValueError('Shape(bathOUT) != 1 in copy_bath_component')
    edi2py.copy_bath_component(bathIN,bathOUT,type)
    return ;




# #FROM ED_AUX_FUNX:
def search_variable(var:float,ntmp:float,converged:bool):
    assert isinstance(var, float),     'type(var)!=float'
    assert isinstance(ntmp, float),    'type(ntmp)!=float'
    assert isinstance(converged, bool),'type(converged)!=bool'
    edi2py.search_variable(var,ntmp,converged)
    return converged;






#FROM ED_IO:
def get_sigma_matsubara(arg):
    DimArg=len(np.shape(arg))
    if(DimArg==5):
        edi2py.get_sigma_matsubara_site(arg)
    elif(DimArg==6):
        edi2py.get_sigma_matsubara_ineq(arg)
    else:
        raise ValueError('Shape(Sigma) != [{Nlat},Nspin,Nspin,Norb,Norb,:]')
    return arg;

def get_sigma_realaxis(arg):
    DimArg=len(np.shape(arg))
    if(DimArg==5):
        edi2py.get_sigma_realaxis_site(arg)
    elif(DimArg==6):
        edi2py.get_sigma_realaxis_ineq(arg)
    else:
        raise ValueError('Shape(Sigma) != [{Nlat},Nspin,Nspin,Norb,Norb,:]')
    return arg;

def get_gimp_matsubara(arg):
    DimArg=len(np.shape(arg))
    if(DimArg==5):
        edi2py.get_gimp_matsubara_site(arg)
    elif(DimArg==6):
        edi2py.get_gimp_matsubara_ineq(arg)
    else:
        raise ValueError('Shape(Gimp) != [{Nlat},Nspin,Nspin,Norb,Norb,:]')
    return arg;

def get_gimp_realaxis(arg):
    DimArg=len(np.shape(arg))
    if(DimArg==5):
        edi2py.get_gimp_realaxis_site(arg)
    elif(DimArg==6):
        edi2py.get_gimp_realaxis_ineq(arg)
    else:
        raise ValueError('Shape(Gimp) != [{Nlat},Nspin,Nspin,Norb,Norb,:]')
    return arg;

def get_dens(arg):
    DimArg=len(np.shape(arg))
    if(DimArg==1):
        edi2py.get_dens_site(arg)
    elif(DimArg==2):
        edi2py.get_dens_ineq(arg)
    else:
        raise ValueError('Shape(Dens) != [{Nlat},Norb]')
    return arg;

def get_mag(arg):
    DimArg=len(np.shape(arg))
    if(DimArg==1):
        edi2py.get_mag_site(arg)
    elif(DimArg==2):
        edi2py.get_mag_ineq(arg)
    else:
        raise ValueError('Shape(Mag) != [{Nlat},Norb]')
    return arg;

def get_docc(arg):
    DimArg=len(np.shape(arg))
    if(DimArg==1):
        edi2py.get_docc_site(arg)
    elif(DimArg==2):
        edi2py.get_docc_ineq(arg)
    else:
        raise ValueError('Shape(Docc) != [{Nlat},Norb]')
    return arg;

def get_eimp(arg):
    DimArg=len(np.shape(arg))
    if(DimArg==1):
        edi2py.get_eimp_site(arg)
    elif(DimArg==2):
        edi2py.get_eimp_ineq(arg)
    else:
        raise ValueError('Shape(Eimp) != [{Nlat},4]')
    return arg;

def get_double(arg):
    DimArg=len(np.shape(arg))
    if(DimArg==1):
        edi2py.get_doubles_site(arg)
    elif(DimArg==2):
        edi2py.get_doubles_ineq(arg)
    else:
        raise ValueError('Shape(Doubles) != [{Nlat},4]')
    return arg;






#FROM ED_BATH_FIT
def chi2_fitgf(func,bath,hloc=None,ispin:int=1,iorb:int=None):
    assert isinstance(ispin, int),    'type(ispin)!=int'
    assert isinstance(iorb, int),   'type(iorb)!=int'
    DimArg=len(np.shape(bath))
    if(DimArg==1):                        #single site(isite,iorb)
        if iorb is None:
            edi2py.chi2_fitgf_site(func,bath,ispin)
        else:
            edi2py.chi2_fitgf_site(func,bath,ispin,iorb)        
    elif(DimArg==2):                      #ineq sites(hloc)
        if hloc is None:
            raise ValueError('chi2_fitgf: missing required Hloc as argument')
        edi2py.chi2_fitgf_ineq(func,bath,hloc,ispin)
    else:
        raise ValueError('Shape(bath) != [{Nlat},Nb] in Chi2_FitGf')
    return bath;




#FROM ED_MAIN:
def init_solver(bath):
    DimArg=len(np.shape(bath))
    if(DimArg==1):
        edi2py.init_solver_site(bath)
    elif(DimArg==2):
        edi2py.init_solver_ineq(bath)
    else:
        raise ValueError('Shape(Bath) != [{Nlat},Nb]')
    return ;

def solve(bath,hloc,sflag:bool=True,mpi_lanc:bool=True):
    assert isinstance(sflag, bool),   'type(sflag)!=bool'
    assert isinstance(mpi_lanc, bool),'type(mpi_lanc)!=bool'    
    DimArg=len(np.shape(bath))
    if(DimArg==1):
        edi2py.solve_site(bath,hloc,sflag)
    elif(DimArg==2):
        edi2py.solve_ineq(bath,hloc,mpi_lanc)
    else:
        raise ValueError('Shape(Bath) != [{Nlat},Nb]')
    return ;






