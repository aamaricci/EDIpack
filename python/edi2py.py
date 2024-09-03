from ctypes import *
import numpy as np
import os,sys
import types

#################################
#link to access global variables
#################################

#dummy class, to be filled
class Link:
    pass

#function that will add a variable to the dummy class, will be called in variable definition
def add_global_variable(obj, dynamic_name, target_object, target_attribute):
    @property
    def getter(self):
        try:
            attrib = getattr(target_object, target_attribute)
            try: #this is for strings
                attrib = attrib.decode()
            except:
                pass
        except: #this is for arrays
            if(len(target_object)>1):
                return [target_object[x] for x in range(len(target_object))]
        return attrib

    @getter.setter
    def setter(self, new_value):
        try: #this is for arrays
            if(len(target_object)>1):
                minlength=min(len(target_object),len(new_value))
                target_object[0:minlength]=new_value[0:minlength]
        except:
            try:
                new_value = new_value.encode()
            except:
                pass
            setattr(target_object, target_attribute, new_value)

    # Dynamically add the property to the class
    setattr(obj.__class__, dynamic_name, getter)
    setattr(obj.__class__, dynamic_name, setter)


######################################
# Load shared library with C-bindings
######################################

system = sys.platform
libext = '.so' 
if(system=='darwin'):
    libext = '.dylib'
libpath = os.path.dirname(os.path.realpath(__file__))
libfile = os.path.join(libpath, 'libedi2py'+libext)
libedi2py = CDLL(libfile)


######################################
# READ_INPUT
######################################

read_input_wrap = libedi2py.read_input
read_input_wrap.argtypes = [c_char_p]  # allows for automatic type conversion
read_input_wrap.restype = None

def read_input(self,input_string):
    c_string = c_char_p(input_string.encode())
    read_input_wrap(c_string)
    
######################################
# ED_MAIN FUNCTIONS
######################################

# Define the function signature for the Fortran function `init_solver_site`.
init_solver_site = libedi2py.init_solver_site
init_solver_site.argtypes = [np.ctypeslib.ndpointer(dtype=float,ndim=1, flags='F_CONTIGUOUS'),c_int]  # allows for automatic type conversion
init_solver_site.restype = None

# Define the function signature for the Fortran function `init_solver_ineq`.
init_solver_ineq = libedi2py.init_solver_ineq
init_solver_ineq.argtypes = [np.ctypeslib.ndpointer(dtype=float,ndim=2, flags='F_CONTIGUOUS'),c_int,c_int]  # allows for automatic type conversion
init_solver_ineq.restype = None

# Define the function signature for the Fortran function `solve_site`.
solve_site = libedi2py.solve_site
solve_site.argtypes = [np.ctypeslib.ndpointer(dtype=float,ndim=1, flags='F_CONTIGUOUS'),c_int,np.ctypeslib.ndpointer(ndim=4, flags='F_CONTIGUOUS'),c_int,c_int,c_int,c_int,c_int]  # allows for automatic type conversion
solve_site.restype = None

# Define the function signature for the Fortran function `solve_ineq`.
solve_ineq = libedi2py.solve_ineq
solve_ineq.argtypes = [np.ctypeslib.ndpointer(dtype=float,ndim=2, flags='F_CONTIGUOUS'),c_int,c_int,np.ctypeslib.ndpointer(ndim=5, flags='F_CONTIGUOUS'),c_int,c_int,c_int,c_int,c_int,c_int]  # allows for automatic type conversion
solve_ineq.restype = None

def init_solver(self,bath):
    dim_bath=np.shape(bath)
    if len(dim_bath)<2:
        init_solver_site(bath,dim_bath[0])
    else:
        init_solver_ineq(bath,dim_bath[0],dim_bath[1])

    
def solve(self,bath,hloc,sflag=True,mpi_lanc=False):
    dim_bath=np.shape(bath)
    dim_hloc=np.shape(hloc)
    if len(dim_bath)<2:
        solve_site(bath,dim_bath[0],hloc,dim_hloc[0],dim_hloc[1],dim_hloc[2],dim_hloc[3],sflag)
    else:
        solve_ineq(bath,dim_bath[0],dim_bath[1],hloc,dim_hloc[0],dim_hloc[1],dim_hloc[2],dim_hloc[3],dim_hloc[4],mpi_lanc)
        

######################################
# ED_BATH FUNCTIONS
######################################

#get_bath_dimension
get_bath_dimension_wrap = libedi2py.get_bath_dimension
get_bath_dimension_wrap.argtypes = None  # allows for automatic type conversion
get_bath_dimension_wrap.restype = c_int

def get_bath_dimension(self):
    return get_bath_dimension_wrap()

#init_hreplica
init_hreplica_direct_nn = libedi2py.init_Hreplica_direct_nn
init_hreplica_direct_nn.argtypes =[np.ctypeslib.ndpointer(dtype=float,ndim=4, flags='F_CONTIGUOUS')]
init_hreplica_direct_nn.restype = None

init_hreplica_direct_so = libedi2py.init_Hreplica_direct_so
init_hreplica_direct_so.argtypes =[np.ctypeslib.ndpointer(dtype=float,ndim=2, flags='F_CONTIGUOUS')]
init_hreplica_direct_so.restype = None

init_hreplica_symmetries_site = libedi2py.init_Hreplica_symmetries_site
init_hreplica_symmetries_site.argtypes =[np.ctypeslib.ndpointer(dtype=float,ndim=5, flags='F_CONTIGUOUS'),np.ctypeslib.ndpointer(dtype=float,ndim=1, flags='F_CONTIGUOUS'),c_int]
init_hreplica_symmetries_site.restype = None

init_hreplica_symmetries_lattice = libedi2py.init_Hreplica_symmetries_ineq
init_hreplica_symmetries_lattice.argtypes =[np.ctypeslib.ndpointer(dtype=float,ndim=5, flags='F_CONTIGUOUS'),np.ctypeslib.ndpointer(dtype=float,ndim=2, flags='F_CONTIGUOUS'),c_int,c_int]
init_hreplica_symmetries_lattice.restype = None
    
def set_Hreplica(self,hloc=None,hvec=None,lambdavec=None):
    
    aux_norb=c_int.in_dll(libedi2py, "Norb").value
    aux_nspin=c_int.in_dll(libedi2py, "Nspin").value

    if hloc is not None:
        shape_hloc=np.shape(hloc)
        if(len(shape_hloc)==4):                       #hloc[Nspin,Nspin,Norb,Norb]
            if shape_hloc == (aux_nspin,aux_nspin,aux_norb,aux_norb):
                init_hreplica_direct_nn(hloc)
            else:
                raise ValueError("Shape of 4-dimensional Hloc != [Nspin,Nspin,Norb,Norb] in set_Hreplica")
        elif(len(shape_hloc)==2):                           #hloc[Nso,Nso]
            if shape_hloc == (aux_nspin*aux_norb,aux_nspin*aux_norb):
                init_hreplica_direct_so(hloc)
            else:
                raise ValueError('Shape of 2-dimensional Hloc != [Nso,Nso] in set_Hreplica')
    elif (hvec is not None) and (lambdavec is not None):
        Dimhvec = np.shape(hvec)
        if(len(Dimhvec) != 5):
            raise ValueError("Dim(hvec) != [Nspin,Nspin,Norb,Norb,Nsym] in set_Hreplica")
        DimLam=np.shape(lambdavec)
        if(len(DimLam)==1):
            init_hreplica_symmetries_site(hvec,lambdavec,DimLam[0])
        elif(len(DimLam)==2):
            init_hreplica_symmetries_lattice(hvec,lambdavec,DimLam[0],DimLam[1])
        else:
             raise ValueError('Shape(lambdavec) != 1 or 2 [Nsym] or [Nineq,Nsym] in set_Hreplica')
    else:
        raise ValueError("set_Hreplica requires either hloc or (hvec,lambdavec) as arguments")
    return ;
    
#break_symmetry_bath
break_symmetry_bath_site = libedi2py.break_symmetry_bath_site
break_symmetry_bath_site.argtypes = [np.ctypeslib.ndpointer(dtype=float,ndim=1, flags='F_CONTIGUOUS'),c_int,c_double,c_double,c_int]
break_symmetry_bath_site.restype = None

break_symmetry_bath_ineq = libedi2py.break_symmetry_bath_ineq
break_symmetry_bath_ineq.argtypes = [np.ctypeslib.ndpointer(dtype=float,ndim=2, flags='F_CONTIGUOUS'),c_int,c_int,c_double,c_double,c_int]
break_symmetry_bath_ineq.restype = None

def break_symmetry_bath(self, bath, field, sign, save=True):
    if save:
        save_int=1
    else:
        save_int=0
    bath_shape = np.shape(bath)
    if (len(bath_shape)) == 1:
        break_symmetry_bath_site(bath,bath_shape[0],field,float(sign),save_int)
    else:
        break_symmetry_bath_ineq(bath,bath_shape[0],bath_shape[1],field,float(sign),save_int)
    return bath
    
#spin_symmetrize_bath
spin_symmetrize_bath_site = libedi2py.spin_symmetrize_bath_site
spin_symmetrize_bath_site.argtypes = [np.ctypeslib.ndpointer(dtype=float,ndim=1, flags='F_CONTIGUOUS'),c_int,c_int]
spin_symmetrize_bath_site.restypes = None

spin_symmetrize_bath_ineq = libedi2py.spin_symmetrize_bath_ineq
spin_symmetrize_bath_ineq.argtypes = [np.ctypeslib.ndpointer(dtype=float,ndim=1, flags='F_CONTIGUOUS'),c_int,c_int,c_int]
spin_symmetrize_bath_ineq.restypes = None
    
def spin_symmetrize_bath(self, bath, save=True):
    if save:
        save_int=1
    else:
        save_int=0
    bath_shape = np.shape(bath)
    if (len(bath_shape)) == 1:
        spin_symmetrize_bath_site(bath,bath_shape[0],save_int)
    else:
        spin_symmetrize_bath_ineq(bath,bath_shape[0],bath_shape[1],save_int)
    return bath
    
#orb_symmetrize_bath
orb_symmetrize_bath_site = libedi2py.orb_symmetrize_bath_site
orb_symmetrize_bath_site.argtypes = [np.ctypeslib.ndpointer(dtype=float,ndim=1, flags='F_CONTIGUOUS'),c_int,c_int]
orb_symmetrize_bath_site.restypes = None

orb_symmetrize_bath_ineq = libedi2py.orb_symmetrize_bath_ineq
orb_symmetrize_bath_ineq.argtypes = [np.ctypeslib.ndpointer(dtype=float,ndim=1, flags='F_CONTIGUOUS'),c_int,c_int,c_int]
orb_symmetrize_bath_ineq.restypes = None
    
def orb_symmetrize_bath(self, bath, save):
    if save:
        save_int=1
    else:
        save_int=0
    bath_shape = np.shape(bath)
    if (len(bath_shape)) == 1:
        orb_symmetrize_bath_site(bath,bath_shape[0],save_int)
    else:
        orb_symmetrize_bath_ineq(bath,bath_shape[0],bath_shape[1],save_int)
    return bath
    
#orb_equality_bath
orb_equality_bath_site = libedi2py.orb_equality_bath_site
orb_equality_bath_site.argtypes = [np.ctypeslib.ndpointer(dtype=float,ndim=1, flags='F_CONTIGUOUS'),c_int,c_int,c_int]
orb_equality_bath_site.restypes = None

orb_equality_bath_ineq = libedi2py.orb_equality_bath_ineq
orb_equality_bath_ineq.argtypes = [np.ctypeslib.ndpointer(dtype=float,ndim=1, flags='F_CONTIGUOUS'),c_int,c_int,c_int,c_int]
orb_equality_bath_ineq.restypes = None
    
def orb_equality_bath(self, bath, indx, save=True):
    aux_norb=c_int.in_dll(libedi2py, "Norb").value
    if save:
        save_int=1
    else:
        save_int=0
    bath_shape = np.shape(bath)
    if (indx < 0) or (indx >= aux_norb):
        raise ValueError("orb_equality_bath: orbital index should be in [0,Norb]")
    else:
        indx = indx + 1 #python to fortran convention 
        if (len(bath_shape)) == 1:
            orb_equality_bath_site(bath,bath_shape[0],indx,save_int)
        else:
            orb_equality_bath_ineq(bath,bath_shape[0],bath_shape[1],indx,save_int)
    return bath
    
    
#ph_symmetrize_bath
ph_symmetrize_bath_site = libedi2py.ph_symmetrize_bath_site
ph_symmetrize_bath_site.argtypes = [np.ctypeslib.ndpointer(dtype=float,ndim=1, flags='F_CONTIGUOUS'),c_int,c_int]
ph_symmetrize_bath_site.restypes = None

ph_symmetrize_bath_ineq = libedi2py.ph_symmetrize_bath_ineq
ph_symmetrize_bath_ineq.argtypes = [np.ctypeslib.ndpointer(dtype=float,ndim=1, flags='F_CONTIGUOUS'),c_int,c_int,c_int]
ph_symmetrize_bath_ineq.restypes = None

def ph_symmetrize_bath(self, bath, save):
    if save:
        save_int=1
    else:
        save_int=0
    bath_shape = np.shape(bath)
    if (len(bath_shape)) == 1:
        ph_symmetrize_bath_site(bath,bath_shape[0],save_int)
    else:
        ph_symmetrize_bath_ineq(bath,bath_shape[0],bath_shape[1],save_int)
    return bath

#ph_trans_bath
ph_trans_bath_site = libedi2py.ph_trans_bath_site
ph_trans_bath_site.argtypes = [np.ctypeslib.ndpointer(dtype=float,ndim=1, flags='F_CONTIGUOUS'),c_int,c_int]
ph_trans_bath_site.restypes = None

ph_trans_bath_ineq = libedi2py.ph_trans_bath_ineq
ph_trans_bath_ineq.argtypes = [np.ctypeslib.ndpointer(dtype=float,ndim=1, flags='F_CONTIGUOUS'),c_int,c_int,c_int]
ph_trans_bath_ineq.restypes = None

def ph_trans_bath(self, bath, save):
    if save:
        save_int=1
    else:
        save_int=0
    bath_shape = np.shape(bath)
    if (len(bath_shape)) == 1:
        ph_trans_bath_site(bath,bath_shape[0],save_int)
    else:
        ph_trans_bath_ineq(bath,bath_shape[0],bath_shape[1],save_int)
    return bath
    
#get_bath_component_dimension
get_bath_component_dimension_wrap = libedi2py.get_bath_component_dimension
get_bath_component_dimension_wrap.argtypes = [c_char_p,np.ctypeslib.ndpointer(dtype=float,ndim=1, flags='F_CONTIGUOUS')]
get_bath_component_dimension_wrap.restypes = None

def get_bath_component_dimension(self,typ="e"):
    ndim=np.zeros(3,order="F")
    get_bath_component_dimension_wrap(c_char_p(typ.encode()),ndim)
    ndim = [int(a) for a in ndim]
    return ndim

#get_bath_component
get_bath_component_wrap = libedi2py.get_bath_component
get_bath_component_wrap.argtypes = [np.ctypeslib.ndpointer(dtype=float,ndim=1, flags='F_CONTIGUOUS'),c_int,c_int,c_int,np.ctypeslib.ndpointer(dtype=float,ndim=1, flags='F_CONTIGUOUS'),c_int,c_char_p]
get_bath_component_wrap.restypes = None

def get_bath_component(self,array,bath,typ="e"):
    DimArray=np.shape(array)
    DimBath=np.shape(bath)
    if not len(DimArray)==3:
        raise ValueError('Shape(array) != 3 in get_bath_component')
    if not len(DimBath)==1:
        raise ValueError('Shape(bath) != 1 in get_bath_component')
    edi2py.get_bath_component(array,DimArray[0],DimArray[1],DimArray[2],bath,DimBath[0],c_char_p(typ.encode()))
    return array

#set_bath_component
set_bath_component_wrap = libedi2py.set_bath_component
set_bath_component_wrap.argtypes = [np.ctypeslib.ndpointer(dtype=float,ndim=1, flags='F_CONTIGUOUS'),c_int,c_int,c_int,np.ctypeslib.ndpointer(dtype=float,ndim=1, flags='F_CONTIGUOUS'),c_int,c_char_p]
set_bath_component_wrap.restypes = None

def set_bath_component(self,array,bath,typ="e"):
    DimArray=np.shape(array)
    DimBath=np.shape(bath)
    if not len(DimArray)==3:
        raise ValueError('Shape(array) != 3 in get_bath_component')
    if not len(DimBath)==1:
        raise ValueError('Shape(bath) != 1 in get_bath_component')
    edi2py.set_bath_component(array,DimArray[0],DimArray[1],DimArray[2],bath,DimBath[0],c_char_p(typ.encode()))
    return bath
 
#copy_bath_component
copy_bath_component_wrap = libedi2py.copy_bath_component
copy_bath_component_wrap.argtypes = [np.ctypeslib.ndpointer(dtype=float,ndim=1, flags='F_CONTIGUOUS'),c_int,np.ctypeslib.ndpointer(dtype=float,ndim=1, flags='F_CONTIGUOUS'),c_int,c_char_p]
copy_bath_component_wrap.restypes = None

def copy_bath_component(bathIN,bathOUT,typ="e"):
    DimBathIn=p.shape(bathIN)
    DimBathOut=np.shape(bathOUT)
    if not DimBathIn==1:
        raise ValueError('Shape(bathIN) != 1 in copy_bath_component')
    if not DimBathOut==1:
        raise ValueError('Shape(bathOUT) != 1 in copy_bath_component')
    edi2py.copy_bath_component(bathIN,bathOUT,c_char_p(typ.encode()))
    return bathOUT


######################################
# ED_IO FUNCTIONS
######################################

#sigma_matsubara
get_sigma_matsubara_site = libedi2py.get_sigma_matsubara_site
get_sigma_matsubara_site.argtypes = [np.ctypeslib.ndpointer(dtype=complex,ndim=5, flags='F_CONTIGUOUS'),c_int,c_int,c_int,c_int,c_int]  # allows for automatic type conversion
get_sigma_matsubara_site.restype = None

get_sigma_matsubara_ineq = libedi2py.get_sigma_matsubara_ineq
get_sigma_matsubara_ineq.argtypes = [np.ctypeslib.ndpointer(dtype=complex,ndim=6, flags='F_CONTIGUOUS'),c_int,c_int,c_int,c_int,c_int,c_int]  # allows for automatic type conversion
get_sigma_matsubara_ineq.restype = None

def get_sigma_matsubara(self,Smats):
    dim_Smats=np.shape(Smats)
    if len(dim_Smats)==5:
        get_sigma_matsubara_site(Smats,dim_Smats[0],dim_Smats[1],dim_Smats[2],dim_Smats[3],dim_Smats[4])
    else:
        get_sigma_matsubara_ineq(Smats,dim_Smats[0],dim_Smats[1],dim_Smats[2],dim_Smats[3],dim_Smats[4],dim_Smats[5])
    return Smats

#sigma_realaxis
get_sigma_realaxis_site = libedi2py.get_sigma_realaxis_site
get_sigma_realaxis_site.argtypes = [np.ctypeslib.ndpointer(dtype=complex,ndim=5, flags='F_CONTIGUOUS'),c_int,c_int,c_int,c_int,c_int]  # allows for automatic type conversion
get_sigma_realaxis_site.restype = None

get_sigma_realaxis_ineq = libedi2py.get_sigma_matsubara_ineq
get_sigma_realaxis_ineq.argtypes = [np.ctypeslib.ndpointer(dtype=complex,ndim=6, flags='F_CONTIGUOUS'),c_int,c_int,c_int,c_int,c_int,c_int]  # allows for automatic type conversion
get_sigma_realaxis_ineq.restype = None

        
def get_sigma_realaxis(self,Sreal):
    dim_Sreal=np.shape(Sreal)
    if len(dim_Sreal)==5:
        get_sigma_realaxis_site(Sreal,dim_Sreal[0],dim_Sreal[1],dim_Sreal[2],dim_Sreal[3],dim_Sreal[4])
    else:
        get_sigma_realaxis_ineq(Sreal,dim_Sreal[0],dim_Sreal[1],dim_Sreal[2],dim_Sreal[3],dim_Sreal[4],dim_Sreal[5])
    return Sreal
        
#gimp_matsubara
get_gimp_matsubara_site = libedi2py.get_gimp_matsubara_site
get_gimp_matsubara_site.argtypes = [np.ctypeslib.ndpointer(dtype=complex,ndim=5, flags='F_CONTIGUOUS'),c_int,c_int,c_int,c_int,c_int]  # allows for automatic type conversion
get_gimp_matsubara_site.restype = None

get_gimp_matsubara_ineq = libedi2py.get_gimp_matsubara_ineq
get_gimp_matsubara_ineq.argtypes = [np.ctypeslib.ndpointer(dtype=complex,ndim=6, flags='F_CONTIGUOUS'),c_int,c_int,c_int,c_int,c_int,c_int]  # allows for automatic type conversion
get_gimp_matsubara_ineq.restype = None

def get_gimp_matsubara(self,Gmats):
    dim_Gmats=np.shape(Gmats)
    if len(dim_Gmats)==5:
        get_gimp_matsubara_site(Gmats,dim_Gmats[0],dim_Gmats[1],dim_Gmats[2],dim_Gmats[3],dim_Gmats[4])
    else:
        get_gimp_matsubara_ineq(Gmats,dim_Gmats[0],dim_Gmats[1],dim_Gmats[2],dim_Gmats[3],dim_Gmats[4],dim_Gmats[5])
    return Gmats

#gimp_realaxis
get_gimp_realaxis_site = libedi2py.get_gimp_realaxis_site
get_gimp_realaxis_site.argtypes = [np.ctypeslib.ndpointer(dtype=complex,ndim=5, flags='F_CONTIGUOUS'),c_int,c_int,c_int,c_int,c_int]  # allows for automatic type conversion
get_gimp_realaxis_site.restype = None

get_gimp_realaxis_ineq = libedi2py.get_gimp_matsubara_ineq
get_gimp_realaxis_ineq.argtypes = [np.ctypeslib.ndpointer(dtype=complex,ndim=6, flags='F_CONTIGUOUS'),c_int,c_int,c_int,c_int,c_int,c_int]  # allows for automatic type conversion
get_gimp_realaxis_ineq.restype = None

        
def get_gimp_realaxis(self,Greal):
    dim_Greal=np.shape(Greal)
    if len(dim_Greal)==5:
        get_gimp_realaxis_site(Greal,dim_Greal[0],dim_Greal[1],dim_Greal[2],dim_Greal[3],dim_Greal[4])
    else:
        get_gimp_realaxis_ineq(Greal,dim_Greal[0],dim_Greal[1],dim_Greal[2],dim_Greal[3],dim_Greal[4],dim_Greal[5])
    return Greal
        
        
#dens
get_dens_site = libedi2py.get_dens_site
get_dens_site.argtypes = [np.ctypeslib.ndpointer(dtype=float,ndim=1, flags='F_CONTIGUOUS'),c_int]  # allows for automatic type conversion
get_dens_site.restype = None

get_dens_ineq = libedi2py.get_dens_ineq
get_dens_ineq.argtypes = [np.ctypeslib.ndpointer(dtype=float,ndim=2, flags='F_CONTIGUOUS'),c_int,c_int,c_int]  # allows for automatic type conversion
get_dens_ineq.restype = None
        
def get_dens(self,Nlat=-1,iorb=-1):
    aux_norb=c_int.in_dll(libedi2py, "Norb").value
    if (Nlat < 0):
        dens = np.zeros((aux_norb),dtype='float',order='F')
        get_dens_site(dens,aux_norb)
        if(iorb > -1):
            return dens[iorb]
        else:
            return dens
    else:
        dens = np.zeros((Nlat,aux_norb),dtype='float',order='F')
        get_dens_ineq(dens,Nlat,aux_norb,Nlat)
        if(iorb > -1):
            return dens[iorb]
        else:
            return dens
        
        
#mag
get_mag_site = libedi2py.get_mag_site
get_mag_site.argtypes = [np.ctypeslib.ndpointer(dtype=float,ndim=1, flags='F_CONTIGUOUS'),c_int]  # allows for automatic type conversion
get_mag_site.restype = None

get_mag_ineq = libedi2py.get_mag_ineq
get_mag_ineq.argtypes = [np.ctypeslib.ndpointer(dtype=float,ndim=2, flags='F_CONTIGUOUS'),c_int,c_int,c_int]  # allows for automatic type conversion
get_mag_ineq.restype = None
        
def get_mag(self,Nlat=-1,iorb=-1):
    aux_norb=c_int.in_dll(libedi2py, "Norb").value
    if (Nlat < 0):
        mag = np.zeros((aux_norb),dtype='float',order='F')
        get_mag_site(mag,aux_norb)
        if(iorb > -1):
            return mag[iorb]
        else:
            return mag
    else:
        mag = np.zeros((Nlat,aux_norb),dtype='float',order='F')
        get_mag_ineq(mag,Nlat,aux_norb,Nlat)
        if(iorb > -1):
            return mag[iorb]
        else:
            return mag

#docc
get_docc_site = libedi2py.get_docc_site
get_docc_site.argtypes = [np.ctypeslib.ndpointer(dtype=float,ndim=1, flags='F_CONTIGUOUS'),c_int]  # allows for automatic type conversion
get_docc_site.restype = None

get_docc_ineq = libedi2py.get_docc_ineq
get_docc_ineq.argtypes = [np.ctypeslib.ndpointer(dtype=float,ndim=2, flags='F_CONTIGUOUS'),c_int,c_int,c_int]  # allows for automatic type conversion
get_docc_ineq.restype = None
        
def get_docc(self,Nlat=-1,iorb=-1):
    aux_norb=c_int.in_dll(libedi2py, "Norb").value
    if (Nlat < 0):
        docc = np.zeros((aux_norb),dtype='float',order='F')
        get_docc_site(docc,aux_norb)
        if(iorb > -1):
            return docc[iorb]
        else:
            return docc
    else:
        docc = np.zeros((Nlat,aux_norb),dtype='float',order='F')
        get_docc_ineq(docc,Nlat,aux_norb,Nlat)
        if(iorb > -1):
            return docc[iorb]
        else:
            return docc
        
#eimp
get_eimp_site = libedi2py.get_eimp_site
get_eimp_site.argtypes = [np.ctypeslib.ndpointer(dtype=float,ndim=1, flags='F_CONTIGUOUS')]  # allows for automatic type conversion
get_eimp_site.restype = None

get_eimp_ineq = libedi2py.get_eimp_ineq
get_eimp_ineq.argtypes = [np.ctypeslib.ndpointer(dtype=float,ndim=2, flags='F_CONTIGUOUS'),c_int]  # allows for automatic type conversion
get_eimp_ineq.restype = None
        
def get_eimp(self,Nlat=-1):
    if (Nlat < 0):
        eimp = np.zeros(4,dtype='float',order='F')
        get_eimp_site(eimp)
        return eimp
    else:
        eimp = np.zeros((Nlat,4),dtype='float',order='F')
        get_eimp_ineq(eimp,Nlat)
        return eimp
        
#doubles
get_doubles_site = libedi2py.get_doubles_site
get_doubles_site.argtypes = [np.ctypeslib.ndpointer(dtype=float,ndim=1, flags='F_CONTIGUOUS')]  # allows for automatic type conversion
get_doubles_site.restype = None

get_doubles_ineq = libedi2py.get_doubles_ineq
get_doubles_ineq.argtypes = [np.ctypeslib.ndpointer(dtype=float,ndim=2, flags='F_CONTIGUOUS'),c_int]  # allows for automatic type conversion
get_doubles_ineq.restype = None
        
def get_doubles(self,Nlat=-1):
    if (Nlat < 0):
        doubles = np.zeros(4,dtype='float',order='F')
        get_doubles_site(doubles)
        return doubles
    else:
        doubles = np.zeros((Nlat,4),dtype='float',order='F')
        get_doubles_ineq(doubles,Nlat)
        return doubles

######################################
# ED_BATH_FIT FUNCTIONS
######################################

chi2_fitgf_site = libedi2py.chi2_fitgf_site
chi2_fitgf_site.argtypes = [np.ctypeslib.ndpointer(dtype=complex,ndim=5, flags='F_CONTIGUOUS'),c_int,c_int,c_int,c_int,c_int,np.ctypeslib.ndpointer(dtype=float,ndim=1, flags='F_CONTIGUOUS'),c_int,c_int,c_int] 
chi2_fitgf_site.restype = None

chi2_fitgf_ineq = libedi2py.chi2_fitgf_ineq
chi2_fitgf_ineq.argtypes = [np.ctypeslib.ndpointer(dtype=complex,ndim=6, flags='F_CONTIGUOUS'),c_int,c_int,c_int,c_int,c_int,c_int,
                            np.ctypeslib.ndpointer(dtype=float,ndim=2, flags='F_CONTIGUOUS'),c_int,c_int,
                            np.ctypeslib.ndpointer(dtype=float,ndim=5, flags='F_CONTIGUOUS'),c_int,c_int,c_int,c_int,c_int,
                            c_int] 
chi2_fitgf_ineq.restype = None

def chi2_fitgf(self,func,bath,hloc=None,ispin=1,iorb=0):
    dim_func=np.shape(func)
    dim_bath=np.shape(bath)
    dim_hloc=np.shape(hloc)
    if(len(dim_func)==5):
        chi2_fitgf_site(func,dim_func[0],dim_func[1],dim_func[2],dim_func[3],dim_func[4],bath,dim_bath[0],ispin,iorb)
    else:
        chi2_fitgf_ineq(func,dim_func[0],dim_func[1],dim_func[2],dim_func[3],dim_func[4],dim_func[5],bath,dim_bath[0],dim_bath[1],hloc,dim_hloc[0],dim_hloc[1],dim_hloc[2],dim_hloc[3],dim_hloc[4],ispin)
        


######################################
# ED_AUX_FUNX FUNCTIONS
######################################

#search_variable
search_variable_wrap = libedi2py.search_variable
search_variable_wrap.argtypes = [np.ctypeslib.ndpointer(dtype=float,ndim=1, flags='F_CONTIGUOUS'),np.ctypeslib.ndpointer(dtype=float,ndim=1, flags='F_CONTIGUOUS'),np.ctypeslib.ndpointer(dtype=int,ndim=1, flags='F_CONTIGUOUS')]
search_variable_wrap.restype = None

def search_variable(self,var,ntmp,converged):
    var = np.asarray([var])
    ntmp = np.asarray([ntmp])
    converged = np.asarray([converged])
    conv_int=int(converged)
    search_variable_wrap(var,ntmp,converged)
    if conv_int[0]==0:
        converged=False
    else:
        converged=True
    return err[0],conv_bool

#check_convergence
check_convergence_wrap = libedi2py.check_convergence
check_convergence_wrap.argtypes = [np.ctypeslib.ndpointer(dtype=complex,ndim=1, flags='F_CONTIGUOUS'),c_int,c_double,c_int,c_int,np.ctypeslib.ndpointer(dtype=float,ndim=1, flags='F_CONTIGUOUS'),np.ctypeslib.ndpointer(dtype=int,ndim=1, flags='F_CONTIGUOUS')]
check_convergence_wrap.restype = None

def check_convergence(self,func,threshold,N1,N2):
    err=np.asarray([1.0])
    converged=np.asarray([0])
    dim_func=np.shape(func)
    check_convergence_wrap(func,dim_func[0],threshold,N1,N2,err,converged)
    if converged[0]==0:
        conv_bool=False
    else:
        conv_bool=True
    return err[0],conv_bool

####################################################################
# Create the global_env class (this is what the python module sees)
####################################################################

global_env=Link()

#variables
add_global_variable(global_env, "Nbath", c_int.in_dll(libedi2py, "Nbath"), "value")
add_global_variable(global_env, "Norb", c_int.in_dll(libedi2py, "Norb"), "value")
add_global_variable(global_env, "Nspin", c_int.in_dll(libedi2py, "Nspin"), "value")
add_global_variable(global_env, "Nloop", c_int.in_dll(libedi2py, "Nloop"), "value")
add_global_variable(global_env, "Nph", c_int.in_dll(libedi2py, "Nph"), "value")
add_global_variable(global_env, "Nsuccess", c_int.in_dll(libedi2py, "Nsuccess"), "value")
add_global_variable(global_env, "Lmats", c_int.in_dll(libedi2py, "Lmats"), "value")
add_global_variable(global_env, "Lreal", c_int.in_dll(libedi2py, "Lreal"), "value")
add_global_variable(global_env, "Ltau", c_int.in_dll(libedi2py, "Ltau"), "value")
add_global_variable(global_env, "Lpos", c_int.in_dll(libedi2py, "Lpos"), "value")
add_global_variable(global_env, "LOGfile", c_int.in_dll(libedi2py, "LOGfile"), "value")

add_global_variable(global_env, "Uloc", ARRAY(c_double, 5).in_dll(libedi2py, "Uloc"), "value")
add_global_variable(global_env, "Ust", c_double.in_dll(libedi2py, "Ust"), "value")
add_global_variable(global_env, "Jh", c_double.in_dll(libedi2py, "Jh"), "value")
add_global_variable(global_env, "Jx", c_double.in_dll(libedi2py, "Jx"), "value")
add_global_variable(global_env, "Jp", c_double.in_dll(libedi2py, "Jp"), "value")
add_global_variable(global_env, "xmu", c_double.in_dll(libedi2py, "xmu"), "value")
add_global_variable(global_env, "beta", c_double.in_dll(libedi2py, "beta"), "value")
add_global_variable(global_env, "dmft_error", c_double.in_dll(libedi2py, "dmft_error"), "value")
add_global_variable(global_env, "eps", c_double.in_dll(libedi2py, "eps"), "value")
add_global_variable(global_env, "wini", c_double.in_dll(libedi2py, "wini"), "value")
add_global_variable(global_env, "wfin", c_double.in_dll(libedi2py, "wfin"), "value")
add_global_variable(global_env, "xmin", c_double.in_dll(libedi2py, "xmin"), "value")
add_global_variable(global_env, "xmax", c_double.in_dll(libedi2py, "xmax"), "value")
add_global_variable(global_env, "sb_field", c_double.in_dll(libedi2py, "sb_field"), "value")
add_global_variable(global_env, "nread", c_double.in_dll(libedi2py, "nread"), "value")


#functions
global_env.read_input = types.MethodType(read_input, global_env)

global_env.init_solver = types.MethodType(init_solver, global_env)
global_env.solve = types.MethodType(solve, global_env)

global_env.get_bath_dimension = types.MethodType(get_bath_dimension, global_env)
global_env.set_Hreplica = types.MethodType(set_Hreplica, global_env)
global_env.break_symmetry_bath = types.MethodType(break_symmetry_bath, global_env)
global_env.spin_symmetrize_bath = types.MethodType(spin_symmetrize_bath, global_env)
global_env.orb_symmetrize_bath = types.MethodType(orb_symmetrize_bath, global_env)
global_env.orb_equality_bath = types.MethodType(orb_equality_bath, global_env)
global_env.ph_symmetrize_bath = types.MethodType(ph_symmetrize_bath, global_env)
global_env.ph_trans_bath = types.MethodType(ph_trans_bath, global_env)
global_env.get_bath_component_dimension = types.MethodType(get_bath_component_dimension, global_env)
global_env.get_bath_component = types.MethodType(get_bath_component, global_env)
global_env.set_bath_component = types.MethodType(set_bath_component, global_env)
global_env.copy_bath_component = types.MethodType(copy_bath_component, global_env)

global_env.get_sigma_matsubara = types.MethodType(get_sigma_matsubara, global_env)
global_env.get_sigma_realaxis = types.MethodType(get_sigma_realaxis, global_env)
global_env.get_gimp_matsubara = types.MethodType(get_gimp_matsubara, global_env)
global_env.get_gimp_realaxis = types.MethodType(get_gimp_realaxis, global_env)
global_env.get_dens = types.MethodType(get_dens, global_env)
global_env.get_mag = types.MethodType(get_mag, global_env)
global_env.get_docc = types.MethodType(get_docc, global_env)
global_env.get_eimp = types.MethodType(get_eimp, global_env)
global_env.get_doubles = types.MethodType(get_doubles, global_env)

global_env.chi2_fitgf = types.MethodType(chi2_fitgf, global_env)

global_env.search_variable = types.MethodType(search_variable, global_env)
global_env.check_convergence = types.MethodType(check_convergence, global_env)
