MODULE ED_FIT_COMMON
  USE SF_CONSTANTS
  USE SF_OPTIMIZE, only:fmin_cg,fmin_cgplus,fmin_cgminimize
  USE SF_LINALG,   only:eye,zeye,inv,inv_her,operator(.x.)
  USE SF_IOTOOLS,  only:reg,free_unit,txtfy
  USE SF_ARRAYS,   only:arange
  USE SF_MISC,     only:assert_shape
  !
  USE ED_INPUT_VARS
  USE ED_VARS_GLOBAL  
  USE ED_AUX_FUNX
  USE ED_BATH
  USE ED_BATH_FUNCTIONS
#ifdef _MPI
  USE MPI
  USE SF_MPI
#endif

  implicit none


  integer                               :: Ldelta
  complex(8),dimension(:,:),allocatable :: Gdelta
  complex(8),dimension(:,:),allocatable :: Fdelta
  real(8),dimension(:),allocatable      :: Xdelta,Wdelta
  integer                               :: totNorb,totNspin,totNso
  integer,dimension(:),allocatable      :: getIorb,getJorb,getIspin,getJspin
  integer                               :: Orb_indx,Spin_indx,Spin_mask
  type(effective_bath)                  :: chi2_bath
  integer                               :: cg_iter_count=0
  logical                               :: para_

  integer                               :: MPI_RANK=0
  integer                               :: MPI_SIZE=1
  logical                               :: MPI_MASTER=.true.
  integer                               :: MPI_IERR

  !This contains the number of the lambda expansion
  !for each replica of the impurity
  integer                              :: Nlambdas
  !
  !This is a dummy object which is used here to point
  !to the replica bath lambdas, i.e. the coefficients
  !of the bath item-th Hamiltonian expansion 
  type nsymm_vector
     real(8),dimension(:),allocatable   :: element          
  end type nsymm_vector


END MODULE ED_FIT_COMMON
