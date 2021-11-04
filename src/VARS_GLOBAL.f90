MODULE ED_VARS_GLOBAL
  USE SF_CONSTANTS
  USE ED_SPARSE_MATRIX
#ifdef _MPI
  USE MPI
  USE SF_MPI
#endif
  implicit none






  !-------------------- EFFECTIVE BATH STRUCTURE ----------------------!
  type H_operator
     real(8),dimension(:,:,:,:),allocatable                  :: O          !Replica hamilt
  end type H_operator

  type effective_bath_component
     real(8)                                                 :: v
     real(8),dimension(:),allocatable                        :: lambda![Nsym]
  end type effective_bath_component

  type effective_bath
     real(8),dimension(:,:,:),allocatable                    :: e     !local energies [Nspin][Norb][Nbath]/[Nspin][1][Nbath]_hybrid
     real(8),dimension(:,:,:),allocatable                    :: v     !spin-keep hyb. [Nspin][Norb][Nbath]
     integer                                                 :: Nbasis  !H Basis dimension  
     type(effective_bath_component),dimension(:),allocatable :: item  ![Nbath] Replica bath components, V included
     logical                                                 :: status=.false.
  end type effective_bath






  !---------------- SECTOR-TO-FOCK SPACE STRUCTURE -------------------!

  type sector_map
     integer,dimension(:),allocatable :: map
     logical                          :: status=.false.
  end type sector_map

  type sector
     integer                                     :: index       !
     type(sector_map),dimension(:),allocatable   :: H
     integer,dimension(:),allocatable            :: DimUps
     integer,dimension(:),allocatable            :: DimDws
     integer                                     :: DimUp
     integer                                     :: DimDw
     integer                                     :: DimEl
     integer                                     :: DimPh
     integer                                     :: Dim
     integer,dimension(:),allocatable            :: Nups
     integer,dimension(:),allocatable            :: Ndws
     integer                                     :: Nup
     integer                                     :: Ndw
     integer                                     :: Nlanc
     logical                                     :: status=.false.
  end type sector



  !------------------ ABTRACT INTERFACES PROCEDURES ------------------!
  !SPARSE MATRIX-VECTOR PRODUCTS USED IN ED_MATVEC
  !dbleMat*dbleVec
  abstract interface
     subroutine dd_sparse_HxV(Nloc,v,Hv)
       integer                 :: Nloc
       real(8),dimension(Nloc) :: v
       real(8),dimension(Nloc) :: Hv
     end subroutine dd_sparse_HxV
  end interface




  !-------------------------- ED  VARIABLES --------------------------!

  !SIZE OF THE PROBLEM
  !=========================================================
  integer,save                                       :: Ns       !Number of levels per spin
  integer,save                                       :: Nsectors !Number of sectors
  integer,save                                       :: Ns_orb
  integer,save                                       :: Ns_ud
  integer,save                                       :: DimPh    !Number of phonon states


  !Some maps between sectors and full Hilbert space (pointers)
  !PRIVATE:
  !=========================================================
  integer,allocatable,dimension(:)                   :: getDim             ! [Nsectors]
  integer,allocatable,dimension(:,:,:)               :: getCsector         ! [1/Norb,2,NSectors]
  integer,allocatable,dimension(:,:,:)               :: getCDGsector       ! [1/Norb,2,NSectors]
  integer,allocatable,dimension(:,:)                 :: getBathStride
  integer,allocatable,dimension(:,:)                 :: impIndex
  logical,allocatable,dimension(:)                   :: twin_mask
  logical,allocatable,dimension(:)                   :: sectors_mask

  !Effective Bath used in the ED code (this is opaque to user)
  !=========================================================
  type(effective_bath)                               :: dmft_bath

  !Replica bath basis set
  !=========================================================
  type(H_operator),dimension(:),allocatable          :: Hreplica_basis
  real(8),dimension(:),allocatable                   :: Hreplica_lambda
  integer                                            :: Hreplica_Nsym=0
  logical                                            :: Hreplica_status=.false.

  !local part of the Hamiltonian
  !=========================================================
  real(8),dimension(:,:,:,:),allocatable             :: impHloc           !local hamiltonian

  !Variables for DIAGONALIZATION
  !PRIVATE
  !=========================================================  
  type(sparse_matrix_csr)                            :: spH0d !diagonal part
  type(sparse_matrix_csr)                            :: spH0nd !non-diagonal part
  type(sparse_matrix_csr),dimension(:),allocatable   :: spH0ups,spH0dws !reduced UP and DW parts
  type(sparse_matrix_csr)                            :: spH0_ph !Hamiltonian for phonons
  type(sparse_matrix_csr)                            :: spH0e_eph, spH0edw_eph, spH0ph_eph !electron-phonon interaction
  !
  procedure(dd_sparse_HxV),pointer                   :: spHtimesV_p=>null()


  !Variables for DIAGONALIZATION
  !PRIVATE
  !=========================================================  
  integer,allocatable,dimension(:)                   :: neigen_sector
  logical                                            :: trim_state_list=.false.

  !Partition function
  !PRIVATE
  !=========================================================
  real(8)                                            :: zeta_function
  real(8)                                            :: gs_energy


  logical                                            :: finiteT             !flag for finite temperature calculation
  

  !Impurity Green's function and Self-Energies: (Nspin,Nspin,Norb,Norb,:)
  !PRIVATE (now public but accessible thru routine)
  !=========================================================
  complex(8),allocatable,dimension(:,:,:,:,:)        :: impSmats
  complex(8),allocatable,dimension(:,:,:,:,:)        :: impSreal
  complex(8),allocatable,dimension(:,:,:,:,:)        :: impGmats
  complex(8),allocatable,dimension(:,:,:,:,:)        :: impGreal
  complex(8),allocatable,dimension(:,:,:,:,:)        :: impG0mats
  complex(8),allocatable,dimension(:,:,:,:,:)        :: impG0real
  complex(8),allocatable,dimension(:)                :: impDmats
  complex(8),allocatable,dimension(:)                :: impDreal

  !Spin Susceptibilities
  !=========================================================
  real(8),allocatable,dimension(:,:,:)               :: spinChi_tau
  complex(8),allocatable,dimension(:,:,:)            :: spinChi_w
  complex(8),allocatable,dimension(:,:,:)            :: spinChi_iv


  !Diagonal/Off-diagonal charge-charge Susceptibilities
  !=========================================================  
  real(8),allocatable,dimension(:,:,:)               :: densChi_tau
  complex(8),allocatable,dimension(:,:,:)            :: densChi_w
  complex(8),allocatable,dimension(:,:,:)            :: densChi_iv

  !Pair-Pair Susceptibilities
  !=========================================================
  real(8),allocatable,dimension(:,:,:)               :: pairChi_tau
  complex(8),allocatable,dimension(:,:,:)            :: pairChi_w
  complex(8),allocatable,dimension(:,:,:)            :: pairChi_iv

  !Exciton Susceptibilities
  !=========================================================
  real(8),allocatable,dimension(:,:,:,:)             :: exctChi_tau ![0:4,:]
  complex(8),allocatable,dimension(:,:,:,:)          :: exctChi_w
  complex(8),allocatable,dimension(:,:,:,:)          :: exctChi_iv


  !Density and double occupancy
  !PRIVATE (now public but accessible thru routines)
  !=========================================================
  real(8),dimension(:),allocatable                   :: ed_dens
  real(8),dimension(:),allocatable                   :: ed_dens_up,ed_dens_dw
  real(8),dimension(:),allocatable                   :: ed_docc
  real(8),dimension(:),allocatable                   :: ed_mag
  real(8)                                            :: ed_Ekin
  real(8)                                            :: ed_Epot
  real(8)                                            :: ed_Eint
  real(8)                                            :: ed_Ehartree
  real(8)                                            :: ed_Eknot
  real(8)                                            :: ed_Dust,ed_Dund,ed_Dse,ed_Dph



  !Frequency and time arrays:
  !=========================================================
  real(8),dimension(:),allocatable                   :: wm,tau,wr,vm,vr




  !--------------- LATTICE WRAP VARIABLES -----------------!
  complex(8),dimension(:,:,:,:,:,:),allocatable,save :: Smats_ineq,Sreal_ineq  ![Nlat,Nspin,Nspin,Norb,Norb,L]
  complex(8),dimension(:,:,:,:,:,:),allocatable,save :: Gmats_ineq,Greal_ineq
  complex(8),dimension(:,:,:,:,:,:),allocatable,save :: G0mats_ineq,G0real_ineq
  complex(8),dimension(:,:),allocatable,save         :: Dmats_ineq,Dreal_ineq
  real(8),dimension(:,:),allocatable,save            :: dens_ineq 
  real(8),dimension(:,:),allocatable,save            :: docc_ineq
  real(8),dimension(:,:),allocatable,save            :: mag_ineq
  real(8),dimension(:,:),allocatable,save            :: dd_ineq,e_ineq
  integer,allocatable,dimension(:,:)                 :: neigen_sector_ineq
  integer,allocatable,dimension(:)                   :: neigen_total_ineq
  real(8),dimension(:,:),allocatable                 :: Hreplica_lambda_ineq

  !File suffixes for printing fine tuning.
  !=========================================================
  character(len=32)                                  :: ed_file_suffix=""       !suffix string attached to the output files.
  character(len=10)                                  :: ineq_site_suffix="_ineq"
  integer                                            :: site_indx_padding=4
  logical                                            :: Jhflag              !spin-exchange and pair-hopping flag.
  logical                                            :: offdiag_gf_flag=.false.
  ! character(len=200)                               :: ed_input_file=""



  logical              :: cg_minimize_ver=.false.
  real(8)              :: cg_minimize_hh=1d-4

  !This is the internal Mpi Communicator and variables.
  !=========================================================
#ifdef _MPI
  integer                                            :: MpiComm_Global=MPI_COMM_NULL
  integer                                            :: MpiComm=MPI_COMM_NULL
#endif
  integer                                            :: MpiGroup_Global=MPI_GROUP_NULL
  integer                                            :: MpiGroup=MPI_GROUP_NULL
  logical                                            :: MpiStatus=.false.
  logical                                            :: MpiMaster=.true.
  integer                                            :: MpiRank=0
  integer                                            :: MpiSize=1
  integer,allocatable,dimension(:)                   :: MpiMembers
  integer                                            :: mpiQup=0
  integer                                            :: mpiRup=0
  integer                                            :: mpiQdw=0
  integer                                            :: mpiRdw=0
  integer                                            :: mpiQ=0
  integer                                            :: mpiR=0
  integer                                            :: mpiIstart
  integer                                            :: mpiIend
  integer                                            :: mpiIshift
  logical                                            :: mpiAllThreads=.true.
  !

contains






  !IF code is compiled with MPI support
  !+  MPI is initialized:
  !THEN this routine setup the internal communicator
  !(inherited from MPI_COMM_WORLD) plus global variables
  !ELSE it does nothing
  !
  !
  subroutine ed_set_MpiComm()
#ifdef _MPI
    integer :: ierr
    if(check_MPI())then
       MpiComm_Global = MPI_COMM_WORLD
       MpiComm        = MPI_COMM_WORLD
       call Mpi_Comm_group(MpiComm_Global,MpiGroup_Global,ierr)
       MpiStatus      = .true.
       MpiSize        = get_Size_MPI(MpiComm_Global)
       MpiRank        = get_Rank_MPI(MpiComm_Global)
       MpiMaster      = get_Master_MPI(MpiComm_Global)
    endif
#endif
  end subroutine ed_set_MpiComm


  !IF code is compiled with MPI support
  !THEN this routine reset global variables to default values (SERIAL)
  subroutine ed_del_MpiComm()
#ifdef _MPI    
    MpiComm_Global = MPI_COMM_NULL
    MpiComm        = MPI_COMM_NULL
    MpiGroup_Global= MPI_GROUP_NULL
    MpiStatus      = .false.
    MpiSize        = 1
    MpiRank        = 0
    MpiMaster      = .true.
#endif
  end subroutine ed_del_MpiComm




END MODULE ED_VARS_GLOBAL
