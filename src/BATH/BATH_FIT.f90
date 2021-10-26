MODULE ED_BATH_FIT
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
  USE ED_FIT_COMMON
  USE ED_FIT_NORMAL
  USE ED_FIT_HYBRID
  USE ED_FIT_REPLICA
#ifdef _MPI
  USE MPI
  USE SF_MPI
#endif


  implicit none
  private

  interface ed_chi2_fitgf
     module procedure chi2_fitgf_generic_normal
     module procedure chi2_fitgf_sites_normal
  end interface ed_chi2_fitgf

  public :: ed_chi2_fitgf


contains


  !+----------------------------------------------------------------------+
  !PURPOSE  : Chi^2 fit of the G0/Delta 
  ! - CHI2_FITGF_GENERIC_NORMAL interface for the normal case 
  !+----------------------------------------------------------------------+
  subroutine chi2_fitgf_generic_normal(fg,bath,ispin,iorb)
    complex(8),dimension(:,:,:,:,:) :: fg ![Nspin][Nspin][Norb][Norb][Niw] 
    real(8),dimension(:)            :: bath
    integer,optional                :: ispin,iorb
    integer                         :: ispin_
    !
#ifdef _MPI
    if(check_MPI())then
       MPI_MASTER=get_Master_MPI(MPI_COMM_WORLD)
    endif
#endif
    !
    ispin_=1;if(present(ispin))ispin_=ispin
    !
    call assert_shape(fg,[Nspin,Nspin,Norb,Norb,size(fg,5)],"chi2_fitgf_generic_normal","fg")
    select case(cg_method)
    case default
       stop "ED Error: cg_method > 1"
    case (0)
       if(ed_verbose>2)write(LOGfile,"(A,I1,A,A)")"master: Chi^2 fit with CG-nr and CG-weight: ",cg_weight," on: ",cg_scheme
    case (1)
       if(ed_verbose>2)write(LOGfile,"(A,I1,A,A)")"master: Chi^2 fit with CG-minimize and CG-weight: ",cg_weight," on: ",cg_scheme
    end select
    !
    if(MPI_MASTER)then
       select case(bath_type)
       case default
          if(present(iorb))then
             call chi2_fitgf_normal_normal(fg(ispin_,ispin_,:,:,:),bath,ispin_,iorb)
          else
             call chi2_fitgf_normal_normal(fg(ispin_,ispin_,:,:,:),bath,ispin_)
          endif
          !
       case ("hybrid")
          call chi2_fitgf_hybrid_normal(fg(ispin_,ispin_,:,:,:),bath,ispin_)
          !
       case ("replica")
          call chi2_fitgf_replica(fg,bath)
       end select
    endif
    !
#ifdef _MPI
    if(check_MPI())then
       call Bcast_MPI(MPI_COMM_WORLD,bath)
    endif
#endif
    if(.not.MPI_MASTER)write(LOGfile,"(A)")"Bath received from master node"
    !
    !set trim_state_list to true after the first fit has been done: this 
    !marks the ends of the cycle of the 1st DMFT loop.
    trim_state_list=.true.
  end subroutine chi2_fitgf_generic_normal



  !+----------------------------------------------------------------------!
  ! PURPOSE: given a number of independent baths, evaluate N independent
  ! Delta/G0 functions and fit them to update the effective baths for ED.
  !+----------------------------------------------------------------------!
  subroutine chi2_fitgf_sites_normal(fg,bath,Hloc,ispin)
    real(8),intent(inout)    :: bath(:,:)
    complex(8),intent(inout) :: fg(size(bath,1),Nspin,Nspin,Norb,Norb,Lmats)
    real(8)                  :: Hloc(size(bath,1),Nspin,Nspin,Norb,Norb)
    integer,optional         :: ispin
    !MPI auxiliary vars
    real(8)                  :: bath_tmp(size(bath,1),size(bath,2))
    integer                  :: ilat,i,iorb,ispin_
    integer                  :: Nsites
    logical                  :: check_dim
    character(len=5)         :: tmp_suffix
    !
#ifdef _MPI
    if(check_MPI())then
       MPI_RANK = get_Rank_MPI(MPI_COMM_WORLD)
       MPI_SIZE = get_Size_MPI(MPI_COMM_WORLD)
       MPI_MASTER=get_Master_MPI(MPI_COMM_WORLD)
    endif
#endif
    !
    ! Check dimensions !
    Nsites=size(bath,1)
    !
    do ilat=1+MPI_RANK,Nsites,MPI_SIZE
       check_dim = check_bath_dimension(bath(ilat,:))
       if(.not.check_dim) stop "init_lattice_bath: wrong bath size dimension 1 or 2 "
    end do
    !
    bath_tmp=0d0
    do ilat = 1+MPI_RANK,Nsites,MPI_SIZE
       bath_tmp(ilat,:)=bath(ilat,:)
       impHloc = Hloc(ilat,:,:,:,:)
       !
       call ed_set_suffix(ilat)
       !
       if(present(ispin))then
          ispin_=ispin
          if(ispin_>Nspin)stop "ed_fit_bath_sites error: required spin index > Nspin"
          call ed_chi2_fitgf(fg(ilat,:,:,:,:,:),bath_tmp(ilat,:),ispin=ispin_)
       else
          call ed_chi2_fitgf(fg(ilat,:,:,:,:,:),bath_tmp(ilat,:))
       end if
    end do
    !
#ifdef _MPI
    if(check_MPI())then
       bath=0d0
       call MPI_ALLREDUCE(bath_tmp(:,:),bath(:,:),size(bath),&
            MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,MPI_IERR)
    else
       bath = bath_tmp
    endif
#else
    bath = bath_tmp
#endif
    !
    call ed_reset_suffix()
  end subroutine chi2_fitgf_sites_normal










end MODULE ED_BATH_FIT
