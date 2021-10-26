module ED_MAIN
  USE SF_IOTOOLS, only: str,reg
  USE SF_TIMER,only: start_timer,stop_timer
  USE ED_INPUT_VARS
  USE ED_VARS_GLOBAL
  USE ED_EIGENSPACE, only: state_list,es_delete_espace
  USE ED_AUX_FUNX
  USE ED_SETUP
  USE ED_BATH
  USE ED_HAMILTONIAN
  USE ED_GREENS_FUNCTIONS
  USE ED_CHI_FUNCTIONS
  USE ED_OBSERVABLES
  USE ED_DIAG

  implicit none
  private
  !
  !>INIT ED SOLVER
  !
  interface ed_init_solver
     module procedure :: ed_init_solver_single
     module procedure :: ed_init_solver_lattice
  end interface ed_init_solver
  public :: ed_init_solver


  !
  !> ED SOLVER
  !
  interface ed_solve
     module procedure :: ed_solve_single
     module procedure :: ed_solve_lattice
  end interface ed_solve
  public :: ed_solve

  !Boolean to select MPI mode for the inequivalent sites solver:
  !T: solve each site SERIALLY, using MPI lanczos, 
  !F: solve each site PARALLEL, using Serial Lanczos   
  logical :: mpi_lanc_=.true.

contains





  !+-----------------------------------------------------------------------------+!
  ! PURPOSE: allocate and initialize one or multiple baths -+!
  !+-----------------------------------------------------------------------------+!
  !                              SINGLE SITE                                      !
  !+-----------------------------------------------------------------------------+!
  subroutine ed_init_solver_single(bath)
    real(8),dimension(:),intent(inout) :: bath
    logical                            :: check 
    logical,save                       :: isetup=.true.
    integer                            :: i
    !
    !SET THE LOCAL MPI COMMUNICATOR
    if(mpi_lanc_)call ed_set_MpiComm()
    !
    write(LOGfile,"(A)")"INIT SOLVER FOR "//trim(ed_file_suffix)
    !
    !Init ED Structure & memory
    if(isetup)call init_ed_structure() 
    !
    check = check_bath_dimension(bath)
    if(.not.check)stop "init_ed_solver_single error: wrong bath dimensions"
    !
    bath = 0d0
    !
    call allocate_dmft_bath(dmft_bath)
    call init_dmft_bath(dmft_bath)
    call get_dmft_bath(dmft_bath,bath)
    !
    if(isetup)then
       call setup_global
    endif
    call deallocate_dmft_bath(dmft_bath)
    isetup=.false.
    !
    call ed_del_MpiComm()
    !
  end subroutine ed_init_solver_single



  !+-----------------------------------------------------------------------------+!
  !                           INEQUVALENT SITES                                   !
  !+-----------------------------------------------------------------------------+!
  subroutine ed_init_solver_lattice(bath)
    real(8),dimension(:,:) :: bath        ![Nlat,Nb]
    integer                :: ilat,Nineq
    logical                :: check_dim
    integer                :: MPI_ERR
    !
    if(allocated(dens_ineq))deallocate(dens_ineq)
    if(allocated(docc_ineq))deallocate(docc_ineq)
    if(allocated(mag_ineq))deallocate(mag_ineq)
    if(allocated(e_ineq))deallocate(e_ineq)
    if(allocated(dd_ineq))deallocate(dd_ineq)
    if(allocated(Smats_ineq))deallocate(Smats_ineq)
    if(allocated(Sreal_ineq))deallocate(Sreal_ineq)
    if(allocated(Gmats_ineq))deallocate(Gmats_ineq)
    if(allocated(Greal_ineq))deallocate(Greal_ineq)
    if(allocated(Dmats_ineq))deallocate(Dmats_ineq)
    if(allocated(Dreal_ineq))deallocate(Dreal_ineq)
    if(allocated(neigen_sector_ineq))deallocate(neigen_sector_ineq)
    if(allocated(neigen_total_ineq))deallocate(neigen_total_ineq)
    !
    Nineq = size(bath,1)
    if(bath_type=='replica' .AND. .not.allocated(Hreplica_lambda_ineq))&
         stop "ERROR ed_init_solver: replica parameters lambda not defined for all sites"
    !
    allocate(dens_ineq(Nineq,Norb))
    allocate(docc_ineq(Nineq,Norb))
    allocate(mag_ineq(Nineq,Norb))
    allocate(e_ineq(Nineq,4))
    allocate(dd_ineq(Nineq,4))
    !
    allocate(Smats_ineq(Nineq,Nspin,Nspin,Norb,Norb,Lmats))
    allocate(Sreal_ineq(Nineq,Nspin,Nspin,Norb,Norb,Lreal))
    !
    allocate(Gmats_ineq(Nineq,Nspin,Nspin,Norb,Norb,Lmats))
    allocate(Greal_ineq(Nineq,Nspin,Nspin,Norb,Norb,Lreal))
    !
    allocate(Dmats_ineq(Nineq,Lmats))
    allocate(Dreal_ineq(Nineq,Lreal))    
    !
    !
    do ilat=1,Nineq
       call ed_set_suffix(ilat)
       if(bath_type=='replica')call Hreplica_site(ilat)
       call ed_init_solver_single(bath(ilat,:))
    enddo
    !
    allocate(neigen_sector_ineq(Nineq,Nsectors))
    allocate(neigen_total_ineq(Nineq))
    do ilat=1,Nineq
       neigen_sector_ineq(ilat,:) = neigen_sector(:)
       neigen_total_ineq(ilat)    = lanc_nstates_total
    end do
    !
    call ed_reset_suffix
  end subroutine ed_init_solver_lattice










  !##################################################################
  !##################################################################
  !##################################################################
  !##################################################################








  !+-----------------------------------------------------------------------------+!
  !PURPOSE: solve the impurity problems for a single or many independent
  ! lattice site using ED. 
  !+-----------------------------------------------------------------------------+!
  !+-----------------------------------------------------------------------------+!
  !                              SINGLE SITE                                      !
  !+-----------------------------------------------------------------------------+!
  subroutine ed_solve_single(bath,Hloc,sflag)
    real(8),dimension(:),intent(in) :: bath
    real(8),optional,intent(in)     :: Hloc(Nspin,Nspin,Norb,Norb)
    logical,optional                :: sflag
    logical                         :: check,iflag
    !
    iflag=.true. ;if(present(sflag))iflag=sflag
    !
    !SET THE LOCAL MPI COMMUNICATOR
    !Only if mpi_lanc_ is TRUE: solve using MPI Lanczos  
    if(mpi_lanc_)call ed_set_MpiComm()
    !
    if(MpiMaster)call save_input_file(str(ed_input_file))
    !
    call set_Himpurity(Hloc)
    !
    check = check_bath_dimension(bath)
    if(.not.check)stop "ED_SOLVE_SINGLE Error: wrong bath dimensions"
    !
    call allocate_dmft_bath(dmft_bath)
    call set_dmft_bath(bath,dmft_bath)
    call write_dmft_bath(dmft_bath,LOGfile)
    if(MpiMaster)call save_dmft_bath(dmft_bath,used=.true.)
    !
    !
    !SOLVE THE QUANTUM IMPURITY PROBLEM:
    call diagonalize_impurity()
    if(iflag)then
       call buildgf_impurity()
       call buildchi_impurity()
    endif
    call observables_impurity()
    call local_energy_impurity()
    !
    call deallocate_dmft_bath(dmft_bath)
    call es_delete_espace(state_list)
    !
    ! !DELETE THE LOCAL MPI COMMUNICATOR:
    call ed_del_MpiComm()
    !
    nullify(spHtimesV_p)
    write(Logfile,"(A)")""
  end subroutine ed_solve_single



  !+-----------------------------------------------------------------------------+!
  !                          INEQUIVALENT SITES                                   !
  !+-----------------------------------------------------------------------------+!
  subroutine ed_solve_lattice(bath,Hloc,mpi_lanc)
    real(8)          :: bath(:,:) ![Nlat][Nb]
    real(8)          :: Hloc(size(bath,1),Nspin,Nspin,Norb,Norb)
    logical,optional :: mpi_lanc
    !MPI  auxiliary vars
    complex(8)       :: Smats_tmp(size(bath,1),Nspin,Nspin,Norb,Norb,Lmats)
    complex(8)       :: Sreal_tmp(size(bath,1),Nspin,Nspin,Norb,Norb,Lreal)
    complex(8)       :: Gmats_tmp(size(bath,1),Nspin,Nspin,Norb,Norb,Lmats)
    complex(8)       :: Greal_tmp(size(bath,1),Nspin,Nspin,Norb,Norb,Lreal)
    complex(8)       :: Dmats_tmp(size(bath,1),Lmats)
    complex(8)       :: Dreal_tmp(size(bath,1),Lreal)
    real(8)          :: dens_tmp(size(bath,1),Norb)
    real(8)          :: docc_tmp(size(bath,1),Norb)
    real(8)          :: mag_tmp(size(bath,1),Norb)
    real(8)          :: e_tmp(size(bath,1),4)
    real(8)          :: dd_tmp(size(bath,1),4)
    integer          :: neigen_sector_tmp(size(bath,1),Nsectors)
    integer          :: neigen_total_tmp(size(bath,1))
    ! 
    integer          :: i,j,ilat,iorb,jorb,ispin,jspin
    integer          :: Nineq
    logical          :: check_dim
    character(len=5) :: tmp_suffix
    !Local MPI 
    integer          :: MPI_ID=0
    integer          :: MPI_SIZE=1
    logical          :: MPI_MASTER=.true.
    !
    integer          :: mpi_err 
    !
#ifdef _MPI
    if(check_MPI())then
       MPI_ID     = get_Rank_MPI(MPI_COMM_WORLD)
       MPI_SIZE   = get_Size_MPI(MPI_COMM_WORLD)
       MPI_MASTER = get_Master_MPI(MPI_COMM_WORLD)
    endif
#endif
    !
    mpi_lanc_=.false.;if(present(mpi_lanc))mpi_lanc_=mpi_lanc
    !
    ! Check dimensions !
    Nineq=size(bath,1)
    !
    if(size(neigen_sector_ineq,1)<Nineq)stop "ed_solve_lattice error: size(neigen_sectorii,1)<Nineq"
    if(size(neigen_total_ineq)<Nineq)stop "ed_solve_lattice error: size(neigen_totalii,1)<Nineq"
    !
    !Check the dimensions of the bath are ok.
    !This can always be done in parallel no issues with mpi_lanc
    do ilat=1+MPI_ID,Nineq,MPI_SIZE
       check_dim = check_bath_dimension(bath(ilat,:))
       if(.not.check_dim) stop "init_lattice_bath: wrong bath size dimension 1 or 2 "
    end do
    !
    Smats_ineq = zero ; Sreal_ineq = zero
    Gmats_ineq = zero ; Greal_ineq = zero
    Dmats_ineq = zero ; Dreal_ineq = zero 
    dens_ineq  = 0d0  ; docc_ineq  = 0d0
    e_ineq     = 0d0  ; dd_ineq    = 0d0 
    mag_ineq   = 0d0
    neigen_sector_ineq=0 ; neigen_total_ineq =0
    !
    Smats_tmp  = zero ; Sreal_tmp  = zero
    Gmats_tmp  = zero ; Greal_tmp  = zero
    Dmats_tmp  = zero ; Dreal_tmp  = zero
    dens_tmp   = 0d0  ; docc_tmp   = 0d0
    mag_tmp    = 0d0  ; 
    e_tmp      = 0d0  ; dd_tmp     = 0d0
    neigen_sector_tmp = 0
    neigen_total_tmp  = 0
    !
    select case(mpi_lanc_)
    case (.false.)             !mpi_lanc=False => solve sites with MPI
       if(MPI_MASTER)call start_timer
       do ilat = 1 + MPI_ID, Nineq, MPI_SIZE
          write(LOGfile,*)str(MPI_ID)//" SOLVES INEQ SITE: "//str(ilat,Npad=site_indx_padding)
          call ed_set_suffix(ilat)
          !
          !Solve the impurity problem for the ilat-th site
          neigen_sector(:)   = neigen_sector_ineq(ilat,:)
          lanc_nstates_total = neigen_total_ineq(ilat)
          !
          !Call ed_solve in SERIAL MODE: This is parallel on the ineq. sites
          call ed_solve_single(bath(ilat,:),Hloc(ilat,:,:,:,:))
          !
          neigen_sector_tmp(ilat,:)  = neigen_sector(:)
          neigen_total_tmp(ilat)     = lanc_nstates_total
          Smats_tmp(ilat,:,:,:,:,:)  = impSmats(:,:,:,:,:)
          Sreal_tmp(ilat,:,:,:,:,:)  = impSreal(:,:,:,:,:)
          Gmats_tmp(ilat,:,:,:,:,:)  = impGmats(:,:,:,:,:)
          Greal_tmp(ilat,:,:,:,:,:)  = impGreal(:,:,:,:,:)
          Dmats_tmp(ilat,:)          = impDmats(:)
          Dreal_tmp(ilat,:)          = impDreal(:)
          dens_tmp(ilat,1:Norb)      = ed_dens(1:Norb)
          docc_tmp(ilat,1:Norb)      = ed_docc(1:Norb)
          mag_tmp(ilat,1:Norb)       = ed_mag(1:Norb)
          e_tmp(ilat,:)              = [ed_Epot,ed_Eint,ed_Ehartree,ed_Eknot]
          dd_tmp(ilat,:)             = [ed_Dust,ed_Dund,ed_Dse,ed_Dph]
       enddo
#ifdef _MPI
       if(check_MPI())then
          call MPI_Barrier(MPI_COMM_WORLD,MPI_ERR)
       endif
#endif
       if(MPI_MASTER)call stop_timer(unit=LOGfile)
       call ed_reset_suffix
       !
#ifdef _MPI
       if(check_MPI())then
          call AllReduce_MPI(MPI_COMM_WORLD,Smats_tmp,Smats_ineq)
          call AllReduce_MPI(MPI_COMM_WORLD,Sreal_tmp,Sreal_ineq)
          call AllReduce_MPI(MPI_COMM_WORLD,Gmats_tmp,Gmats_ineq)
          call AllReduce_MPI(MPI_COMM_WORLD,Greal_tmp,Greal_ineq)
          call AllReduce_MPI(MPI_COMM_WORLD,Dmats_tmp,Dmats_ineq)
          call AllReduce_MPI(MPI_COMM_WORLD,Dreal_tmp,Dreal_ineq)
          call AllReduce_MPI(MPI_COMM_WORLD,dens_tmp,dens_ineq)
          call AllReduce_MPI(MPI_COMM_WORLD,docc_tmp,docc_ineq)
          call AllReduce_MPI(MPI_COMM_WORLD,mag_tmp,mag_ineq)
          call AllReduce_MPI(MPI_COMM_WORLD,e_tmp,e_ineq)
          call AllReduce_MPI(MPI_COMM_WORLD,dd_tmp,dd_ineq)
          neigen_sector_ineq=0
          neigen_total_ineq=0
          call AllReduce_MPI(MPI_COMM_WORLD,neigen_sector_tmp,neigen_sector_ineq)
          call AllReduce_MPI(MPI_COMM_WORLD,neigen_total_tmp,neigen_total_ineq)
       else
          Smats_ineq = Smats_tmp
          Sreal_ineq = Sreal_tmp
          Gmats_ineq = Gmats_tmp
          Greal_ineq = Greal_tmp
          Dmats_ineq = Dmats_tmp
          Dreal_ineq = Dreal_tmp
          dens_ineq  = dens_tmp
          docc_ineq  = docc_tmp
          mag_ineq   = mag_tmp
          e_ineq     = e_tmp
          dd_ineq    = dd_tmp
          neigen_sector_ineq = neigen_sector_tmp
          neigen_total_ineq  = neigen_total_tmp
       endif
#else
       Smats_ineq = Smats_tmp
       Sreal_ineq = Sreal_tmp
       Gmats_ineq = Gmats_tmp
       Greal_ineq = Greal_tmp
       Dmats_ineq = Dmats_tmp
       Dreal_ineq = Dreal_tmp
       dens_ineq  = dens_tmp
       docc_ineq  = docc_tmp
       mag_ineq   = mag_tmp
       e_ineq     = e_tmp
       dd_ineq    = dd_tmp
       neigen_sector_ineq = neigen_sector_tmp
       neigen_total_ineq  = neigen_total_tmp
#endif
       !       
    case(.true.)                !solve sites serial, Lanczos with MPI
       if(MPI_MASTER)call start_timer
       do ilat = 1, Nineq
          write(LOGfile,*)" SOLVES INEQ SITE: "//str(ilat,Npad=site_indx_padding)
          call ed_set_suffix(ilat)
          !
          !Solve the impurity problem for the ilat-th site
          neigen_sector(:)   = neigen_sector_ineq(ilat,:)
          lanc_nstates_total = neigen_total_ineq(ilat)
          !
          call ed_solve_single(bath(ilat,:),Hloc(ilat,:,:,:,:))
          !
          neigen_sector_ineq(ilat,:)  = neigen_sector(:)
          neigen_total_ineq(ilat)     = lanc_nstates_total
          Smats_ineq(ilat,:,:,:,:,:)  = impSmats(:,:,:,:,:)
          Sreal_ineq(ilat,:,:,:,:,:)  = impSreal(:,:,:,:,:)
          Gmats_ineq(ilat,:,:,:,:,:)  = impGmats(:,:,:,:,:)
          Greal_ineq(ilat,:,:,:,:,:)  = impGreal(:,:,:,:,:)
          Dmats_ineq(ilat,:)          = impDmats(:)
          Dreal_ineq(ilat,:)          = impDreal(:)
          dens_ineq(ilat,1:Norb)      = ed_dens(1:Norb)
          docc_ineq(ilat,1:Norb)      = ed_docc(1:Norb)
          mag_ineq(ilat,1:Norb)       = ed_mag(1:Norb)
          e_ineq(ilat,:)              = [ed_Epot,ed_Eint,ed_Ehartree,ed_Eknot]
          dd_ineq(ilat,:)             = [ed_Dust,ed_Dund,ed_Dse,ed_Dph]
       enddo
       if(MPI_MASTER)call stop_timer(unit=LOGfile)
       call ed_reset_suffix
    end select
    !
  end subroutine ed_solve_lattice







end module ED_MAIN







