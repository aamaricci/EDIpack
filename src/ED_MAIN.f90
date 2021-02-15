module ED_MAIN
  USE SF_IOTOOLS, only: str,reg
  USE SF_TIMER,only: start_timer,stop_timer
  USE ED_INPUT_VARS
  USE ED_VARS_GLOBAL
  USE ED_EIGENSPACE, only: state_list,es_delete_espace
  USE ED_AUX_FUNX
  USE ED_HLOC_DECOMPOSITION
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
#ifdef _MPI
     module procedure :: ed_init_solver_single_mpi
     module procedure :: ed_init_solver_lattice_mpi
#endif
  end interface ed_init_solver
  !>
  public :: ed_init_solver


  !
  !> ED SOLVER
  !
  interface ed_solve
     module procedure :: ed_solve_single
     module procedure :: ed_solve_lattice
#ifdef _MPI
     module procedure :: ed_solve_single_mpi
     module procedure :: ed_solve_lattice_mpi
#endif
  end interface ed_solve
  !>
  public :: ed_solve

  character(len=64)                                  :: suffix



contains





  !+-----------------------------------------------------------------------------+!
  ! PURPOSE: allocate and initialize one or multiple baths -+!
  !+-----------------------------------------------------------------------------+!
  !                              SINGLE SITE                                      !
  !+-----------------------------------------------------------------------------+!
  subroutine ed_init_solver_single(bath,Hloc)
    real(8),dimension(:),intent(inout) :: bath
    real(8),intent(in),optional        :: Hloc(Nspin,Nspin,Norb,Norb)
    logical                            :: check 
    logical,save                       :: isetup=.true.
    integer                            :: i
    logical                            :: MPI_MASTER=.true.
    integer                            :: MPI_ERR
    !
    write(LOGfile,"(A)")"INIT SOLVER FOR "//trim(ed_file_suffix)
    !
    !Init ED Structure & memory
    if(isetup)call init_ed_structure()
    !
    !Init bath:
    if(present(Hloc))then
       if(bath_type/="replica")call set_Hloc(Hloc)
    else
       if(.not.allocated(impHloc))then
          print*,"ed_init ERROR: impHloc not allocated. requires calling set_Hloc befor ed_init"
          stop
       endif
    endif
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
  end subroutine ed_init_solver_single

#ifdef _MPI
  subroutine ed_init_solver_single_mpi(MpiComm,bath,Hloc)
    integer                            :: MpiComm
    real(8),dimension(:),intent(inout) :: bath
    real(8),intent(in),optional        :: Hloc(Nspin,Nspin,Norb,Norb)
    logical                            :: check 
    logical,save                       :: isetup=.true.
    integer                            :: i
    !
    !
    !SET THE LOCAL MPI COMMUNICATOR :
    call ed_set_MpiComm(MpiComm)
    !
    write(LOGfile,"(A)")"INIT SOLVER FOR "//trim(ed_file_suffix)
    !
    !Init ED Structure & memory
    if(isetup)call init_ed_structure()
    !
    !Init bath:
    if(present(Hloc))then
       if(bath_type/="replica")call set_Hloc(Hloc)
    else
       if(.not.allocated(impHloc))then
          print*,"ed_init ERROR: impHloc not allocated. requires calling set_Hloc befor ed_init"
          stop
       endif
    endif
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
  end subroutine ed_init_solver_single_mpi
#endif








  !+-----------------------------------------------------------------------------+!
  !                           INEQUVALENT SITES                                   !
  !+-----------------------------------------------------------------------------+!
  subroutine ed_init_solver_lattice(bath,Hloc)
    real(8),dimension(:,:)         :: bath ![Nlat][:]
    real(8),intent(in)             :: Hloc(size(bath,1),Nspin,Nspin,Norb,Norb)
    integer                        :: ilat,Nineq
    logical                        :: check_dim
    character(len=5)               :: tmp_suffix
    integer                        :: MPI_ERR
    !
    !
    Nineq = size(bath,1)
    !
    do ilat=1,Nineq             !all nodes check the bath, u never know...
       !
       ed_file_suffix=reg(ineq_site_suffix)//str(ilat,site_indx_padding)
       !
       call ed_init_solver_single(bath(ilat,:),Hloc(ilat,:,:,:,:))
       !
    end do
    !
    call ineq_memory_allocation(Nineq)
    !
    do ilat=1,Nineq             !all nodes check the bath, u never know...
       neigen_sector_ineq(ilat,:) = neigen_sector(:)
       neigen_total_ineq(ilat)    = lanc_nstates_total
    end do
    !
    ed_file_suffix=""
    !
  end subroutine ed_init_solver_lattice

#ifdef _MPI
  subroutine ed_init_solver_lattice_mpi(MpiComm,bath,Hloc)
    integer                        :: MpiComm
    real(8),dimension(:,:)         :: bath ![Nlat][:]
    real(8),intent(in)             :: Hloc(size(bath,1),Nspin,Nspin,Norb,Norb)
    integer                        :: ilat,Nineq,Nsect
    logical                        :: check_dim
    character(len=5)               :: tmp_suffix
    integer                        :: MPI_ERR
    !
    !
    Nineq = size(bath,1)
    !
    do ilat=1,Nineq             !all nodes check the bath, u never know...
       !
       ed_file_suffix=reg(ineq_site_suffix)//str(ilat,site_indx_padding)
       !
       call ed_init_solver_single(bath(ilat,:),Hloc(ilat,:,:,:,:))
       !
    end do
    !
    call ineq_memory_allocation(Nineq)
    !
    do ilat=1,Nineq             !all nodes check the bath, u never know...
       neigen_sector_ineq(ilat,:) = neigen_sector(:)
       neigen_total_ineq(ilat)    = lanc_nstates_total
    end do
    !
    call MPI_Barrier(MpiComm,MPI_ERR)
    !
    ed_file_suffix=""
    !
  end subroutine ed_init_solver_lattice_mpi
#endif
















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
    if(MpiMaster)call save_input_file(str(ed_input_file))
    !
    if(present(Hloc).AND.(bath_type/="replica"))call set_Hloc(Hloc)
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
    call diagonalize_impurity()         !find target states by digonalization of Hamiltonian
    if(iflag)then
       call buildgf_impurity()             !build the one-particle impurity Green's functions  & Self-energy
       call buildchi_impurity()            !build the local susceptibilities (spin [todo charge])
    endif
    call observables_impurity()         !obtain impurity observables as thermal averages.          
    call local_energy_impurity()        !obtain the local energy of the effective impurity problem
    !
    call deallocate_dmft_bath(dmft_bath)
    call es_delete_espace(state_list)
    !
    nullify(spHtimesV_p)
  end subroutine ed_solve_single



#ifdef _MPI
  !+-----------------------------------------------------------------------------+!
  !                              SINGLE SITE                                      !
  !+-----------------------------------------------------------------------------+!
  subroutine ed_solve_single_mpi(MpiComm,bath,Hloc,sflag)
    integer                         :: MpiComm
    real(8),dimension(:),intent(in) :: bath
    real(8),optional,intent(in)     :: Hloc(Nspin,Nspin,Norb,Norb)
    logical,optional                :: sflag
    logical                         :: check,iflag
    !
    iflag=.true. ;if(present(sflag))iflag=sflag
    !
    !SET THE LOCAL MPI COMMUNICATOR :
    call ed_set_MpiComm(MpiComm)
    !
    if(MpiMaster)call save_input_file(str(ed_input_file))
    !
    if(present(Hloc).AND.(bath_type/="replica"))call set_Hloc(Hloc)
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
    call diagonalize_impurity()         !find target states by digonalization of Hamiltonian
    if(iflag)then
       call buildgf_impurity()             !build the one-particle impurity Green's functions  & Self-energy
       call buildchi_impurity() !build the local susceptibilities (spin [todo charge])
    endif
    call observables_impurity()         !obtain impurity observables as thermal averages.
    call local_energy_impurity()        !obtain the local energy of the effective impurity problem
    !
    call deallocate_dmft_bath(dmft_bath)
    call es_delete_espace(state_list)
    !
    !DELETE THE LOCAL MPI COMMUNICATOR:
    call ed_del_MpiComm()
    !
    nullify(spHtimesV_p)
  end subroutine ed_solve_single_mpi
#endif













  !+-----------------------------------------------------------------------------+!
  !                          INEQUIVALENT SITES                                   !
  !+-----------------------------------------------------------------------------+!
  subroutine ed_solve_lattice(bath,Hloc,Uloc_ii,Ust_ii,Jh_ii,Jp_ii,Jx_ii)
    !inputs
    real(8)          :: bath(:,:) ![Nlat][Nb]
    real(8)          :: Hloc(size(bath,1),Nspin,Nspin,Norb,Norb)
    real(8),optional :: Uloc_ii(size(bath,1),Norb)
    real(8),optional :: Ust_ii(size(bath,1))
    real(8),optional :: Jh_ii(size(bath,1))
    real(8),optional :: Jp_ii(size(bath,1))
    real(8),optional :: Jx_ii(size(bath,1))
    ! 
    integer          :: i,j,ilat,iorb,jorb,ispin,jspin
    integer          :: Nineq
    logical          :: check_dim
    character(len=5) :: tmp_suffix
    !
    logical          :: MPI_MASTER=.true.
    !
    ! Check dimensions !
    Nineq=size(bath,1)
    !
    if(size(neigen_sector_ineq,1)<Nineq)stop "ed_solve_lattice error: size(neigen_sector_ineq,1)<Nineq"
    if(size(neigen_total_ineq)<Nineq)stop "ed_solve_lattice error: size(neigen_total_ineq,1)<Nineq"
    !
    !Check the dimensions of the bath are ok:
    do ilat=1,Nineq
       check_dim = check_bath_dimension(bath(ilat,:))
       if(.not.check_dim) stop "init_lattice_bath: wrong bath size dimension 1 or 2 "
    end do
    Smats_ineq    = zero 
    Sreal_ineq    = zero 
    Gmats_ineq    = zero 
    Greal_ineq    = zero 
    Dmats_ph_ineq = zero 
    Dreal_ph_ineq = zero 
    dens_ineq     = 0d0  
    docc_ineq     = 0d0  
    mag_ineq      = 0d0  
    e_ineq        = 0d0  
    dd_ineq       = 0d0 
    imp_density_matrix_ineq = zero
    !
    call start_timer
    !
    do ilat = 1, Nineq
       write(LOGfile,*)" SOLVING INEQ SITE: "//str(ilat,Npad=4)
       !
       ed_file_suffix=reg(ineq_site_suffix)//str(ilat,site_indx_padding)
       !
       !If required set the local value of U per each site
       if(present(Uloc_ii))Uloc(1:Norb) = Uloc_ii(ilat,1:Norb)
       if(present(Ust_ii)) Ust = Ust_ii(ilat)
       if(present(Jh_ii))  Jh  = Jh_ii(ilat)
       if(present(Jp_ii))  Jp  = Jp_ii(ilat)
       if(present(Jx_ii))  Jx  = Jx_ii(ilat)
       !
       !Solve the impurity problem for the ilat-th site
       neigen_sector(:)   = neigen_sector_ineq(ilat,:)
       lanc_nstates_total = neigen_total_ineq(ilat)
       !
       !Call ed_solve in SERIAL MODE!! This is parallel on the ineq. sites
       call ed_solve_single(bath(ilat,:),Hloc(ilat,:,:,:,:))
       !
       neigen_sector_ineq(ilat,:)  = neigen_sector(:)
       neigen_total_ineq(ilat)     = lanc_nstates_total
       Smats_ineq(ilat,:,:,:,:,:)  = impSmats(:,:,:,:,:)
       Sreal_ineq(ilat,:,:,:,:,:)  = impSreal(:,:,:,:,:)
       Gmats_ineq(ilat,:,:,:,:,:)  = impGmats(:,:,:,:,:)
       Greal_ineq(ilat,:,:,:,:,:)  = impGreal(:,:,:,:,:)
       Dmats_ph_ineq(ilat,:)       = impDmats_ph(:)
       Dreal_ph_ineq(ilat,:)       = impDreal_ph(:)
       dens_ineq(ilat,1:Norb)      = ed_dens(1:Norb)
       docc_ineq(ilat,1:Norb)      = ed_docc(1:Norb)
       mag_ineq(ilat,1:Norb)       = ed_mag(1:Norb)
       e_ineq(ilat,:)              = [ed_Epot,ed_Eint,ed_Ehartree,ed_Eknot]
       dd_ineq(ilat,:)             = [ed_Dust,ed_Dund,ed_Dse,ed_Dph]
       imp_density_matrix_ineq(ilat,:,:,:,:) = imp_density_matrix(:,:,:,:)
       !
    enddo
    !
    call stop_timer(unit=LOGfile)
    !
    ed_file_suffix=""
    !
  end subroutine ed_solve_lattice


  !FALL BACK: DO A VERSION THAT DOES THE SITES IN PARALLEL USING SERIAL ED CODE
#ifdef _MPI
  subroutine ed_solve_lattice_mpi(MpiComm,bath,Hloc,mpi_lanc,Uloc_ii,Ust_ii,Jh_ii,Jp_ii,Jx_ii)
    integer          :: MpiComm
    !inputs
    real(8)          :: bath(:,:) ![Nlat][Nb]
    real(8)          :: Hloc(size(bath,1),Nspin,Nspin,Norb,Norb)
    logical,optional :: mpi_lanc
    real(8),optional :: Uloc_ii(size(bath,1),Norb)
    real(8),optional :: Ust_ii(size(bath,1))
    real(8),optional :: Jh_ii(size(bath,1))
    real(8),optional :: Jp_ii(size(bath,1))
    real(8),optional :: Jx_ii(size(bath,1))
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
    !    
    complex(8)       :: imp_density_matrix_tmp(size(bath,1),Nspin,Nspin,Norb,Norb)
    !
    integer          :: neigen_sectortmp(size(bath,1),Nsectors)
    integer          :: neigen_totaltmp(size(bath,1))
    ! 
    integer          :: i,j,ilat,iorb,jorb,ispin,jspin
    integer          :: Nineq
    logical          :: check_dim,mpi_lanc_
    character(len=5) :: tmp_suffix
    !
    integer          :: MPI_ID=0
    integer          :: MPI_SIZE=1
    logical          :: MPI_MASTER=.true.
    !
    integer          :: mpi_err 
    !
    MPI_ID     = get_Rank_MPI(MpiComm)
    MPI_SIZE   = get_Size_MPI(MpiComm)
    MPI_MASTER = get_Master_MPI(MpiComm)
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
    Smats_ineq    = zero ; Sreal_ineq    = zero
    Gmats_ineq    = zero ; Greal_ineq    = zero
    Dmats_ph_ineq = zero ; Dreal_ph_ineq = zero 
    dens_ineq     = 0d0  ; docc_ineq     = 0d0
    mag_ineq      = 0d0
    e_ineq        = 0d0  ; dd_ineq       = 0d0 
    neigen_sector_ineq=0 ; neigen_total_ineq =0
    imp_density_matrix_ineq = zero
    !
    Smats_tmp  = zero ; Sreal_tmp  = zero
    Gmats_tmp  = zero ; Greal_tmp  = zero
    Dmats_tmp  = zero ; Dreal_tmp  = zero
    dens_tmp   = 0d0  ; docc_tmp   = 0d0
    mag_tmp    = 0d0  ; 
    e_tmp      = 0d0  ; dd_tmp     = 0d0
    neigen_sectortmp = 0
    neigen_totaltmp  = 0
    imp_density_matrix_tmp = zero
    !
    select case(mpi_lanc_)
    case default              !mpi_lanc=False => solve sites with MPI
       if(MPI_MASTER)call start_timer
       do ilat = 1 + MPI_ID, Nineq, MPI_SIZE
          write(LOGfile,*)str(MPI_ID)//" SOLVES INEQ SITE: "//str(ilat,Npad=4)
          !
          ed_file_suffix=reg(ineq_site_suffix)//str(ilat,site_indx_padding)
          !
          !If required set the local value of U per each site
          if(present(Uloc_ii))Uloc(1:Norb) = Uloc_ii(ilat,1:Norb)
          if(present(Ust_ii)) Ust = Ust_ii(ilat)
          if(present(Jh_ii))  Jh  = Jh_ii(ilat)
          if(present(Jp_ii))  Jp  = Jp_ii(ilat)
          if(present(Jx_ii))  Jx  = Jx_ii(ilat)
          !
          !Solve the impurity problem for the ilat-th site
          neigen_sector(:)   = neigen_sector_ineq(ilat,:)
          lanc_nstates_total = neigen_total_ineq(ilat)
          !
          !Call ed_solve in SERIAL MODE!! This is parallel on the ineq. sites
          call ed_solve_single(bath(ilat,:),Hloc(ilat,:,:,:,:))
          !
          neigen_sectortmp(ilat,:)   = neigen_sector(:)
          neigen_totaltmp(ilat)      = lanc_nstates_total
          Smats_tmp(ilat,:,:,:,:,:)  = impSmats(:,:,:,:,:)
          Sreal_tmp(ilat,:,:,:,:,:)  = impSreal(:,:,:,:,:)
          Gmats_tmp(ilat,:,:,:,:,:)  = impGmats(:,:,:,:,:)
          Greal_tmp(ilat,:,:,:,:,:)  = impGreal(:,:,:,:,:)
          Dmats_tmp(ilat,:)          = impDmats_ph(:)
          Dreal_tmp(ilat,:)          = impDreal_ph(:)
          dens_tmp(ilat,1:Norb)      = ed_dens(1:Norb)
          docc_tmp(ilat,1:Norb)      = ed_docc(1:Norb)
          mag_tmp(ilat,1:Norb)       = ed_mag(1:Norb)
          e_tmp(ilat,:)              = [ed_Epot,ed_Eint,ed_Ehartree,ed_Eknot]
          dd_tmp(ilat,:)             = [ed_Dust,ed_Dund,ed_Dse,ed_Dph]
          imp_density_matrix_tmp(ilat,:,:,:,:) = imp_density_matrix(:,:,:,:)
       enddo
       call MPI_Barrier(MpiComm,MPI_ERR)
       if(MPI_MASTER)call stop_timer(unit=LOGfile)
       ed_file_suffix=""
       !
       call AllReduce_MPI(MpiComm,Smats_tmp,Smats_ineq)
       call AllReduce_MPI(MpiComm,Sreal_tmp,Sreal_ineq)
       call AllReduce_MPI(MpiComm,Gmats_tmp,Gmats_ineq)
       call AllReduce_MPI(MpiComm,Greal_tmp,Greal_ineq)
       call AllReduce_MPI(MpiComm,Dmats_tmp,Dmats_ph_ineq)
       call AllReduce_MPI(MpiComm,Dreal_tmp,Dreal_ph_ineq)
       call AllReduce_MPI(MpiComm,dens_tmp,dens_ineq)
       call AllReduce_MPI(MpiComm,docc_tmp,docc_ineq)
       call AllReduce_MPI(MpiComm,mag_tmp,mag_ineq)
       call AllReduce_MPI(MpiComm,e_tmp,e_ineq)
       call AllReduce_MPI(MpiComm,dd_tmp,dd_ineq)
       call AllReduce_MPI(MpiComm,imp_density_matrix_tmp,imp_density_matrix_ineq)
       call AllReduce_MPI(MpiComm,neigen_sectortmp,neigen_sector_ineq)
       call AllReduce_MPI(MpiComm,neigen_totaltmp,neigen_total_ineq)
       !
       !       
    case(.true.)                !solve sites serial, Lanczos with MPI
       if(MPI_MASTER)call start_timer
       do ilat = 1, Nineq
          write(LOGfile,*)" SOLVES INEQ SITE: "//str(ilat,Npad=4)
          ed_file_suffix=reg(ineq_site_suffix)//str(ilat,site_indx_padding)
          !
          !If required set the local value of U per each site
          if(present(Uloc_ii))Uloc(1:Norb) = Uloc_ii(ilat,1:Norb)
          if(present(Ust_ii)) Ust = Ust_ii(ilat)
          if(present(Jh_ii))  Jh  = Jh_ii(ilat)
          if(present(Jp_ii))  Jp  = Jp_ii(ilat)
          if(present(Jx_ii))  Jx  = Jx_ii(ilat)
          !
          !Solve the impurity problem for the ilat-th site
          neigen_sector(:)   = neigen_sector_ineq(ilat,:)
          lanc_nstates_total = neigen_total_ineq(ilat)
          !
          !Call ed_solve in SERIAL MODE!! This is parallel on the ineq. sites
          call ed_solve_single_mpi(MpiComm,bath(ilat,:),Hloc(ilat,:,:,:,:))
          !
          neigen_sector_ineq(ilat,:)  = neigen_sector(:)
          neigen_total_ineq(ilat)     = lanc_nstates_total
          Smats_ineq(ilat,:,:,:,:,:)  = impSmats(:,:,:,:,:)
          Sreal_ineq(ilat,:,:,:,:,:)  = impSreal(:,:,:,:,:)
          Gmats_ineq(ilat,:,:,:,:,:)  = impGmats(:,:,:,:,:)
          Greal_ineq(ilat,:,:,:,:,:)  = impGreal(:,:,:,:,:)
          Dmats_ph_ineq(ilat,:)       = impDmats_ph(:)
          Dreal_ph_ineq(ilat,:)       = impDreal_ph(:)
          dens_ineq(ilat,1:Norb)      = ed_dens(1:Norb)
          docc_ineq(ilat,1:Norb)      = ed_docc(1:Norb)
          mag_ineq(ilat,1:Norb)       = ed_mag(1:Norb)
          e_ineq(ilat,:)              = [ed_Epot,ed_Eint,ed_Ehartree,ed_Eknot]
          dd_ineq(ilat,:)             = [ed_Dust,ed_Dund,ed_Dse,ed_Dph]
          imp_density_matrix_ineq(ilat,:,:,:,:) = imp_density_matrix(:,:,:,:)
       enddo
       if(MPI_MASTER)call stop_timer(unit=LOGfile)
       ed_file_suffix=""
    end select
    !
  end subroutine ed_solve_lattice_mpi
#endif






  subroutine ineq_memory_allocation(Nineq)
    integer :: Nineq
    if(allocated(dens_ineq))deallocate(dens_ineq)
    if(allocated(docc_ineq))deallocate(docc_ineq)
    if(allocated(mag_ineq))deallocate(mag_ineq)
    if(allocated(e_ineq))deallocate(e_ineq)
    if(allocated(dd_ineq))deallocate(dd_ineq)
    if(allocated(Smats_ineq))deallocate(Smats_ineq)
    if(allocated(Sreal_ineq))deallocate(Sreal_ineq)
    if(allocated(Gmats_ineq))deallocate(Gmats_ineq)
    if(allocated(Greal_ineq))deallocate(Greal_ineq)
    if(allocated(Dmats_ph_ineq))deallocate(Dmats_ph_ineq)
    if(allocated(Dreal_ph_ineq))deallocate(Dreal_ph_ineq)
    if(allocated(imp_density_matrix_ineq))deallocate(imp_density_matrix_ineq)
    if(allocated(neigen_sector_ineq))deallocate(neigen_sector_ineq)
    if(allocated(neigen_total_ineq))deallocate(neigen_total_ineq)
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
    allocate(Dmats_ph_ineq(Nineq,Lmats))
    allocate(Dreal_ph_ineq(Nineq,Lreal))    
    !
    allocate(imp_density_matrix_ineq(Nineq,Nspin,Nspin,Norb,Norb))
    !
    allocate(neigen_sector_ineq(Nineq,Nsectors))
    allocate(neigen_total_ineq(Nineq))
  end subroutine ineq_memory_allocation

end module ED_MAIN







