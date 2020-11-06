MODULE ED_SETUP
  USE ED_INPUT_VARS
  USE ED_VARS_GLOBAL
  USE ED_AUX_FUNX
  USE SF_TIMER
  USE SF_IOTOOLS, only:free_unit,reg,file_length
#ifdef _MPI
  USE MPI
  USE SF_MPI
#endif
  implicit none
  private

  interface map_allocate
     module procedure :: map_allocate_scalar
     module procedure :: map_allocate_vector
  end interface map_allocate

  interface map_deallocate
     module procedure :: map_deallocate_scalar
     module procedure :: map_deallocate_vector
  end interface map_deallocate



  public :: init_ed_structure
  public :: setup_global
  !
  public :: build_sector
  ! public :: build_sector_
  public :: delete_sector
  ! public :: delete_sector_
  ! public :: map_allocate
  ! public :: map_deallocate
  !
  public :: apply_op_C
  public :: apply_op_CDG
  public :: apply_op_Sz
  public :: apply_op_N

  public :: get_Sector
  public :: get_QuantumNumbers
  public :: get_Nup
  public :: get_Ndw
  public :: get_DimUp
  public :: get_DimDw
  !
  public :: indices2state
  public :: state2indices
  public :: iup_index
  public :: idw_index
  !
  public :: bdecomp
  public :: breorder
  public :: bjoin
  !
  public :: c,cdg
  !    
  public :: twin_sector_order
  public :: get_twin_sector
  public :: flip_state
  !
  public :: binary_search

#ifdef _MPI
  public :: scatter_vector_MPI
  public :: scatter_basis_MPI
  public :: gather_vector_MPI
  public :: allgather_vector_MPI
#endif



  interface print_state_vector
     module procedure print_state_vector_ivec
     module procedure print_state_vector_ivec_ud
     module procedure print_state_vector_int
  end interface print_state_vector
  public :: print_state_vector

contains

  subroutine ed_checks_global
    !
    if(Lfit>Lmats)Lfit=Lmats
    if(Nspin>2)stop "ED ERROR: Nspin > 2 is currently not supported"
    if(Norb>5)stop "ED ERROR: Norb > 5 is currently not supported"
    !
    if(.not.ed_total_ud)then
       if(bath_type=="hybrid")stop "ED ERROR: ed_total_ud=F can not be used with bath_type=hybrid"
       if(Jhflag)stop "ED ERROR: ed_total_ud=F can not be used with Jx!=0 OR Jp!=0"
       if(ph_type==2)stop "ED_ERROR: ed_total_ud=F can not be used with ph_type=2"
       ! !<ACTHUNG:
       ! lanc_dim_threshold=2
    endif
    !
    if(Nspin>1.AND.ed_twin.eqv..true.)then
       write(LOGfile,"(A)")"WARNING: using twin_sector with Nspin>1"
       call sleep(1)
    end if
    !
    if(lanc_method=="lanczos")then
       if(lanc_nstates_total>1)stop "ED ERROR: lanc_method==lanczos available only for lanc_nstates_total==1, T=0"
       if(lanc_nstates_sector>1)stop "ED ERROR: lanc_method==lanczos available only for lanc_nstates_sector==1, T=0"
    endif
    !
    if(lanc_method=="dvdson".AND.MpiStatus)then
       if(mpiSIZE>1)stop "ED ERROR: lanc_method=Dvdson + MPIsize>1: not possible at the moment"       
    endif
    !
    if(ed_diag_type=="full".AND.MpiStatus)then
       if(mpiSIZE>1)stop "ED ERROR: ed_diag_type=FULL + MPIsize>1: not possible at the moment"
    end if
    !
    if(ed_diag_type=="lanc")then
       if(ed_finite_temp)then
          if(lanc_nstates_total==1)stop "ED ERROR: ed_diag_type==lanc + ed_finite_temp=T *but* lanc_nstates_total==1 => T=0. Increase lanc_nstates_total"
       else
          if(lanc_nstates_total>1)print*, "ED WARNING: ed_diag_type==lanc + ed_finite_temp=F, T=0 *AND* lanc_nstates_total>1. re-Set lanc_nstates_total=1"
          lanc_nstates_total=1
       endif
    endif
    !
  end subroutine ed_checks_global


  !+------------------------------------------------------------------+
  !PURPOSE  : Setup Dimensions of the problem
  ! Norb    = # of impurity orbitals
  ! Nbath   = # of bath levels (depending on bath_type)
  ! Ns      = # of levels (per spin)
  ! Nlevels = 2*Ns = Total # of levels (counting spin degeneracy 2) 
  !+------------------------------------------------------------------+
  subroutine ed_setup_dimensions()
    select case(bath_type)
    case default
       Ns = (Nbath+1)*Norb
    case ('hybrid')
       Ns = Nbath+Norb
       if(.not.ed_total_ud)stop "ed_setup_dimension: bath_type==hybrid AND .NOT.ed_total_ud"
    case ('replica')
       Ns = Norb*(Nbath+1)
    end select
    !
    select case(ed_total_ud)
    case (.true.)
       Ns_Orb = Ns
       Ns_Ud  = 1
    case (.false.)
       Ns_Orb = Ns/Norb
       Ns_Ud  = Norb
    end select
    !
    DimPh = Nph+1
    Nsectors = ((Ns_Orb+1)*(Ns_Orb+1))**Ns_Ud
  end subroutine ed_setup_dimensions



  !+------------------------------------------------------------------+
  !PURPOSE  : Init ED structure and calculation
  !+------------------------------------------------------------------+
  subroutine init_ed_structure()
    logical                          :: control
    integer                          :: i,iud,iorb,jorb,ispin,jspin
    integer,dimension(:),allocatable :: DimUps,DimDws
    !
    Jhflag=.FALSE.
    if(Norb>1.AND.(Jx/=0d0.OR.Jp/=0d0))Jhflag=.TRUE.
    !
    call ed_checks_global
    !
    call ed_setup_dimensions
    !
    !
    allocate(DimUps(Ns_Ud))
    allocate(DimDws(Ns_Ud))
    do iud=1,Ns_Ud
       DimUps(iud) = get_sector_dimension(Ns_Orb,Ns_Orb/2)
       DimDws(iud) = get_sector_dimension(Ns_Orb,Ns_Orb-Ns_Orb/2)
    enddo
    if(MpiMaster)then
       write(LOGfile,"(A)")"Summary:"
       write(LOGfile,"(A)")"--------------------------------------------"
       write(LOGfile,"(A,I15)")'# of levels/spin      = ',Ns
       write(LOGfile,"(A,I15)")'Total size            = ',2*Ns
       write(LOGfile,"(A,I15)")'# of impurities       = ',Norb
       write(LOGfile,"(A,I15)")'# of bath/impurity    = ',Nbath
       write(LOGfile,"(A,I15)")'# of Bath levels/spin = ',Ns-Norb
       write(LOGfile,"(A,I15)")'# of  sectors         = ',Nsectors
       write(LOGfile,"(A,I15)")'Ns_Orb                = ',Ns_Orb
       write(LOGfile,"(A,I15)")'Ns_Ud                 = ',Ns_Ud
       write(LOGfile,"(A,I15)")'Nph                   = ',Nph
       write(LOGfile,"(A,"//str(Ns_Ud)//"I8,2X,"//str(Ns_Ud)//"I8,I8,I20)")&
            'Largest Sector(s)     = ',DimUps,DimDws,DimPh,product(DimUps)*product(DimDws)*DimPh
       write(LOGfile,"(A)")"--------------------------------------------"
    endif
    call sleep(1)
    !
    ! allocate(impHloc(Nspin,Nspin,Norb,Norb))
    ! impHloc=zero
    !
    allocate(spH0ups(Ns_Ud))
    allocate(spH0dws(Ns_Ud))
    !
    !Allocate indexing arrays
    allocate(getCsector(Ns_Ud,2,Nsectors))  ;getCsector  =0
    allocate(getCDGsector(Ns_Ud,2,Nsectors));getCDGsector=0
    !
    allocate(getDim(Nsectors));getDim=0
    !
    allocate(getBathStride(Norb,Nbath));getBathStride=0
    allocate(twin_mask(Nsectors))
    allocate(sectors_mask(Nsectors))
    allocate(neigen_sector(Nsectors))
    !
    !
    finiteT = ed_finite_temp
    !
    if(ed_diag_type=="lanc")then
       if(finiteT)then
          if(mod(lanc_nstates_sector,2)/=0)then
             lanc_nstates_sector=lanc_nstates_sector+1
             write(LOGfile,"(A,I10)")"Increased Lanc_nstates_sector:",lanc_nstates_sector
          endif
          if(mod(lanc_nstates_total,2)/=0)then
             lanc_nstates_total=lanc_nstates_total+1
             write(LOGfile,"(A,I10)")"Increased Lanc_nstates_total:",lanc_nstates_total
          endif
          write(LOGfile,"(A,I3)")"Nstates x Sector = ", lanc_nstates_sector
          write(LOGfile,"(A,I3)")"Nstates   Total  = ", lanc_nstates_total
          !
          write(LOGfile,"(A)")"Lanczos FINITE temperature calculation:"
          call sleep(1)
       else
          write(LOGfile,"(A)")"Lanczos ZERO temperature calculation:"
          call sleep(1)
       endif
    else
       if(finiteT)then
          write(LOGfile,"(A)")"Full ED finite T calculation"
          call sleep(1)
       else
          ed_diag_type='lanc'
          lanc_nstates_total=1
          lanc_dim_threshold=product(DimUps)*product(DimDws)*Dimph
          write(LOGfile,"(A)")"Full ED T=0 calculation. Set LANC_DIM_THRESHOLD to "//str(lanc_dim_threshold)
          if(lanc_dim_threshold>2**13)stop "Full ED T=0: LANC_DIM_THRESHOLD > 2**13=8192!"
          call sleep(1)
       endif
    endif
    !
    !
    offdiag_gf_flag=ed_solve_offdiag_gf
    if(bath_type/="normal")offdiag_gf_flag=.true.
    if(.not.ed_total_ud.AND.offdiag_gf_flag)then
       write(LOGfile,"(A)")"ED WARNING: can not do offdiag_gf_flag=T.AND.ed_total_ud=F. Set to F."
       offdiag_gf_flag=.false.
    endif
    !
    !
    if(nread/=0.d0)then
       i=abs(floor(log10(abs(nerr)))) !modulus of the order of magnitude of nerror
       niter=nloop/3
    endif
    !
    !
    !allocate functions
    allocate(impSmats(Nspin,Nspin,Norb,Norb,Lmats))
    allocate(impSreal(Nspin,Nspin,Norb,Norb,Lreal))
    impSmats=zero
    impSreal=zero
    !
    allocate(impGmats(Nspin,Nspin,Norb,Norb,Lmats))
    allocate(impGreal(Nspin,Nspin,Norb,Norb,Lreal))
    impGmats=zero
    impGreal=zero
    !
    allocate(impG0mats(Nspin,Nspin,Norb,Norb,Lmats))
    allocate(impG0real(Nspin,Nspin,Norb,Norb,Lreal))
    impG0mats=zero
    impG0real=zero
    !
    allocate(impDmats_ph(0:Lmats))
    allocate(impDreal_ph(Lreal))
    impDmats_ph=zero
    impDreal_ph=zero
    !
    !allocate observables
    allocate(ed_dens(Norb),ed_docc(Norb),ed_dens_up(Norb),ed_dens_dw(Norb))
    ed_dens=0d0
    ed_docc=0d0
    ed_dens_up=0d0
    ed_dens_dw=0d0
    !
    allocate(spinChi_tau(Norb,Norb,0:Ltau))
    allocate(spinChi_w(Norb,Norb,Lreal))
    allocate(spinChi_iv(Norb,Norb,0:Lmats))
    !
    allocate(densChi_tau(Norb,Norb,0:Ltau))
    allocate(densChi_w(Norb,Norb,Lreal))
    allocate(densChi_iv(Norb,Norb,0:Lmats))
    !
    allocate(pairChi_tau(Norb,Norb,0:Ltau))
    allocate(pairChi_w(Norb,Norb,Lreal))
    allocate(pairChi_iv(Norb,Norb,0:Lmats))
    !
    allocate(exctChi_tau(0:2,Norb,Norb,0:Ltau))
    allocate(exctChi_w(0:2,Norb,Norb,Lreal))
    allocate(exctChi_iv(0:2,Norb,Norb,0:Lmats))
    !
  end subroutine init_ed_structure





  !+------------------------------------------------------------------+
  !PURPOSE: SETUP THE GLOBAL POINTERS FOR THE ED CALCULAIONS.
  !+------------------------------------------------------------------+
  subroutine setup_global
    integer                          :: DimUp,DimDw
    integer                          :: DimUps(Ns_Ud),DimDws(Ns_Ud)
    integer                          :: Indices(2*Ns_Ud),Jndices(2*Ns_Ud)
    integer                          :: Nups(Ns_ud),Ndws(Ns_ud)
    integer                          :: Jups(Ns_ud),Jdws(Ns_ud)
    integer                          :: i,iud,iorb
    integer                          :: isector,jsector,gsector,ksector,lsector
    integer                          :: unit,status,istate,ishift,isign
    logical                          :: IOfile
    integer                          :: list_len
    integer,dimension(:),allocatable :: list_sector
    type(sector) :: sectorI,sectorJ,sectorK,sectorG,sectorL
    !
    !Store full dimension of the sectors:
    do isector=1,Nsectors
       call get_DimUp(isector,DimUps)
       call get_DimDw(isector,DimDws)
       DimUp = product(DimUps)
       DimDw = product(DimDws)  
       getDim(isector)  = DimUp*DimDw*DimPh
    enddo
    !
    !
    inquire(file="state_list"//reg(ed_file_suffix)//".restart",exist=IOfile)
    if(IOfile)then
       write(LOGfile,"(A)")"Restarting from a state_list file:"
       list_len=file_length("state_list"//reg(ed_file_suffix)//".restart")
       allocate(list_sector(list_len))
       !
       open(free_unit(unit),file="state_list"//reg(ed_file_suffix)//".restart",status="old")
       status=0
       do while(status>=0)
          read(unit,*,iostat=status)istate,isector,indices
          list_sector(istate)=isector
          call get_Nup(isector,Nups)
          call get_Ndw(isector,Ndws)
          if(any(Indices /= [Nups,Ndws]))&
               stop "setup_global error: nups!=nups(isector).OR.ndws!=ndws(isector)"
       enddo
       close(unit)
       !
       lanc_nstates_total = list_len
       do isector=1,Nsectors
          neigen_sector(isector) = max(1,count(list_sector==isector))
       enddo
    else
       do isector=1,Nsectors
          neigen_sector(isector) = min(getDim(isector),lanc_nstates_sector) !init every sector to required eigenstates
       enddo
    endif
    !
    twin_mask=.true.
    if(ed_twin)then
       do isector=1,Nsectors
          call get_Nup(isector,Nups)
          call get_Ndw(isector,Ndws)
          if(any(Nups < Ndws))twin_mask(isector)=.false.
       enddo
       write(LOGfile,"(A,I6,A,I9)")"Looking into ",count(twin_mask)," sectors out of ",Nsectors
       call sleep(1)
    endif
    !
    select case(bath_type)
    case default
       do i=1,Nbath
          do iorb=1,Norb
             getBathStride(iorb,i) = Norb + (iorb-1)*Nbath + i
          enddo
       enddo
    case ('hybrid')
       do i=1,Nbath
          getBathStride(:,i)       = Norb + i
       enddo
    case ('replica')
       do i=1,Nbath
          do iorb=1,Norb
             getBathStride(iorb,i) = iorb + i*Norb 
          enddo
       enddo
    end select
    !
    getCsector  = 0
    getCDGsector= 0
    do isector=1,Nsectors
       call get_Nup(isector,Nups)
       call get_Ndw(isector,Ndws)
       !
       !UPs:
       !DEL:
       do iud=1,Ns_Ud
          Jups=Nups
          Jdws=Ndws 
          Jups(iud)=Jups(iud)-1; if(Jups(iud) < 0)cycle
          call get_Sector([Jups,Jdws],Ns_Orb,jsector)
          getCsector(iud,1,isector)=jsector
       enddo
       !ADD
       do iud=1,Ns_Ud
          Jups=Nups
          Jdws=Ndws
          Jups(iud)=Jups(iud)+1; if(Jups(iud) > Ns_Orb)cycle
          call get_Sector([Jups,Jdws],Ns_Orb,jsector)
          getCDGsector(iud,1,isector)=jsector
       enddo
       !
       !DWs:
       !DEL
       do iud=1,Ns_Ud
          Jups=Nups
          Jdws=Ndws 
          Jdws(iud)=Jdws(iud)-1; if(Jdws(iud) < 0)cycle
          call get_Sector([Jups,Jdws],Ns_Orb,jsector)
          getCsector(iud,2,isector)=jsector
       enddo
       !DEL
       do iud=1,Ns_Ud
          Jups=Nups
          Jdws=Ndws 
          Jdws(iud)=Jdws(iud)+1; if(Jdws(iud) > Ns_Orb)cycle
          call get_Sector([Jups,Jdws],Ns_Orb,jsector)
          getCDGsector(iud,2,isector)=jsector
       enddo
    enddo
    return
  end subroutine setup_global













  !##################################################################
  !##################################################################
  !AUXILIARY PROCEDURES - Sectors,Nup,Ndw,DimUp,DimDw,...
  !##################################################################
  !##################################################################
  !+------------------------------------------------------------------+
  !PURPOSE  : input a state |i> and output a vector ivec(Nlevels)
  !with its binary decomposition
  !(corresponds to the decomposition of the number i-1)
  !+------------------------------------------------------------------+
  function bdecomp(i,Ntot) result(ivec)
    integer :: Ntot,ivec(Ntot),l,i
    logical :: busy
    !this is the configuration vector |1,..,Ns,Ns+1,...,Ntot>
    !obtained from binary decomposition of the state/number i\in 2^Ntot
    do l=0,Ntot-1
       busy=btest(i,l)
       ivec(l+1)=0
       if(busy)ivec(l+1)=1
    enddo
  end function bdecomp


  !+------------------------------------------------------------------+
  ! Reorder a binary decomposition so to have a state of the form:
  ! default: |(1:Norb),([1:Nbath]_1, [1:Nbath]_2, ... ,[1:Nbath]_Norb)>_spin
  ! hybrid:  |(1:Norb),([1:Nbath])_spin
  ! replica: |(1:Norb),([1:Norb]_1, [1:Norb]_2, ...  , [1:Norb]_Nbath)>_spin
  !
  !> case (ed_total_ud):
  !   (T): Ns_Ud=1, Ns_Orb=Ns.
  !        bdecomp is already of the form above [1:Ns]
  !   (F): Ns_Ud=Norb, Ns_Orb=Ns/Norb==1+Nbath
  !        bdecomp is
  !        |( [1:1+Nbath]_1,...,[1:1+Nbath]_Norb)>_spin
  !+------------------------------------------------------------------+
  function breorder(Nins) result(Ivec)
    integer,intent(in),dimension(Ns_Ud,Ns_Orb) :: Nins ![1,Ns] - [Norb,1+Nbath]
    integer,dimension(Ns)                      :: Ivec ![Ns]
    integer                                    :: iud,ibath,indx
    select case (ed_total_ud)
    case (.true.)
       Ivec = Nins(1,:)
    case (.false.)
       do iud=1,Ns_Ud           ![1:Norb]
          Ivec(iud) = Nins(iud,1)
          do ibath=1,Nbath
             indx = getBathStride(iud,ibath) !take care of normal/
             Ivec(indx) = Nins(iud,1+ibath)
          enddo
       enddo
    end select
  end function breorder


  !+------------------------------------------------------------------+
  !PURPOSE  : input a vector ib(Nlevels) with the binary sequence 
  ! and output the corresponding state |i>
  !(corresponds to the recomposition of the number i-1)
  !+------------------------------------------------------------------+
  function bjoin(ib,Ntot) result(i)
    integer                 :: Ntot
    integer,dimension(Ntot) :: ib
    integer                 :: i,j
    i=0
    do j=0,Ntot-1
       i=i+ib(j+1)*2**j
    enddo
  end function bjoin


  elemental function get_sector_dimension(n,np) result(dim)
    integer,intent(in) :: n,np
    integer            :: dim
    dim = binomial(n,np)
  end function get_sector_dimension


  subroutine get_Sector(QN,N,isector)
    integer,dimension(:) :: QN
    integer              :: N
    integer              :: isector
    integer              :: i,Nind,factor
    Nind = size(QN)
    Factor = N+1
    isector = 1
    do i=Nind,1,-1
       isector = isector + QN(i)*(Factor)**(Nind-i)
    enddo
  end subroutine get_Sector


  subroutine get_QuantumNumbers(isector,N,QN)
    integer                          :: isector,N
    integer,dimension(:)             :: QN
    integer                          :: i,count,Dim
    integer,dimension(size(QN)) :: QN_
    !
    Dim = size(QN)
    if(mod(Dim,2)/=0)stop "get_QuantumNumbers error: Dim%2 != 0"
    count=isector-1
    do i=1,Dim
       QN_(i) = mod(count,N+1)
       count      = count/(N+1)
    enddo
    QN = QN_(Dim:1:-1)
  end subroutine get_QuantumNumbers


  subroutine get_Nup(isector,Nup)
    integer                   :: isector,Nup(Ns_Ud)
    integer                   :: i,count
    integer,dimension(2*Ns_Ud)  :: indices_
    count=isector-1
    do i=1,2*Ns_Ud
       indices_(i) = mod(count,Ns_Orb+1)
       count      = count/(Ns_Orb+1)
    enddo
    Nup = indices_(2*Ns_Ud:Ns_Ud+1:-1)
  end subroutine get_Nup


  subroutine get_Ndw(isector,Ndw)
    integer                   :: isector,Ndw(Ns_Ud)
    integer                   :: i,count
    integer,dimension(2*Ns_Ud) :: indices_
    count=isector-1
    do i=1,2*Ns_Ud
       indices_(i) = mod(count,Ns_Orb+1)
       count      = count/(Ns_Orb+1)
    enddo
    Ndw = indices_(Ns_Ud:1:-1)
  end subroutine get_Ndw


  subroutine  get_DimUp(isector,DimUps)
    integer                :: isector,DimUps(Ns_Ud)
    integer                :: Nups(Ns_Ud),iud
    call get_Nup(isector,Nups)
    do iud=1,Ns_Ud
       DimUps(iud) = binomial(Ns_Orb,Nups(iud))
    enddo
  end subroutine get_DimUp


  subroutine get_DimDw(isector,DimDws)
    integer                :: isector,DimDws(Ns_Ud)
    integer                :: Ndws(Ns_Ud),iud
    call get_Ndw(isector,Ndws)
    do iud=1,Ns_Ud
       DimDws(iud) = binomial(Ns_Orb,Ndws(iud))
    enddo
  end subroutine get_DimDw


  subroutine indices2state(ivec,Nvec,istate)
    integer,dimension(:)          :: ivec
    integer,dimension(size(ivec)) :: Nvec
    integer                       :: istate,i
    istate=ivec(1)
    do i=2,size(ivec)
       istate = istate + (ivec(i)-1)*product(Nvec(1:i-1))
    enddo
  end subroutine indices2state

  subroutine state2indices(istate,Nvec,ivec)
    integer                       :: istate
    integer,dimension(:)          :: Nvec
    integer,dimension(size(Nvec)) :: Ivec
    integer                       :: i,count,N
    count = istate-1
    N     = size(Nvec)
    do i=1,N
       Ivec(i) = mod(count,Nvec(i))+1
       count   = count/Nvec(i)
    enddo
  end subroutine state2indices


  function iup_index(i,DimUp) result(iup)
    integer :: i
    integer :: DimUp
    integer :: iup
    iup = mod(i,DimUp);if(iup==0)iup=DimUp
  end function iup_index


  function idw_index(i,DimUp) result(idw)
    integer :: i
    integer :: DimUp
    integer :: idw
    idw = (i-1)/DimUp+1
  end function idw_index





  !##################################################################
  !##################################################################
  !BUILD SECTORS
  !##################################################################
  !##################################################################
  subroutine build_sector(isector,self)
    integer,intent(in)                  :: isector
    type(sector)                        :: self
    integer                             :: iup,idw
    integer                             :: nup_,ndw_
    integer                             :: dim,iud
    !
    if(self%status)call delete_sector(self)
    !
    self%index = isector
    !
    allocate(self%H(2*Ns_Ud))
    allocate(self%DimUps(Ns_Ud))
    allocate(self%DimDws(Ns_Ud))
    allocate(self%Nups(Ns_Ud))
    allocate(self%Ndws(Ns_Ud))
    !
    call get_Nup(isector,self%Nups);self%Nup=sum(self%Nups)
    call get_Ndw(isector,self%Ndws);self%Ndw=sum(self%Ndws)
    call get_DimUp(isector,self%DimUps);self%DimUp=product(self%DimUps)
    call get_DimDw(isector,self%DimDws);self%DimDw=product(self%DimDws)
    self%DimEl=self%DimUp*self%DimDw
    self%DimPh=Nph+1
    self%Dim=self%DimEl*self%DimPh
    !
    call map_allocate(self%H,[self%DimUps,self%DimDws])
    do iud=1,Ns_Ud
       !UP    
       dim=0
       do iup=0,2**Ns_Orb-1
          nup_ = popcnt(iup)
          if(nup_ /= self%Nups(iud))cycle
          dim  = dim+1
          self%H(iud)%map(dim) = iup
       enddo
       !DW
       dim=0
       do idw=0,2**Ns_Orb-1
          ndw_= popcnt(idw)
          if(ndw_ /= self%Ndws(iud))cycle
          dim = dim+1
          self%H(iud+Ns_Ud)%map(dim) = idw
       enddo
    enddo
    !
    self%Nlanc = min(self%Dim,lanc_nGFiter)
    self%status=.true.
  end subroutine build_sector


  subroutine delete_sector(self)
    type(sector) :: self
    call map_deallocate(self%H)
    if(allocated(self%H))deallocate(self%H)
    if(allocated(self%DimUps))deallocate(self%DimUps)
    if(allocated(self%DimDws))deallocate(self%DimDws)
    if(allocated(self%Nups))deallocate(self%Nups)
    if(allocated(self%Ndws))deallocate(self%Ndws)
    self%index=0
    self%DimUp=0
    self%DimDw=0
    self%Dim=0
    self%Nup=0
    self%Ndw=0
    self%Nlanc=0
    self%status=.false.
  end subroutine delete_sector


  subroutine map_allocate_scalar(H,N)
    type(sector_map) :: H
    integer          :: N
    if(H%status) call map_deallocate_scalar(H)
    allocate(H%map(N))
    H%status=.true.
  end subroutine map_allocate_scalar
  !
  subroutine map_allocate_vector(H,N)
    type(sector_map),dimension(:)       :: H
    integer,dimension(size(H))          :: N
    integer                             :: i
    do i=1,size(H)
       call map_allocate_scalar(H(i),N(i))
    enddo
  end subroutine map_allocate_vector



  subroutine map_deallocate_scalar(H)
    type(sector_map) :: H
    if(.not.H%status)then
       write(*,*) "WARNING map_deallocate_scalar: H is not allocated"
       return
    endif
    if(allocated(H%map))deallocate(H%map)
    H%status=.false.
  end subroutine map_deallocate_scalar
  !
  subroutine map_deallocate_vector(H)
    type(sector_map),dimension(:) :: H
    integer                       :: i
    do i=1,size(H)
       call map_deallocate_scalar(H(i))
    enddo
  end subroutine map_deallocate_vector




  subroutine apply_op_C(i,j,sgn,ipos,ialfa,ispin,sectorI,sectorJ) 
    integer, intent(in)         :: i,ipos,ialfa,ispin
    type(sector),intent(in)     :: sectorI,sectorJ
    integer,intent(out)         :: j
    real(8),intent(out)         :: sgn
    integer                     :: ibeta
    integer                     :: r
    integer                     :: iph,i_el
    integer,dimension(2*Ns_Ud)  :: Indices
    integer,dimension(2*Ns_Ud)  :: Jndices
    integer,dimension(2,Ns_Orb) :: Nud !Nbits(Ns_Orb)
    integer,dimension(2)        :: Iud
    !
    j=0
    sgn=0d0
    !
    ibeta  = ialfa + (ispin-1)*Ns_Ud
    iph = (i-1)/(sectorI%DimEl) + 1
    i_el = mod(i-1,sectorI%DimEl) + 1
    !
    call state2indices(i_el,[sectorI%DimUps,sectorI%DimDws],Indices)
    iud(1)   = sectorI%H(ialfa)%map(Indices(ialfa))
    iud(2)   = sectorI%H(ialfa+Ns_Ud)%map(Indices(ialfa+Ns_Ud))
    nud(1,:) = Bdecomp(iud(1),Ns_Orb)
    nud(2,:) = Bdecomp(iud(2),Ns_Orb)
    if(Nud(ispin,ipos)/=1)return
    call c(ipos,iud(ispin),r,sgn)
    Jndices        = Indices
    Jndices(ibeta) = binary_search(sectorJ%H(ibeta)%map,r)
    call indices2state(Jndices,[sectorJ%DimUps,sectorJ%DimDws],j)
    !
    j = j + (iph-1)*sectorJ%DimEl
  end subroutine apply_op_C


  subroutine apply_op_CDG(i,j,sgn,ipos,ialfa,ispin,sectorI,sectorJ) 
    integer, intent(in)         :: i,ipos,ialfa,ispin
    type(sector),intent(in)     :: sectorI,sectorJ
    integer,intent(out)         :: j
    real(8),intent(out)         :: sgn
    integer                     :: ibeta
    integer                     :: r
    integer                     :: iph,i_el
    integer,dimension(2*Ns_Ud)  :: Indices
    integer,dimension(2*Ns_Ud)  :: Jndices
    integer,dimension(2,Ns_Orb) :: Nud !Nbits(Ns_Orb)
    integer,dimension(2)        :: Iud
    !
    j=0
    sgn=0d0
    !
    ibeta  = ialfa + (ispin-1)*Ns_Ud
    iph = (i-1)/(sectorI%DimEl) + 1
    i_el = mod(i-1,sectorI%DimEl) + 1
    !
    call state2indices(i_el,[sectorI%DimUps,sectorI%DimDws],Indices)
    iud(1)   = sectorI%H(ialfa)%map(Indices(ialfa))
    iud(2)   = sectorI%H(ialfa+Ns_Ud)%map(Indices(ialfa+Ns_Ud))
    nud(1,:) = Bdecomp(iud(1),Ns_Orb)
    nud(2,:) = Bdecomp(iud(2),Ns_Orb)
    if(Nud(ispin,ipos)/=0)return
    call cdg(ipos,iud(ispin),r,sgn)
    Jndices        = Indices
    Jndices(ibeta) = binary_search(sectorJ%H(ibeta)%map,r)
    call indices2state(Jndices,[sectorJ%DimUps,sectorJ%DimDws],j)
    !
    j = j + (iph-1)*sectorJ%DimEl
  end subroutine apply_op_CDG


  subroutine apply_op_Sz(i,sgn,ipos,ialfa,sectorI) 
    integer, intent(in)         :: i,ipos,ialfa
    type(sector),intent(in)     :: sectorI
    real(8),intent(out)         :: sgn
    integer                     :: iph,i_el
    integer,dimension(2*Ns_Ud)  :: Indices
    integer,dimension(2*Ns_Ud)  :: Jndices
    integer,dimension(2,Ns_Orb) :: Nud !Nbits(Ns_Orb)
    integer,dimension(2)        :: Iud
    !
    sgn=0d0
    !
    iph = (i-1)/(sectorI%DimEl) + 1
    i_el = mod(i-1,sectorI%DimEl) + 1
    !
    call state2indices(i_el,[sectorI%DimUps,sectorI%DimDws],Indices)
    iud(1)   = sectorI%H(ialfa)%map(Indices(ialfa))
    iud(2)   = sectorI%H(ialfa+Ns_Ud)%map(Indices(ialfa+Ns_Ud))
    nud(1,:) = Bdecomp(iud(1),Ns_Orb)
    nud(2,:) = Bdecomp(iud(2),Ns_Orb)
    !
    sgn = dble(nud(1,ipos))-dble(nud(2,ipos))
    sgn = sgn/2d0
  end subroutine apply_op_Sz



  subroutine apply_op_N(i,sgn,ipos,ialfa,sectorI) 
    integer, intent(in)         :: i,ipos,ialfa
    type(sector),intent(in)     :: sectorI
    real(8),intent(out)         :: sgn
    integer                     :: iph,i_el
    integer,dimension(2*Ns_Ud)  :: Indices
    integer,dimension(2*Ns_Ud)  :: Jndices
    integer,dimension(2,Ns_Orb) :: Nud !Nbits(Ns_Orb)
    integer,dimension(2)        :: Iud
    !
    sgn=0d0
    !
    iph = (i-1)/(sectorI%DimEl) + 1
    i_el = mod(i-1,sectorI%DimEl) + 1
    !
    call state2indices(i_el,[sectorI%DimUps,sectorI%DimDws],Indices)
    iud(1)   = sectorI%H(ialfa)%map(Indices(ialfa))
    iud(2)   = sectorI%H(ialfa+Ns_Ud)%map(Indices(ialfa+Ns_Ud))
    nud(1,:) = Bdecomp(iud(1),Ns_Orb)
    nud(2,:) = Bdecomp(iud(2),Ns_Orb)
    !
    sgn = dble(nud(1,ipos))+dble(nud(2,ipos))
  end subroutine apply_op_N



  !##################################################################
  !##################################################################
  !CREATION / DESTRUCTION OPERATORS
  !##################################################################
  !##################################################################
  !+-------------------------------------------------------------------+
  !PURPOSE: input state |in> of the basis and calculates 
  !   |out>=C_pos|in>  OR  |out>=C^+_pos|in> ; 
  !   the sign of |out> has the phase convention, pos labels the sites
  !+-------------------------------------------------------------------+
  subroutine c(pos,in,out,fsgn)
    integer,intent(in)    :: pos
    integer,intent(in)    :: in
    integer,intent(inout) :: out
    real(8),intent(inout) :: fsgn    
    integer               :: l
    if(.not.btest(in,pos-1))stop "C error: C_i|...0_i...>"
    fsgn=1d0
    do l=1,pos-1
       if(btest(in,l-1))fsgn=-fsgn
    enddo
    out = ibclr(in,pos-1)
  end subroutine c

  subroutine cdg(pos,in,out,fsgn)
    integer,intent(in)    :: pos
    integer,intent(in)    :: in
    integer,intent(inout) :: out
    real(8),intent(inout) :: fsgn    
    integer               :: l
    if(btest(in,pos-1))stop "C^+ error: C^+_i|...1_i...>"
    fsgn=1d0
    do l=1,pos-1
       if(btest(in,l-1))fsgn=-fsgn
    enddo
    out = ibset(in,pos-1)
  end subroutine cdg






  !##################################################################
  !##################################################################
  !TWIN SECTORS ROUTINES:
  !##################################################################
  !##################################################################

  !+------------------------------------------------------------------+
  !PURPOSE  : Build the re-ordering map to go from sector A(nup,ndw)
  ! to its twin sector B(ndw,nup), with nup!=ndw.
  !
  !- build the map from the A-sector to \HHH
  !- get the list of states in \HHH corresponding to sector B twin of A
  !- return the ordering of B-states in \HHH with respect to those of A
  !+------------------------------------------------------------------+
  subroutine twin_sector_order(isector,order)
    integer                             :: isector
    integer,dimension(:)                :: order
    type(sector)                        :: sectorH
    type(sector_map),dimension(2*Ns_Ud) :: H
    integer,dimension(2*Ns_Ud)          :: Indices,Istates
    integer,dimension(Ns_Ud)            :: DimUps,DimDws
    integer                             :: Dim,DimUp,DimDw
    integer                             :: i,iud,iph,i_el
    !
    Dim = GetDim(isector)
    if(size(Order)/=Dim)stop "twin_sector_order error: wrong dimensions of *order* array"
    call get_DimUp(isector,DimUps)
    call get_DimDw(isector,DimDws)
    DimUp = product(DimUps)
    DimDw = product(DimDws)
    !
    call build_sector(isector,sectorH)
    do i=1,sectorH%Dim
       iph = (i-1)/(sectorH%DimEl) + 1   !find number of phonons
       i_el = mod(i-1,sectorH%DimEl) + 1 !electronic index
       call state2indices(i_el,[sectorH%DimUps,sectorH%DimDws],Indices)
       forall(iud=1:2*Ns_Ud)Istates(iud) = sectorH%H(iud)%map(Indices(iud))
       Order(i) = flip_state( Istates ) + (iph-1)*2**(2*Ns) !flipped electronic state (GLOBAL state number {1:2^2Ns}) + phononic contribution
    enddo
    call delete_sector(sectorH)
    !
    call sort_array(Order) !sorted and changed the values from the global state numbers to the ones of the sector {1:DimUp*DimDw*DimPh}
    !
  end subroutine twin_sector_order



  !+------------------------------------------------------------------+
  !PURPOSE  : Flip an Hilbert space state m=|{up}>|{dw}> into:
  !
  ! normal: j=|{dw}>|{up}>  , nup --> ndw
  !+------------------------------------------------------------------+
  function flip_state(istate) result(j)
    integer,dimension(2*Ns_Ud) :: istate
    integer                    :: j
    integer,dimension(Ns_Ud)   :: jups,jdws
    integer,dimension(2*Ns_Ud) :: dims
    !
    jups = istate(Ns_Ud+1:2*Ns_Ud)
    jdws = istate(1:Ns_Ud)
    dims = 2**Ns_Orb
    call indices2state([jups,jdws],Dims,j)
    !
  end function flip_state


  !+------------------------------------------------------------------+
  !PURPOSE  : get the twin of a given sector (the one with opposite 
  ! quantum numbers): 
  ! nup,ndw ==> ndw,nup (spin-exchange)
  !+------------------------------------------------------------------+
  function get_twin_sector(isector) result(jsector)
    integer,intent(in)       :: isector
    integer                  :: jsector
    integer,dimension(Ns_Ud) :: Iups,Idws
    call get_Nup(isector,iups)
    call get_Ndw(isector,idws)
    call get_Sector([idws,iups],Ns_Orb,jsector)
  end function get_twin_sector











  !##################################################################
  !##################################################################
  !AUXILIARY COMPUTATIONAL ROUTINES ARE HERE BELOW:
  !##################################################################
  !##################################################################
#ifdef _MPI
  !! Scatter V into the arrays Vloc on each thread: sum_threads(size(Vloc)) must be equal to size(v)
  subroutine scatter_vector_MPI(MpiComm,v,vloc)
    integer                          :: MpiComm
    real(8),dimension(:)             :: v    !size[N]
    real(8),dimension(:)             :: vloc !size[Nloc]
    integer                          :: i,iph,irank,Nloc,N
    integer                          :: v_start,v_end,vloc_start,vloc_end
    integer,dimension(:),allocatable :: Counts,Offset
    integer                          :: MpiSize,MpiIerr
    logical                          :: MpiMaster
    !
    if( MpiComm == MPI_UNDEFINED .OR. MpiComm == Mpi_Comm_Null )return
    ! stop "scatter_vector_MPI error: MpiComm == MPI_UNDEFINED"
    !
    MpiSize   = get_size_MPI(MpiComm)
    MpiMaster = get_master_MPI(MpiComm)
    !
    Nloc = size(Vloc)
    N = 0
    call AllReduce_MPI(MpiComm,Nloc,N)
    if(MpiMaster.AND.N /= size(V)) stop "scatter_vector_MPI error: size(V) != Mpi_Allreduce(Nloc)"
    !
    allocate(Counts(0:MpiSize-1)) ; Counts=0
    allocate(Offset(0:MpiSize-1)) ; Offset=0
    !
    !Get Counts;
    call MPI_AllGather(Nloc/Dimph,1,MPI_INTEGER,Counts,1,MPI_INTEGER,MpiComm,MpiIerr)
    !
    !Get Offset:
    Offset(0)=0
    do i=1,MpiSize-1
       Offset(i) = Offset(i-1) + Counts(i-1)
    enddo
    !
    Vloc=0d0
    do iph=1,Dimph
       if(MpiMaster)then
          v_start = 1 + (iph-1)*(N/Dimph)
          v_end = iph*(N/Dimph)
       else
          v_start = 1
          v_end = 1
       endif
       vloc_start = 1 + (iph-1)*(Nloc/Dimph)
       vloc_end = iph*(Nloc/Dimph)
       call MPI_Scatterv(V(v_start:v_end),Counts,Offset,MPI_DOUBLE_PRECISION,&
            Vloc(vloc_start:vloc_end),Nloc/DimPh,MPI_DOUBLE_PRECISION,0,MpiComm,MpiIerr)
    enddo
    !
    return
  end subroutine scatter_vector_MPI


  subroutine scatter_basis_MPI(MpiComm,v,vloc)
    integer                :: MpiComm
    real(8),dimension(:,:) :: v    !size[N,N]
    real(8),dimension(:,:) :: vloc !size[Nloc,Neigen]
    integer                :: N,Nloc,Neigen,i
    N      = size(v,1)
    Nloc   = size(vloc,1)
    Neigen = size(vloc,2)
    if( size(v,2) < Neigen ) stop "error scatter_basis_MPI: size(v,2) < Neigen"
    !
    do i=1,Neigen
       call scatter_vector_MPI(MpiComm,v(:,i),vloc(:,i))
    end do
    !
    return
  end subroutine scatter_basis_MPI


  !! AllGather Vloc on each thread into the array V: sum_threads(size(Vloc)) must be equal to size(v)
  subroutine gather_vector_MPI(MpiComm,vloc,v)
    integer                          :: MpiComm
    real(8),dimension(:)             :: vloc !size[Nloc]
    real(8),dimension(:)             :: v    !size[N]
    integer                          :: i,iph,irank,Nloc,N
    integer                          :: v_start,v_end,vloc_start,vloc_end
    integer,dimension(:),allocatable :: Counts,Offset
    integer                          :: MpiSize,MpiIerr
    logical                          :: MpiMaster
    !
    if(  MpiComm == MPI_UNDEFINED .OR. MpiComm == Mpi_Comm_Null ) return
    !stop "gather_vector_MPI error: MpiComm == MPI_UNDEFINED"
    !
    MpiSize   = get_size_MPI(MpiComm)
    MpiMaster = get_master_MPI(MpiComm)
    !
    Nloc = size(Vloc)
    N = 0
    call AllReduce_MPI(MpiComm,Nloc,N)
    if(MpiMaster.AND.N /= size(V)) stop "gather_vector_MPI error: size(V) != Mpi_Allreduce(Nloc)"
    !
    allocate(Counts(0:MpiSize-1)) ; Counts=0
    allocate(Offset(0:MpiSize-1)) ; Offset=0
    !
    !Get Counts;
    call MPI_AllGather(Nloc/Dimph,1,MPI_INTEGER,Counts,1,MPI_INTEGER,MpiComm,MpiIerr)
    !
    !Get Offset:
    Offset(0)=0
    do i=1,MpiSize-1
       Offset(i) = Offset(i-1) + Counts(i-1)
    enddo
    !
    do iph=1,Dimph
       if(MpiMaster)then
          v_start = 1 + (iph-1)*(N/Dimph)
          v_end = iph*(N/Dimph)
       else
          v_start = 1
          v_end = 1
       endif
       vloc_start = 1 + (iph-1)*(Nloc/Dimph)
       vloc_end = iph*(Nloc/Dimph)
       !
       call MPI_Gatherv(Vloc(vloc_start:vloc_end),Nloc/DimPh,MPI_DOUBLE_PRECISION,&
            V(v_start:v_end),Counts,Offset,MPI_DOUBLE_PRECISION,0,MpiComm,MpiIerr)
    enddo
    !
    return
  end subroutine gather_vector_MPI


  !! AllGather Vloc on each thread into the array V: sum_threads(size(Vloc)) must be equal to size(v)
  subroutine allgather_vector_MPI(MpiComm,vloc,v)
    integer                          :: MpiComm
    real(8),dimension(:)             :: vloc !size[Nloc]
    real(8),dimension(:)             :: v    !size[N]
    integer                          :: i,iph,irank,Nloc,N
    integer                          :: v_start,v_end,vloc_start,vloc_end
    integer,dimension(:),allocatable :: Counts,Offset
    integer                          :: MpiSize,MpiIerr
    logical                          :: MpiMaster
    !
    if(  MpiComm == MPI_UNDEFINED .OR. MpiComm == Mpi_Comm_Null ) return
    ! stop "gather_vector_MPI error: MpiComm == MPI_UNDEFINED"
    !
    MpiSize   = get_size_MPI(MpiComm)
    MpiMaster = get_master_MPI(MpiComm)
    !
    Nloc = size(Vloc)
    N    = 0
    call AllReduce_MPI(MpiComm,Nloc,N)
    if(MpiMaster.AND.N /= size(V)) stop "allgather_vector_MPI error: size(V) != Mpi_Allreduce(Nloc)"
    !
    allocate(Counts(0:MpiSize-1)) ; Counts=0
    allocate(Offset(0:MpiSize-1)) ; Offset=0
    !
    !Get Counts;
    call MPI_AllGather(Nloc/Dimph,1,MPI_INTEGER,Counts,1,MPI_INTEGER,MpiComm,MpiIerr)
    !
    !Get Offset:
    Offset(0)=0
    do i=1,MpiSize-1
       Offset(i) = Offset(i-1) + Counts(i-1)
    enddo
    !
    V = 0d0
    do iph=1,Dimph
       v_start = 1 + (iph-1)*(N/Dimph)
       v_end = iph*(N/Dimph)
       vloc_start = 1 + (iph-1)*(Nloc/Dimph)
       vloc_end = iph*(Nloc/Dimph)
       call MPI_AllGatherv(Vloc(vloc_start:vloc_end),Nloc/DimPh,MPI_DOUBLE_PRECISION,&
            V(v_start:v_end),Counts,Offset,MPI_DOUBLE_PRECISION,MpiComm,MpiIerr)
    enddo
    !
    return
  end subroutine Allgather_vector_MPI
#endif





  !+------------------------------------------------------------------+
  !PURPOSE  : calculate the factorial of an integer N!=1.2.3...(N-1).N
  !+------------------------------------------------------------------+
  recursive function factorial(n) result(f)
    integer            :: f
    integer,intent(in) :: n
    if(n<=0)then
       f=1
    else
       f=n*factorial(n-1)
    end if
  end function factorial



  !+------------------------------------------------------------------+
  !PURPOSE  : calculate the binomial factor n1 over n2
  !+------------------------------------------------------------------+
  elemental function binomial(n1,n2) result(nchoos)
    integer,intent(in) :: n1,n2
    real(8)            :: xh
    integer            :: i
    integer nchoos
    xh = 1.d0
    if(n2<0) then
       nchoos = 0
       return
    endif
    if(n2==0) then
       nchoos = 1
       return
    endif
    do i = 1,n2
       xh = xh*dble(n1+1-i)/dble(i)
    enddo
    nchoos = int(xh + 0.5d0)
  end function binomial



  !+------------------------------------------------------------------+
  !PURPOSE : binary search of a value in an array
  !+------------------------------------------------------------------+
  recursive function binary_search(a,value) result(bsresult)
    integer,intent(in) :: a(:), value
    integer            :: bsresult, mid
    mid = size(a)/2 + 1
    if (size(a) == 0) then
       bsresult = 0        ! not found
       !stop "binary_search error: value not found"
    else if (a(mid) > value) then
       bsresult= binary_search(a(:mid-1), value)
    else if (a(mid) < value) then
       bsresult = binary_search(a(mid+1:), value)
       if (bsresult /= 0) then
          bsresult = mid + bsresult
       end if
    else
       bsresult = mid      ! SUCCESS!!
    end if
  end function binary_search







  !+------------------------------------------------------------------+
  !PURPOSE : sort array of integer using random algorithm
  !+------------------------------------------------------------------+
  subroutine sort_array(array)
    integer,dimension(:),intent(inout)      :: array
    integer,dimension(size(array))          :: order
    integer                                 :: i
    forall(i=1:size(array))order(i)=i
    call qsort_sort( array, order, 1, size(array) )
    array=order
  contains
    recursive subroutine qsort_sort( array, order, left, right )
      integer, dimension(:)                 :: array
      integer, dimension(:)                 :: order
      integer                               :: left
      integer                               :: right
      integer                               :: i
      integer                               :: last
      if ( left .ge. right ) return
      call qsort_swap( order, left, qsort_rand(left,right) )
      last = left
      do i = left+1, right
         if ( compare(array(order(i)), array(order(left)) ) .lt. 0 ) then
            last = last + 1
            call qsort_swap( order, last, i )
         endif
      enddo
      call qsort_swap( order, left, last )
      call qsort_sort( array, order, left, last-1 )
      call qsort_sort( array, order, last+1, right )
    end subroutine qsort_sort
    !---------------------------------------------!
    subroutine qsort_swap( order, first, second )
      integer, dimension(:)                 :: order
      integer                               :: first, second
      integer                               :: tmp
      tmp           = order(first)
      order(first)  = order(second)
      order(second) = tmp
    end subroutine qsort_swap
    !---------------------------------------------!
    function qsort_rand( lower, upper )
      implicit none
      integer                               :: lower, upper
      real(8)                               :: r
      integer                               :: qsort_rand
      call random_number(r)
      qsort_rand =  lower + nint(r * (upper-lower))
    end function qsort_rand
    function compare(f,g)
      integer                               :: f,g
      integer                               :: compare
      compare=1
      if(f<g)compare=-1
    end function compare
  end subroutine sort_array






  !+------------------------------------------------------------------+
  !PURPOSE  : print a state vector |{up}>|{dw}>
  !+------------------------------------------------------------------+
  subroutine print_state_vector_ivec(ivec,unit)
    integer,intent(in) :: ivec(:)
    integer,optional   :: unit
    integer            :: unit_
    integer            :: i,j,Ntot
    character(len=2)   :: fbt
    character(len=16)  :: fmt
    unit_=6;if(present(unit))unit_=unit
    Ntot = size(ivec)
    write(fbt,'(I2.2)')Ntot
    fmt="(B"//adjustl(trim(fbt))//"."//adjustl(trim(fbt))//")"
    i= bjoin(ivec,Ntot)
    write(unit_,"(I9,1x,A1)",advance="no")i,"|"
    write(unit_,"(10I1)",advance="no")(ivec(j),j=1,Ntot)
    write(unit_,"(A4)",advance="no")"> - "
    write(unit_,fmt,advance="yes")i
  end subroutine print_state_vector_ivec
  !
  subroutine  print_state_vector_ivec_ud(ivec,jvec,unit)
    integer,intent(in) :: ivec(:),jvec(size(ivec))
    integer,optional   :: unit
    integer            :: unit_
    integer            :: i,j,iup,idw,Ntot
    character(len=2)   :: fbt
    character(len=20)  :: fmt
    unit_=6;if(present(unit))unit_=unit
    Ntot = size(ivec)
    write(fbt,'(I2.2)')Ntot
    fmt="(B"//adjustl(trim(fbt))//"."//adjustl(trim(fbt))//",1x,B"//adjustl(trim(fbt))//"."//adjustl(trim(fbt))//")"
    iup = bjoin(ivec,Ntot)
    idw = bjoin(jvec,Ntot)
    i = bjoin([ivec,jvec],2*Ntot)
    write(unit_,"(I9,1x,I4,1x,A1)",advance="no")i,iup,"|"
    write(unit_,"(10I1)",advance="no")(ivec(j),j=1,Ntot)
    write(unit_,"(A1,I4,A2)",advance="no")">",idw," |"
    write(unit_,"(10I1)",advance="no")(jvec(j),j=1,Ntot)
    write(unit_,"(A4)",advance="no")"> - "
    write(unit_,fmt,advance="yes")ibits(i,0,Ntot),ibits(i,Ntot,2*Ntot)
  end subroutine print_state_vector_ivec_ud
  !
  subroutine print_state_vector_int(i,Ntot,unit)
    integer,intent(in) :: i
    integer,intent(in) :: Ntot
    integer,optional   :: unit
    integer            :: unit_
    integer            :: j
    integer            :: ivec(Ntot)
    character(len=2)   :: fbt
    character(len=16)  :: fmt
    unit_=6;if(present(unit))unit_=unit
    write(fbt,'(I2.2)')Ntot
    fmt="(B"//adjustl(trim(fbt))//"."//adjustl(trim(fbt))//")"
    ivec = bdecomp(i,Ntot)
    write(unit_,"(I9,1x,A1)",advance="no")i,"|"
    write(unit_,"(10I1)",advance="no")(ivec(j),j=1,Ntot)
    write(unit_,"(A4)",advance="no")"> - "
    write(unit_,fmt,advance="yes")i
  end subroutine print_state_vector_int

end MODULE ED_SETUP














! getCsector=0
! do isector=1,Nsectors
!    call get_Nup(isector,Nup)
!    call get_Ndw(isector,Ndw)
!    !
!    jup=nup-1; jdw=ndw; if(jup < 0)cycle
!    !
!    call get_Sector([jup,jdw],Ns,jsector)
!    getCsector(1,1,isector)=jsector
! enddo
! !
! !
! !
! do isector=1,Nsectors
!    call get_Nup(isector,Nup)
!    call get_Ndw(isector,Ndw)
!    !
!    jup=nup;jdw=ndw-1;if(jdw < 0)cycle
!    !
!    call get_Sector([jup,jdw],Ns,jsector)
!    getCsector(1,2,isector)=jsector
! enddo
! !
! !
! !
! getCDGsector=0
! do isector=1,Nsectors
!    call get_Nup(isector,Nup)
!    call get_Ndw(isector,Ndw)
!    !
!    jup=nup+1;jdw=ndw;if(jup > Ns)cycle
!    !
!    call get_Sector([jup,jdw],Ns,jsector)
!    getCDGsector(1,1,isector)=jsector
! enddo
! !
! !
! !
! do isector=1,Nsectors
!    call get_Nup(isector,Nup)
!    call get_Ndw(isector,Ndw)
!    !
!    jup=nup;jdw=ndw+1;if(jdw > Ns)cycle
!    !
!    call get_Sector([jup,jdw],Ns,jsector)
!    getCDGsector(1,2,isector)=jsector
! enddo
