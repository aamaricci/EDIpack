MODULE ED_SETUP
  USE ED_INPUT_VARS
  USE ED_VARS_GLOBAL
  USE ED_AUX_FUNX
  USE ED_SECTOR
  USE SF_TIMER
  USE SF_IOTOOLS, only:free_unit,reg,file_length
#ifdef _MPI
  USE MPI
  USE SF_MPI
#endif
  implicit none
  private



  public :: init_ed_structure
  public :: setup_global

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
    endif
    !
    if(Nspin>1.AND.ed_twin.eqv..true.)then
       write(LOGfile,"(A)")"WARNING: using twin_sector with Nspin>1"
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
    if(ed_finite_temp)then
       if(lanc_nstates_total==1)stop "ED ERROR: ed_finite_temp=T *but* lanc_nstates_total==1 => T=0. Increase lanc_nstates_total"
    else
       if(lanc_nstates_total>1)print*, "ED WARNING: ed_finite_temp=F, T=0 *AND* lanc_nstates_total>1. re-Set lanc_nstates_total=1"
       lanc_nstates_total=1
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
    allocate(ed_dens(Norb),ed_docc(Norb),ed_dens_up(Norb),ed_dens_dw(Norb),ed_mag(Norb))
    ed_dens=0d0
    ed_docc=0d0
    ed_dens_up=0d0
    ed_dens_dw=0d0
    ed_mag =0d0
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
  !SECTOR PROCEDURES - Sectors,Nup,Ndw,DimUp,DimDw,...
  !##################################################################
  !##################################################################
  elemental function get_sector_dimension(n,np) result(dim)
    integer,intent(in) :: n,np
    integer            :: dim
    dim = binomial(n,np)
  end function get_sector_dimension






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



end MODULE ED_SETUP
