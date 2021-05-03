MODULE ED_SECTOR
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



  public :: build_sector
  public :: delete_sector
  !
  public :: apply_op_C
  public :: apply_op_CDG
  public :: apply_op_Sz
  public :: apply_op_N
  public :: build_op_Ns

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
  public :: twin_sector_order
  public :: get_twin_sector
  public :: flip_state

contains


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



  subroutine build_op_Ns(i,Nup,Ndw,sectorI) 
    integer, intent(in)             :: i
    type(sector),intent(in)         :: sectorI
    integer,dimension(Ns)           :: Nup,Ndw  ![Ns]
    integer                         :: iph,i_el,ii,iorb
    integer,dimension(2*Ns_Ud)      :: Indices
    integer,dimension(Ns_Ud,Ns_Orb) :: Nups,Ndws  ![1,Ns]-[Norb,1+Nbath]
    integer,dimension(2*Ns)         :: Ib
    integer,dimension(2)            :: Iud
    !
    iph = (i-1)/(sectorI%DimEl) + 1
    i_el = mod(i-1,sectorI%DimEl) + 1
    !
    call state2indices(i_el,[sectorI%DimUps,sectorI%DimDws],Indices)
    do ii=1,Ns_Ud
       iud(1) = sectorI%H(ii)%map(Indices(ii))
       iud(2) = sectorI%H(ii+Ns_Ud)%map(Indices(ii+Ns_ud))
       Nups(ii,:) = Bdecomp(iud(1),Ns_Orb) ![Norb,1+Nbath]
       Ndws(ii,:) = Bdecomp(iud(2),Ns_Orb)
    enddo
    Nup = Breorder(Nups)
    Ndw = Breorder(Ndws)
    !
  end subroutine build_op_Ns










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

end MODULE ED_SECTOR














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
