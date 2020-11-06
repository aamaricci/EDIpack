module ED_EIGENSPACE
  USE ED_VARS_GLOBAL
  USE ED_SETUP
  implicit none
  private


  type full_espace
     real(8),dimension(:),pointer   :: e
     real(8),dimension(:,:),pointer :: M
  end type full_espace




  type sparse_estate
     integer                      :: sector        !index of the sector
     real(8)                      :: e             !energy of the eigen-state
     real(8),dimension(:),pointer :: cvec=>null()  !double precision eigen-vector
     logical                      :: itwin=.false. !twin sector label
     type(sparse_estate),pointer  :: twin=>null()  !link to twin box 
     type(sparse_estate),pointer  :: next=>null()  !link to next box (chain)
  end type sparse_estate

  type sparse_espace
     integer                      :: size
     real(8)                      :: emax,emin
     logical                      :: status=.false.
     type(sparse_estate),pointer  :: root=>null()       !head/root of the list\== list itself
  end type sparse_espace


  interface es_insert_state
     module procedure :: es_insert_state_c
  end interface es_insert_state

  interface es_add_state
     module procedure :: es_add_state_c
  end interface es_add_state

  interface es_return_cvector
     module procedure :: es_return_cvector_default
#ifdef _MPI
     module procedure :: es_return_cvector_mpi
#endif
  end interface es_return_cvector


  public :: sparse_estate
  public :: sparse_espace
  !
  public :: es_init_espace      !init the espace                 !checked
  public :: es_delete_espace    !del the espace                  !checked
  public :: es_free_espace      !free the espace                 !checked
  public :: es_print_espace     !print the espace                !checked
  !
  public :: es_insert_state     !insert a state                  !checked
  public :: es_add_state        !add a state w/ costraint        !checked
  public :: es_pop_state        !pop a state                     !checked
  !
  public :: es_return_sector       !get the sector of a state       !checked
  public :: es_return_energy       !get the energy of a state       !checked
  public :: es_return_cvector      !get the vector of a state       !checked
  public :: es_return_gs_degeneracy!get the number of degenerate GS !checked
  !
  public :: setup_eigenspace
  public :: delete_eigenspace
  !  
  !
  type(sparse_espace)                        :: state_list
  type(full_espace),dimension(:),allocatable :: espace
  public :: state_list
  public :: espace





contains        !some routine to perform simple operation on the lists





  !+------------------------------------------------------------------+
  !PURPOSE  : Setting up the Full ED eigen-Space
  !+------------------------------------------------------------------+
  subroutine setup_eigenspace
    integer :: isector,dim,jsector
    if(allocated(espace)) deallocate(espace)
    allocate(espace(1:Nsectors))
    do isector=1,Nsectors
       dim=GetDim(isector);if(dim==0)stop "setup_eigenspace: dim==0!"
       allocate(espace(isector)%e(dim))
       allocate(espace(isector)%M(dim,dim))
    enddo
  end subroutine setup_eigenspace





  !+------------------------------------------------------------------+
  !PURPOSE  : Deleting the Full ED eigen-Space (free the memory)
  !+------------------------------------------------------------------+
  subroutine delete_eigenspace
    integer :: isector
    if(allocated(espace))then
       do isector=1,size(espace)
          deallocate(espace(isector)%e)
          deallocate(espace(isector)%M)
       end do
       deallocate(espace)
    endif
  end subroutine delete_eigenspace





  !##################################################################
  !##################################################################
  !##################################################################
  !##################################################################




  !+------------------------------------------------------------------+
  !PURPOSE  : initialize the list of states
  !+------------------------------------------------------------------+
  function es_init_espace() result(space)
    type(sparse_espace) :: space
    allocate(space%root)
    space%status=.true.
    space%root%next => null()
    space%size=0
    space%emax=-huge(1d0)
    space%emin= huge(1d0)
  end function es_init_espace



  !+------------------------------------------------------------------+
  !PURPOSE  : destroy the list of states
  !+------------------------------------------------------------------+
  subroutine es_delete_espace(space)
    type(sparse_espace),intent(inout) :: space
    type(sparse_estate),pointer       :: p,c
    if(.not.space%status)return
    do
       p => space%root
       c => p%next
       if(.not.associated(c))exit  !empty list
       p%next => c%next !
       c%next=>null()
       if(associated(c%cvec))deallocate(c%cvec)
       if(associated(c%twin))c%twin=>null()
       deallocate(c)
    end do
    deallocate(space%root)
    space%status=.false.
    if(associated(p))nullify(p)
  end subroutine es_delete_espace




  !+------------------------------------------------------------------+
  !PURPOSE  : empty the list of states
  !+------------------------------------------------------------------+
  subroutine es_free_espace(space)
    type(sparse_espace),intent(inout) :: space
    type(sparse_estate),pointer       :: p,c
    do
       p => space%root
       c => p%next
       if(.not.associated(c))exit  !empty list
       p%next => c%next            !
       c%next=>null()
       if(associated(c%cvec))deallocate(c%cvec)
       if(associated(c%twin))c%twin=>null()
       deallocate(c)
    end do
    space%size=0
    space%emax=-huge(1.d0)
    space%emin=huge(1.d0)
    if(associated(p))nullify(p)
  end subroutine es_free_espace







  !+------------------------------------------------------------------+
  !PURPOSE  : insert a state into the list using ener,vector,sector
  !+------------------------------------------------------------------+
  subroutine es_add_state_c(espace,e,cvec,sector,twin,size,verbose)
    type(sparse_espace),intent(inout) :: espace
    real(8),intent(in)                :: e
    real(8),dimension(:),intent(in)   :: cvec
    integer,intent(in)                :: sector
    integer,intent(in),optional       :: size
    logical,intent(in),optional       :: verbose
    logical,intent(in),optional       :: twin
    logical                           :: twin_
    twin_=.false.;if(present(twin))twin_=twin
    if(present(size))then !if present size add respecting the size constraint.
       if(espace%size<size)then
          call es_insert_state_c(espace,e,cvec,sector,twin_)
       else
          if(e < es_return_energy(espace))then
             if(present(verbose).AND.(verbose.eqv..true.))print*,"found a new state:"
             call es_pop_state(espace)
             call es_insert_state_c(espace,e,cvec,sector,twin_)
          endif
       endif
    else                      !else add normally
       call es_insert_state_c(espace,e,cvec,sector,twin_)
    endif
  end subroutine es_add_state_c




  !+------------------------------------------------------------------+
  !PURPOSE  : insert a state into the list using ener,vector,sector
  !+------------------------------------------------------------------+
  subroutine es_insert_state_c(space,e,vec,sector,twin)
    type(sparse_espace),intent(inout) :: space
    real(8),intent(in)                :: e
    real(8),dimension(:),intent(in)   :: vec
    integer,intent(in)                :: sector
    logical                           :: twin
    type(sparse_estate),pointer       :: p,c
    p => space%root
    c => p%next
    do                            !traverse the list until e < value (ordered list)
       if(.not.associated(c))exit
       if(e <= c%e)exit
       p => c
       c => c%next
    end do
    !
    allocate(p%next)                !Create a new element in the list
    p%next%e = e
    if(e > space%emax)space%emax=e !update the max energy (corresponds to the top entry)
    if(e < space%emin)space%emin=e !update the min energy (corresponds to the first entry)
    allocate(p%next%cvec(size(vec)))
    p%next%cvec = vec
    p%next%itwin=.false.
    p%next%sector=sector
    space%size = space%size+1
    if(twin)then !Create a twin element in the list with same energy, no vector and twin flag T
       allocate(p%next%next)
       p%next%next%e = e    
       p%next%next%itwin=.true.
       p%next%next%sector=get_twin_sector(sector)
       p%next%next%twin => p%next      ! wiggled arrow of the twin_wout_vector points to its twin_w_vector
       p%next%twin      => p%next%next ! wiggled arrow of the twin_w_vector points to its twin_wout_vector
       space%size = space%size+1
    endif
    if(.not.associated(c))then !end of the list special case (current=>current%next)
       if(.not.twin)then
          p%next%next  => null()
       else
          p%next%next%next  => null()
       endif
    else
       if(.not.twin)then
          p%next%next  => c      !the %next of the new node come to current
       else
          p%next%next%next  => c
       endif
    end if
    if(associated(p))nullify(p)
    if(associated(c))nullify(c)
  end subroutine es_insert_state_c






  !+------------------------------------------------------------------+
  !PURPOSE  : remove last element from the list, if +n is given remove 
  ! the n-th element, if +e is given remove the state with state%e=e
  ! hint: CIRCLE = twin state (a state flagged with itwin=T bearing no vector)
  !       SQUARE = normal state (itwin=F bearing vector)
  !+------------------------------------------------------------------+
  subroutine es_pop_state(space,n)
    type(sparse_espace),intent(inout) :: space
    integer,optional,intent(in)       :: n
    integer                           :: i,pos
    type(sparse_estate),pointer       :: pp,p,c
    !
    pos= space%size ; if(present(n))pos=n
    !
    if(pos>space%size)stop "es_pop_state: pos > espace.size"
    if(space%size==0)stop "es_pop_state: empty list"
    pp => null()
    p  => null()
    c  => space%root
    do i=1,pos
       pp => p
       p  => c
       c  => c%next
       if(.not.associated(c))return !empty or end of the list
    end do
    !c is a circle, so the prev/next are necessarily squares: remove c and p (twins)
    if(c%itwin)then             
       pp%next => c%next
       !delete C
       if(associated(c%cvec))deallocate(c%cvec)
       if(associated(c%twin))c%twin=>null()
       deallocate(c)
       !delete P
       if(associated(p%cvec))deallocate(p%cvec)
       if(associated(p%twin))p%twin=>null()
       deallocate(p)
       p => pp
       space%size=space%size-2
    else
       !c is a square:
       !if c%next is associated:
       ! if it is a circle: delete c and c%next (twins)
       ! if it is a square: delete c
       !else c%next is not associated: delete c
       if(associated(c%next))then
          if(c%next%itwin)then
             p%next => c%next%next 
             !delete C
             if(associated(c%cvec))deallocate(c%cvec)
             if(associated(c%twin))c%twin=>null()
             deallocate(c)
             !delete C%NEXT
             if(associated(c%next%cvec))deallocate(c%next%cvec)
             if(associated(c%next%twin))c%next%twin=>null()
             deallocate(c%next)
             space%size=space%size-2
          else
             p%next => c%next      
             if(associated(c%cvec))deallocate(c%cvec)
             if(associated(c%twin))c%twin=>null()
             deallocate(c)
             space%size=space%size-1
          endif
       else
          p%next => c%next      
          if(associated(c%cvec))deallocate(c%cvec)
          if(associated(c%twin))c%twin=>null()
          deallocate(c)
          space%size=space%size-1
       endif
    endif
    if(space%size>0)then
       space%emax = es_return_energy(space,space%size)
       space%emin = es_return_energy(space,1)
    endif
    if(associated(pp))nullify(pp)
    if(associated(p))nullify(p)
    if(associated(c))nullify(c)
  end subroutine es_pop_state






  !+------------------------------------------------------------------+
  !PURPOSE  : 
  !+------------------------------------------------------------------+
  function es_return_gs_degeneracy(space,gsthreshold) result(numzero)
    type(sparse_espace),intent(in) :: space
    real(8),optional               :: gsthreshold
    real(8)                        :: gsthreshold_
    integer                        :: numzero,pos
    type(sparse_estate),pointer    :: c
    real(8)                        :: oldzero,enemin
    gsthreshold_=1.d-9;if(present(gsthreshold))gsthreshold_=gsthreshold
    if(.not.space%status) stop "es_return_gs_degeneracy: espace not allocated"
    oldzero=1000.d0
    numzero=0
    c => space%root
    pos=0
    do 
       c => c%next
       pos=pos+1
       if(.not.associated(c))exit !end of the list
       enemin = c%e
       if (enemin < oldzero-10.d0*gsthreshold_) then
          numzero=1
          oldzero=enemin
       elseif(abs(enemin-oldzero) <= gsthreshold_)then
          numzero=numzero+1
          oldzero=min(oldzero,enemin)
       endif
    end do
    if(associated(c))nullify(c)
  end function es_return_gs_degeneracy




  !+------------------------------------------------------------------+
  !PURPOSE  : 
  !+------------------------------------------------------------------+
  function es_return_sector(space,n) result(sector)
    type(sparse_espace),intent(in) :: space
    integer,optional,intent(in)    :: n
    integer                        :: sector
    type(sparse_estate),pointer    :: c
    integer                        :: i,pos
    if(.not.space%status) stop "es_return_sector: espace not allocated"
    pos= space%size ; if(present(n))pos=n
    if(pos>space%size)      stop "es_return_sector: n > espace.size"
    sector=0
    c => space%root
    do i=1,pos
       c => c%next
       if(.not.associated(c))exit
    end do
    if(space%size==0)return
    sector = c%sector
    if(associated(c))nullify(c)
  end function es_return_sector






  !+------------------------------------------------------------------+
  !PURPOSE  : 
  !+------------------------------------------------------------------+
  function es_return_energy(space,n) result(egs)
    type(sparse_espace),intent(in) :: space
    integer,optional,intent(in)    :: n
    real(8)                        :: egs
    type(sparse_estate),pointer    :: c
    integer                        :: i,pos
    if(.not.space%status) stop "es_return_energy: espace not allocated"
    pos= space%size ; if(present(n))pos=n
    if(pos>space%size)    stop "es_return_energy: n > espace.size"
    c => space%root
    egs=space%emax
    do i=1,pos
       c => c%next
       if(.not.associated(c))exit
    end do
    if(space%size==0)return
    if(.not.c%itwin)then
       egs = c%e
    else
       egs = c%twin%e
    endif
    if(associated(c))nullify(c)
  end function es_return_energy



  !+------------------------------------------------------------------+
  !PURPOSE  : 
  !+------------------------------------------------------------------+
  function es_return_cvector_default(space,n) result(vector)
    type(sparse_espace),intent(in)   :: space
    integer,optional,intent(in)      :: n
    real(8),dimension(:),pointer     :: vector
    type(sparse_estate),pointer      :: c
    integer                          :: i,pos
    integer                          :: dim
    integer,dimension(:),allocatable :: order
    !
    if(.not.space%status) stop "es_return_cvector ERRROR: espace not allocated"
    pos= space%size ; if(present(n))pos=n
    if(pos>space%size)      stop "es_return_cvector ERRROR: n > espace.size"
    if(space%size==0)stop "es_return_cvector ERRROR: espace emtpy"
    !
    c => space%root
    do i=1,pos
       c => c%next
       if(.not.associated(c))exit
    end do
    if(.not.c%itwin)then
       vector => c%cvec
    else
       dim = getdim(c%sector)
       allocate(Order(dim))
       call twin_sector_order(c%twin%sector,Order)
       allocate(vector(dim))
       do i=1,dim
          vector(i) = c%twin%cvec(Order(i))
       enddo
       deallocate(order)
    endif
    if(associated(c))nullify(c)
  end function es_return_cvector_default

#ifdef _MPI
  function es_return_cvector_mpi(MpiComm,space,n) result(vector)
    integer                          :: MpiComm
    type(sparse_espace),intent(in)   :: space
    integer,optional,intent(in)      :: n
    real(8),dimension(:),pointer     :: vtmp
    real(8),dimension(:),pointer     :: vector
    type(sparse_estate),pointer      :: c
    integer                          :: i,pos,Nloc,Ndim
    integer                          :: dim,ierr
    logical                          :: MpiMaster
    integer,dimension(:),allocatable :: order
    !
    if(MpiComm==MPI_COMM_NULL)return
    if(MpiComm==MPI_UNDEFINED)stop "es_return_cvector ERRROR: MpiComm = MPI_UNDEFINED"
    !
    if(.not.space%status) stop "es_return_cvector ERRROR: espace not allocated"
    pos= space%size ; if(present(n))pos=n
    if(pos>space%size)      stop "es_return_cvector ERRROR: n > espace.size"
    if(space%size==0)stop "es_return_cvector ERRROR: espace emtpy"
    !
    c => space%root
    do i=1,pos
       c => c%next
       if(.not.associated(c))exit
    end do
    !
    if(.not.c%itwin)then
       Nloc = size(c%cvec)
    else
       Nloc = size(c%twin%cvec)
    endif
    !
    !Ensure that the sum of the dimension of all vector chunks equals the sector dimension.
    Dim  = getdim(c%sector)
    Ndim = 0
    call Allreduce_MPI(MpiComm,Nloc,Ndim)
    if(Dim/=Ndim)stop "es_return_cvector ERROR: Dim != Ndim from v chunks"
    !
    MpiMaster = get_master_MPI(MpiComm)
    !
    if(.not.c%itwin)then
       if(MpiMaster)then
          allocate(Vector(Ndim))
       else
          allocate(Vector(1))
       endif
       Vector = 0d0
       call gather_vector_MPI(MpiComm,c%cvec,Vector)
    else
       !
       if(MpiMaster)then
          allocate(Vtmp(Ndim))
          allocate(Order(Dim))
          call twin_sector_order(c%twin%sector,Order)
       else
          allocate(Vtmp(1))
       endif
       Vtmp = 0d0
       call gather_vector_MPI(MpiComm,c%twin%cvec,Vtmp)
       if(MpiMaster)then
          allocate(Vector(Ndim))
          forall(i=1:Dim)Vector(i) = Vtmp(Order(i))
          deallocate(Order)
       else
          allocate(Vector(1))
          Vector = 0d0
       endif
       deallocate(Vtmp)
    endif
    if(associated(c))nullify(c)
  end function es_return_cvector_mpi
#endif





  !+------------------------------------------------------------------+
  !PURPOSE  : pretty print the list of states
  !+------------------------------------------------------------------+
  subroutine es_print_espace(space,unit,wvec)
    type(sparse_espace),intent(in) :: space
    type(sparse_estate),pointer    :: c
    integer                        :: counter,i
    integer,optional               :: unit
    integer                        :: unit_
    logical,optional               :: wvec
    logical                        :: wvec_
    wvec_=.false.;if(present(wvec))wvec_=wvec
    unit_=6;if(present(unit))unit_=unit
    write(*,"(A,I3)")"Print sparse espace unit ->",unit_
    c => space%root%next   !assume is associated,ie list exists
    counter = 0
    if(space%size>0)then
       do
          if(.not.associated(c))exit
          counter=counter+1
          write(unit_,"(A10,I5)")   "Index   : ",counter
          write(unit_,"(A10,I5)")   "Sector  : ",c%sector
          write(unit_,"(A10,3L3)")  "Twin    : ",c%itwin,associated(c%cvec)
          write(unit_,"(A10,I5)")   "Size    : ",getdim(c%sector)!size(c%vec)
          write(unit_,"(A10,f18.9)")"Energy  : ",c%e
          if(wvec_)then
             write(unit_,"(A10)")"Vec     : "
             do i=1,size(c%cvec)
                write(unit_,*)c%cvec(i)
             enddo
          endif
          c => c%next  !traverse list
          write(unit_,*)""
       end do
    else
       write(unit_,*)"Empty space"
       return
    endif
    if(associated(c))nullify(c)
  end subroutine es_print_espace



end module ED_EIGENSPACE
