!>>> ACTHUNG: IN THIS VERSION THE LOCAL (in the memory) PART IS DISABLED <<<<
MODULE ED_SPARSE_MATRIX  !THIS VERSION CONTAINS ONLY DBLE ELEMENT: (SYMMETRIC MATRIX) 
  USE SF_IOTOOLS, only: str,free_unit
#ifdef _MPI
  USE SF_MPI
  USE MPI
#endif
  implicit none
  private



  type sparse_row_csr
     integer                                   :: size !actual 
     real(8),dimension(:),allocatable          :: vals
     integer,dimension(:),allocatable          :: cols
  end type sparse_row_csr

  type sparse_matrix_csr
     type(sparse_row_csr),dimension(:),pointer :: row
     integer                                   :: Nrow
     integer                                   :: Ncol
     logical                                   :: status=.false.
#ifdef _MPI
     integer                                   :: istart=0 !global start index for MPI storage
     integer                                   :: iend=0
     integer                                   :: ishift=0
     logical                                   :: mpi=.false.
#endif
  end type sparse_matrix_csr



  !INIT SPARSE MATRICES 
  interface sp_init_matrix
     module procedure :: sp_init_matrix_csr
#ifdef _MPI
     module procedure :: mpi_sp_init_matrix_csr
#endif
  end interface sp_init_matrix



  !DELETE SPARSE MATRIX 
  interface sp_delete_matrix
     module procedure :: sp_delete_matrix_csr
#ifdef _MPI
     module procedure :: mpi_sp_delete_matrix_csr
#endif
  end interface sp_delete_matrix


  !INSERT ELEMENTS
  interface sp_insert_element
     module procedure :: sp_insert_element_csr
#ifdef _MPI
     module procedure :: mpi_sp_insert_element_csr
#endif
  end interface sp_insert_element




  !DUMP SPARSE MATRIX INTO STANDARD MATRIX
  interface sp_dump_matrix
     module procedure :: sp_dump_matrix_csr
#ifdef _MPI
     module procedure :: mpi_sp_dump_matrix_csr
#endif
  end interface sp_dump_matrix



#ifdef _MPI  
  interface sp_set_mpi_matrix
     module procedure :: sp_set_mpi_matrix_csr
  end interface sp_set_mpi_matrix
#endif

  !Linked-List Sparse Matrix
  public :: sparse_matrix_csr

  public :: sp_init_matrix      !init the sparse matrix   !checked
  public :: sp_delete_matrix    !delete the sparse matrix !checked
  public :: sp_insert_element   !insert an element        !checked
  public :: sp_dump_matrix      !dump sparse into array   !checked
#ifdef _MPI
  public :: sp_set_mpi_matrix
#endif


  

  interface add_to
     module procedure :: add_to_I
     module procedure :: add_to_D
     module procedure :: add_to_Z
  end interface add_to



  integer :: MpiIerr


contains       


  !+------------------------------------------------------------------+
  !PURPOSE:  initialize the sparse matrix list
  !+------------------------------------------------------------------+
  subroutine sp_init_matrix_csr(sparse,N,N1)
    type(sparse_matrix_csr),intent(inout) :: sparse
    integer                               :: N
    integer,optional                      :: N1
    integer                               :: i
    !
    !put here a delete statement to avoid problems
    if(sparse%status)stop "sp_init_matrix: already allocated can not init"
    !
    sparse%Nrow=N
    sparse%Ncol=N 
    if(present(N1))sparse%Ncol=N1
    !
    allocate(sparse%row(N))
    do i=1,N
       sparse%row(i)%size=0
       allocate(sparse%row(i)%vals(0)) !empty array
       allocate(sparse%row(i)%cols(0)) !empty array
    end do
    !
    sparse%status=.true.
    !
  end subroutine sp_init_matrix_csr



#ifdef _MPI
  subroutine mpi_sp_init_matrix_csr(MpiComm,sparse,N,N1)
    integer                              :: MpiComm
    type(sparse_matrix_csr),intent(inout) :: sparse
    integer                              :: N
    integer,optional                     :: N1
    integer                              :: i,Ncol,Nloc
    !
    if(MpiComm==Mpi_Comm_Null)return
    !
    call sp_test_matrix_mpi(MpiComm,sparse,"mpi_sp_init_matrix_csr")
    !
    Ncol = N
    if(present(N1))Ncol=N1
    !
    Nloc = sparse%iend-sparse%istart+1
    !
    call sp_init_matrix_csr(sparse,Nloc,Ncol)
    !
  end subroutine mpi_sp_init_matrix_csr
#endif







  !+------------------------------------------------------------------+
  !PURPOSE: delete an entire sparse matrix
  !+------------------------------------------------------------------+
  subroutine sp_delete_matrix_csr(sparse)    
    type(sparse_matrix_csr),intent(inout) :: sparse
    integer                               :: i
    type(sparse_row_csr),pointer          :: row
    !
    if(.not.sparse%status)return !stop "Error SPARSE/sp_delete_matrix: sparse not allocated already."
    !
    do i=1,sparse%Nrow
       deallocate(sparse%row(i)%vals)
       deallocate(sparse%row(i)%cols)
       sparse%row(i)%Size  = 0
    enddo
    deallocate(sparse%row)
    !
    sparse%Nrow=0
    sparse%Ncol=0
    sparse%status=.false.
#ifdef _MPI
    sparse%istart = 0
    sparse%iend   = 0
    sparse%ishift = 0
    sparse%mpi    = .false.
#endif 
  end subroutine sp_delete_matrix_csr



#ifdef _MPI
  subroutine mpi_sp_delete_matrix_csr(MpiComm,sparse)
    integer                               :: MpiComm
    type(sparse_matrix_csr),intent(inout) :: sparse
    integer                               :: i
    !
    if(MpiComm==Mpi_Comm_Null)return
    !
    if(.not.sparse%status)return !stop "Error SPARSE/mpi_sp_delete_matrix: sparse is not allocated."
    !
    do i=1,sparse%Nrow
       deallocate(sparse%row(i)%vals)
       deallocate(sparse%row(i)%cols)
       sparse%row(i)%Size  = 0
    enddo
    deallocate(sparse%row)
    !
    sparse%Nrow=0
    sparse%Ncol=0
    sparse%status=.false.
    !
    sparse%istart=0
    sparse%iend=0
    sparse%ishift=0
    sparse%mpi=.false.
    !
  end subroutine mpi_sp_delete_matrix_csr
#endif    













  !+------------------------------------------------------------------+
  !PURPOSE: insert an element value at position (i,j) in the sparse matrix
  !+------------------------------------------------------------------+
  subroutine sp_insert_element_csr(sparse,value,i,j)
    type(sparse_matrix_csr),intent(inout) :: sparse
    real(8),intent(in)                    :: value
    integer,intent(in)                    :: i,j
    type(sparse_row_csr),pointer          :: row
    integer                               :: column,pos
    logical                               :: iadd
    !
    column = j
    !
    row => sparse%row(i)
    !
    iadd = .false.                          !check if column already exist
    if(any(row%cols == column))then         !
       pos = binary_search(row%cols,column) !find the position  column in %cols        
       iadd=.true.                          !set Iadd to true
    endif
    !
    if(iadd)then                            !this column exists so just sum it up       
       row%vals(pos)=row%vals(pos) + value  !add up value to the current one in %vals
    else                                    !this column is new. increase counter and store it 
       call add_to(row%vals,value)
       call add_to(row%cols,column)
       row%Size = row%Size + 1
    endif
    !
    if(row%Size > sparse%Ncol)stop "sp_insert_element_csr ERROR: row%Size > sparse%Ncol"
    !
  end subroutine sp_insert_element_csr

#ifdef _MPI
  subroutine mpi_sp_insert_element_csr(MpiComm,sparse,value,i,j)
    integer                               :: MpiComm
    type(sparse_matrix_csr),intent(inout) :: sparse
    real(8),intent(in)                    :: value
    integer,intent(in)                    :: i,j
    type(sparse_row_csr),pointer          :: row
    integer                               :: column,pos
    logical                               :: iadd
    !
    if(MpiComm==Mpi_Comm_Null)return
    !
    call sp_test_matrix_mpi(MpiComm,sparse," mpi_sp_insert_element_csr")
    !
    column = j
    !
    row => sparse%row(i-sparse%Ishift)
    !
    iadd = .false.                          !check if column already exist
    if(any(row%cols == column))then         !
       pos = binary_search(row%cols,column) !find the position  column in %cols        
       iadd=.true.                          !set Iadd to true
    endif
    !
    if(iadd)then                            !this column exists so just sum it up       
       row%vals(pos)=row%vals(pos) + value  !add up value to the current one in %vals
    else                                    !this column is new. increase counter and store it 
       ! row%vals = [row%vals,value]
       ! row%cols = [row%cols,column]
       call add_to(row%vals,value)
       call add_to(row%cols,column)
       row%Size = row%Size + 1
    endif
    !
    if(row%Size > sparse%Ncol)stop "mpi_sp_insert_element_csr ERROR: row%Size > sparse%Ncol"
    !
  end subroutine mpi_sp_insert_element_csr
#endif











  !+------------------------------------------------------------------+
  !PURPOSE: dump a sparse matrix into a regular 2dim array
  !+------------------------------------------------------------------+
  subroutine sp_dump_matrix_csr(sparse,matrix)
    type(sparse_matrix_csr),intent(in)   :: sparse
    real(8),dimension(:,:),intent(inout) :: matrix
    integer                              :: i,j,Ndim1,Ndim2
    !
    Ndim1=size(matrix,1)
    Ndim2=size(matrix,2)
    !
    if(sparse%Nrow/=Ndim1 .OR. sparse%Ncol/=Ndim2)stop "Warning SPARSE/dump_matrix: dimensions error"
    !
    do i=1,Ndim1
       do j=1,sparse%row(i)%Size
          matrix(i,sparse%row(i)%cols(j)) = matrix(i,sparse%row(i)%cols(j)) + sparse%row(i)%vals(j)
       enddo
    enddo
  end subroutine sp_dump_matrix_csr

#ifdef _MPI
  subroutine mpi_sp_dump_matrix_csr(MpiComm,sparse,matrix)
    integer                              :: MpiComm
    type(sparse_matrix_csr),intent(in)   :: sparse
    real(8),dimension(:,:),intent(inout) :: matrix
    real(8),dimension(:,:),allocatable   :: matrix_tmp
    integer                              :: i,impi,j,N1_,N2_,Ndim1,Ndim2,Nrow,Ncol
    !
    if(MpiComm==Mpi_Comm_Null)return
    !
    call sp_test_matrix_mpi(MpiComm,sparse," mpi_sp_dump_matrix_csr")
    !
    Ndim1=size(matrix,1)
    Ndim2=size(matrix,2)
    !
    N1_  = sparse%Nrow
    N2_  = sparse%Ncol
    Nrow = 0
    Ncol = 0
    call MPI_AllReduce(N1_,Nrow,1,MPI_Integer,MPI_SUM,MpiComm,MpiIerr)
    call MPI_AllReduce(N2_,Ncol,1,MPI_Integer,MPI_MAX,MpiComm,MpiIerr)
    !
    if(Nrow>Ndim1 .OR. Ncol>Ndim2)stop "Warning SPARSE/mpi_dump_matrix: dimensions error"
    !
    allocate(matrix_tmp(Ndim1,Ndim2));matrix_tmp=0d0
    do i=sparse%Istart,sparse%Iend
       impi = i - sparse%Ishift
       do j=1,sparse%row(impi)%Size
          matrix_tmp(i,sparse%row(impi)%cols(j))=matrix_tmp(i,sparse%row(impi)%cols(j))+sparse%row(impi)%vals(j)
       enddo
    enddo
    !
    ! Matrix=0d0
    call AllReduce_Mpi(MpiComm,Matrix_tmp,Matrix)
    !   
  end subroutine mpi_sp_dump_matrix_csr
#endif








#ifdef _MPI
  subroutine sp_set_mpi_matrix_csr(MpiComm,sparse,istart,iend,ishift)
    integer                              :: MpiComm
    type(sparse_matrix_csr),intent(inout) :: sparse
    integer                              :: istart,iend,ishift
    !
    if(MpiComm==Mpi_Comm_Null)return
    !
    sparse%istart = istart
    sparse%iend   = iend
    sparse%ishift = ishift
    sparse%mpi    = .true.
  end subroutine sp_set_mpi_matrix_csr

  subroutine sp_test_matrix_mpi(MpiComm,sparse,text)
    integer                              :: MpiComm
    type(sparse_matrix_csr),intent(in)    :: sparse
    character(len=*)                     :: text
    integer                              :: MpiRank
    !
    if(MpiComm==Mpi_Comm_Null)stop "sp_test_matrix_mpi ERROR: called in with MpiComm = Mpi_Comm_Null"
    !
    MpiRank = get_Rank_MPI(MpiComm)
    if(.not.sparse%mpi)then
       print*,"Rank, Error in "//trim(text)//": mpi no set"
       stop
    endif
  end subroutine sp_test_matrix_mpi
#endif







  !##################################################################
  !##################################################################
  !              AUXILIARY COMPUTATIONAL ROUTINES
  !##################################################################
  !##################################################################
  recursive function binary_search(Ain,value) result(bsresult)
    integer,intent(in)           :: Ain(:), value
    integer                      :: bsresult, mid
    integer,dimension(size(Ain)) :: A,Order
    !
    a = ain
    call sort_array(a,Order)
    !
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
    !
    bsresult = Order(bsresult)
    !
  end function binary_search





  subroutine add_to_I(vec,val)
    integer,dimension(:),allocatable,intent(inout) :: vec
    integer,intent(in)                             :: val  
    integer,dimension(:),allocatable               :: tmp
    integer                                        :: n
    !
    if (allocated(vec)) then
       n = size(vec)
       allocate(tmp(n+1))
       tmp(:n) = vec
       call move_alloc(tmp,vec)
       n = n + 1
    else
       n = 1
       allocate(vec(n))
    end if
    !
    !Put val as last entry:
    vec(n) = val
    !
    if(allocated(tmp))deallocate(tmp)
  end subroutine add_to_I

  subroutine add_to_D(vec,val)
    real(8),dimension(:),allocatable,intent(inout) :: vec
    real(8),intent(in)                             :: val  
    real(8),dimension(:),allocatable               :: tmp
    integer                                        :: n
    !
    if (allocated(vec)) then
       n = size(vec)
       allocate(tmp(n+1))
       tmp(:n) = vec
       call move_alloc(tmp,vec)
       n = n + 1
    else
       n = 1
       allocate(vec(n))
    end if
    !
    !Put val as last entry:
    vec(n) = val
    !
    if(allocated(tmp))deallocate(tmp)
  end subroutine add_to_D

  subroutine add_to_Z(vec,val)
    complex(8),dimension(:),allocatable,intent(inout) :: vec
    complex(8),intent(in)                             :: val  
    complex(8),dimension(:),allocatable               :: tmp
    integer                                           :: n
    !
    if (allocated(vec)) then
       n = size(vec)
       allocate(tmp(n+1))
       tmp(:n) = vec
       call move_alloc(tmp,vec)
       n = n + 1
    else
       n = 1
       allocate(vec(n))
    end if
    !
    !Put val as last entry:
    vec(n) = val
    !
    if(allocated(tmp))deallocate(tmp)
  end subroutine add_to_Z








  !+------------------------------------------------------------------+
  !PURPOSE  : Sort an array, gives the new ordering of the label.
  !+------------------------------------------------------------------+
  subroutine sort_array(array,order)
    implicit none
    integer,dimension(:)                    :: array
    integer,dimension(size(array))          :: order
    integer,dimension(size(array))          :: backup
    integer                                 :: i
    forall(i=1:size(array))order(i)=i
    call qsort_sort(array, order,1, size(array))
    do i=1,size(array)
       backup(i)=array(order(i))
    enddo
    array=backup
  contains
    recursive subroutine qsort_sort( array, order, left, right )
      integer, dimension(:) :: array
      integer, dimension(:) :: order
      integer               :: left
      integer               :: right
      integer               :: i
      integer               :: last
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
      integer, dimension(:) :: order
      integer               :: first, second
      integer               :: tmp
      tmp           = order(first)
      order(first)  = order(second)
      order(second) = tmp
    end subroutine qsort_swap
    !---------------------------------------------!
    integer function qsort_rand( lower, upper )
      integer               :: lower, upper
      real(8)               :: r
      call random_number(r)
      qsort_rand =  lower + nint(r * (upper-lower))
    end function qsort_rand
    !---------------------------------------------!
    function compare(f,g)
      implicit none
      integer               :: f,g
      integer               :: compare
      if(f<g) then
         compare=-1
      else
         compare=1
      endif
    end function compare
  end subroutine sort_array



end module ED_SPARSE_MATRIX






