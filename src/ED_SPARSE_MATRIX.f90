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


  !LOAD STANDARD MATRIX INTO SPARSE MATRICES
  interface sp_load_matrix
     module procedure :: sp_load_matrix_csr
#ifdef _MPI
     module procedure :: mpi_sp_load_matrix_csr
#endif
  end interface sp_load_matrix


  !DUMP SPARSE MATRIX INTO STANDARD MATRIX
  interface sp_dump_matrix
     module procedure :: sp_dump_matrix_csr
#ifdef _MPI
     module procedure :: mpi_sp_dump_matrix_csr
#endif
  end interface sp_dump_matrix


  !SPY PRINT SPARSE MATRIX
  interface sp_spy_matrix
     module procedure :: sp_spy_matrix_csr
#ifdef _MPI
     module procedure :: mpi_sp_spy_matrix_csr
#endif
  end interface sp_spy_matrix

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
  public :: sp_load_matrix      !create sparse from array !checked
  public :: sp_dump_matrix      !dump sparse into array   !checked
  public :: sp_spy_matrix       !
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
    if(.not.sparse%status)return !stop "Warning SPARSE/sp_delete_matrix: sparse not allocated already."
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
       ! row%vals = [row%vals,value]
       ! row%cols = [row%cols,column]
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
  !PURPOSE: load a regular matrix (2dim array) into a sparse matrix
  !+------------------------------------------------------------------+
  subroutine sp_load_matrix_csr(matrix,sparse)
    real(8),dimension(:,:),intent(in) :: matrix
    type(sparse_matrix_csr),intent(inout) :: sparse    
    integer                           :: i,j,Ndim1,Ndim2
    !
    Ndim1=size(matrix,1)
    Ndim2=size(matrix,2)   
    !
    if(sparse%status)call sp_delete_matrix_csr(sparse)
    call sp_init_matrix_csr(sparse,Ndim1,Ndim2)
    !
    do i=1,Ndim1
       do j=1,Ndim2
          if(matrix(i,j)/=0.d0)call sp_insert_element_csr(sparse,matrix(i,j),i,j)
       enddo
    enddo
  end subroutine sp_load_matrix_csr

#ifdef _MPI
  subroutine mpi_sp_load_matrix_csr(MpiComm,matrix,sparse)
    integer                               :: MpiComm
    real(8),dimension(:,:),intent(in)     :: matrix
    type(sparse_matrix_csr),intent(inout) :: sparse    
    integer                               :: i,j,Ndim1,Ndim2
    !
    if(MpiComm==Mpi_Comm_Null)return
    !
    call sp_test_matrix_mpi(MpiComm,sparse," mpi_sp_load_matrix_csr")
    !
    Ndim1=size(matrix,1)
    Ndim2=size(matrix,2)
    !
    if(sparse%status)call sp_delete_matrix_csr(sparse)
    call mpi_sp_init_matrix_csr(MpiComm,sparse,Ndim1,Ndim2)
    !
    do i=sparse%Istart,sparse%Iend
       do j=1,Ndim2
          if(matrix(i,j)/=0.d0)call mpi_sp_insert_element_csr(MpiComm,sparse,matrix(i,j),i,j)
       enddo
    enddo
  end subroutine mpi_sp_load_matrix_csr
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









  !+------------------------------------------------------------------+
  !PURPOSE: pretty print a sparse matrix on a given unit using format fmt
  !+------------------------------------------------------------------+  
  subroutine sp_spy_matrix_csr(sparse,header)
    type(sparse_matrix_csr)          :: sparse
    character ( len = * )           :: header
    integer                         :: N1,N2
    character ( len = 255 )         :: command_filename
    integer                         :: command_unit
    character ( len = 255 )         :: data_filename
    integer                         :: data_unit
    integer                         :: i, j
    character ( len = 6 )           :: n1_s,n2_s,n1_i,n2_i
    integer                         :: nz_num
    character ( len = 255 )         :: png_filename
    !
    !  Create data file.
    !
    !
    N1 = sparse%Nrow
    N2 = sparse%Ncol
    data_filename = trim ( header ) // '_data.dat'
    open (unit=free_unit(data_unit), file = data_filename, status = 'replace' )
    nz_num = 0
    do i=1,N1
       do j=1,sparse%row(i)%size
          write(data_unit,'(2x,i6,2x,i6)') sparse%row(i)%cols(j),i
          nz_num = nz_num + 1
       enddo
    enddo
    close(data_unit)
    !
    !  Create command file.
    !
    command_filename = "plot_"//str(header)//'_commands.gp'
    open(unit = free_unit(command_unit), file = command_filename, status = 'replace' )
    write(command_unit,'(a)') '#unset key'
    write(command_unit,'(a)') 'set terminal postscript eps enhanced color font "Times-Roman,16"'
    write(command_unit,'(a)') 'set output "|ps2pdf -sEPSCrop - '//str(header)//".pdf"//'"'
    write(command_unit,'(a)') 'set size ratio -1'
    write(command_unit,'(a)') 'set xlabel "<--- J --->"'
    write(command_unit,'(a)') 'set ylabel "<--- I --->"'
    write(command_unit,'(a,i6,a)')'set title "',nz_num,' nonzeros for '//str(header)//'"'
    write(command_unit,'(a)') 'set timestamp'
    write(command_unit,'(a)' )'plot [x=1:'//str(N1)//'] [y='//str(N2)//':1] "'//&
         str(data_filename)//'" w p pt 5 ps 0.4 lc rgb "red"'
    close ( unit = command_unit )
    return
  end subroutine sp_spy_matrix_csr



#ifdef _MPI
  subroutine mpi_sp_spy_matrix_csr(MpiComm,sparse,header)
    integer                         :: MpiComm
    type(sparse_matrix_csr)          :: sparse
    character ( len = * )           :: header
    integer                         :: N1,N2,N1_,N2_
    character ( len = 255 )         :: command_filename
    integer                         :: command_unit
    character ( len = 255 )         :: data_filename(1)
    integer                         :: data_unit
    integer                         :: i, j
    character ( len = 6 )           :: n1_s,n2_s,n1_i,n2_i
    integer                         :: nz_num,mpirank
    character ( len = 255 )         :: png_filename
    !
    if(MpiComm==Mpi_Comm_Null)return
    !
    call sp_test_matrix_mpi(MpiComm,sparse," mpi_sp_spy_matrix_csr")
    !
    MpiRank  = get_Rank_MPI(MpiComm)
    !
    !  Create data file.
    !
    N1_=sparse%Nrow
    N2_=sparse%Ncol
    N1=0
    N2=0
    call MPI_AllReduce(N1_,N1,1,MPI_Integer,MPI_SUM,MpiComm,MpiIerr)
    call MPI_AllReduce(N2_,N2,1,MPI_Integer,MPI_MAX,MpiComm,MpiIerr)
    !
    nz_num = 0
    !
    data_filename(1) = trim(header)//"_rank"//str(MpiRank,4)//'_matrix.dat'
    open(unit=free_unit(data_unit),file=data_filename(1), status = 'replace' )
    do i=1,sparse%Nrow
       do j=1,sparse%row(i)%Size
          write(data_unit,'(2x,i6,2x,i6)') sparse%row(i)%cols(j),i+sparse%Ishift
          nz_num = nz_num + 1
       enddo
    enddo
    write(data_unit,'(2x,i6,2x,i6)')
    close(data_unit)
    !
    !
    call MPI_Barrier(MpiComm,MpiIerr)
    !
    !  Create command file.
    !
    command_filename = "plot_"//trim(header)//"_rank"//str(MpiRank,4)//'_commands.gp'
    open(unit = free_unit(command_unit), file = command_filename, status = 'replace' )
    write(command_unit,'(a)') '#unset key'
    write(command_unit,'(a)') 'set terminal postscript eps enhanced color font "Times-Roman,16"'
    write(command_unit,'(a)') 'set output "|ps2pdf -sEPSCrop - '//str(header)//"_rank"//str(MpiRank,4)//".pdf"//'"'
    write(command_unit,'(a)') 'set size ratio -1'
    write(command_unit,'(a)') 'set xlabel "<--- J --->"'
    write(command_unit,'(a)') 'set ylabel "<--- I --->"'
    write(command_unit,'(a,i6,a)' ) &
         'set title "',nz_num,' nonzeros for '//str(header)//"_rank"//str(MpiRank,4)//'"'
    write(command_unit,'(a)') 'set timestamp'
    write(command_unit,'(a)' )'plot [x=1:'//str(N1)//'] [y='//str(N2)//':1] "'//&
         str(data_filename(1))//'" w p pt 5 ps 0.4 lc rgb "red"'
    close ( unit = command_unit )
    return
  end subroutine mpi_sp_spy_matrix_csr
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






