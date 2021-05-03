MODULE ED_AUX_FUNX
  USE ED_INPUT_VARS
  USE ED_VARS_GLOBAL
  USE SF_TIMER
  USE SF_LINALG
  USE SF_MISC, only: assert_shape
  USE SF_IOTOOLS, only:free_unit,reg,txtfy
  USE SF_ARRAYS,  only: arange,linspace
  implicit none
  private


  interface lso2nnn_reshape
     module procedure d_nlso2nnn
     module procedure c_nlso2nnn
  end interface lso2nnn_reshape

  interface so2nn_reshape
     module procedure d_nso2nn
     module procedure c_nso2nn
  end interface so2nn_reshape

  interface nnn2lso_reshape
     module procedure d_nnn2nlso
     module procedure c_nnn2nlso
  end interface nnn2lso_reshape

  interface nn2so_reshape
     module procedure d_nn2nso
     module procedure c_nn2nso
  end interface nn2so_reshape

  interface ed_set_suffix
     module procedure :: ed_set_suffix_i
     module procedure :: ed_set_suffix_d
     module procedure :: ed_set_suffix_c
  end interface ed_set_suffix



  !FERMIONIC OPERATORS IN BITWISE OPERATIONS
  public :: c,cdg
  !AUX BIT OPERATIONS
  public :: bdecomp
  public :: breorder
  public :: bjoin
  !BINARY SEARCH
  public :: binary_search
  !MPI PROCEDURES
#ifdef _MPI
  public :: scatter_vector_MPI
  public :: scatter_basis_MPI
  public :: gather_vector_MPI
  public :: allgather_vector_MPI
#endif
  !AUX RESHAPE FUNCTIONS (internal use)
  public :: index_stride_so
  public :: lso2nnn_reshape
  public :: so2nn_reshape
  public :: nnn2lso_reshape
  public :: nn2so_reshape
  !
  !SEARCH CHEMICAL POTENTIAL, this should go into DMFT_TOOLS I GUESS
  public :: ed_search_variable
  public :: search_chemical_potential
  !ALLOCATE/DEALLOCATE GRIDS
  public :: allocate_grids
  public :: deallocate_grids

  public :: ed_set_suffix
  public :: ed_reset_suffix



contains


  subroutine ed_reset_suffix()
    ed_file_suffix=''
  end subroutine ed_reset_suffix

  subroutine ed_set_suffix_i(indx)
    integer :: indx
    ed_file_suffix=reg(ineq_site_suffix)//str(indx,site_indx_padding)
  end subroutine ed_set_suffix_i
  subroutine ed_set_suffix_d(indx)
    real(8) :: indx
    ed_file_suffix=reg(ineq_site_suffix)//str(indx)
  end subroutine ed_set_suffix_d
  subroutine ed_set_suffix_c(indx)
    character(len=*) :: indx
    ed_file_suffix=reg(ineq_site_suffix)//str(indx)
  end subroutine ed_set_suffix_c



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
  !               BITWISE OPERATIONS
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






  !##################################################################
  !                   RESHAPE ROUTINES
  !##################################################################
  !> Get stride position in the one-particle many-body space 
  function index_stride_so(ispin,iorb) result(indx)
    integer :: iorb
    integer :: ispin
    integer :: indx
    indx = iorb  + (ispin-1)*Norb
  end function index_stride_so


  !+-----------------------------------------------------------------------------+!
  !PURPOSE: 
  ! reshape a matrix from the [Nlso][Nlso] shape
  ! from/to the [Nlat][Nspin][Nspin][Norb][Norb] shape.
  ! _nlso2nnn : from [Nlso][Nlso] to [Nlat][Nspin][Nspin][Norb][Norb]  !
  ! _nso2nn   : from [Nso][Nso]   to [Nspin][Nspin][Norb][Norb]
  !+-----------------------------------------------------------------------------+!
  function d_nlso2nnn(Hlso,Nlat,Nspin,Norb) result(Hnnn)
    integer                                            :: Nlat,Nspin,Norb
    real(8),dimension(Nlat*Nspin*Norb,Nlat*Nspin*Norb) :: Hlso
    real(8),dimension(Nlat,Nspin,Nspin,Norb,Norb)      :: Hnnn
    integer                                            :: iorb,ispin,ilat,is
    integer                                            :: jorb,jspin,js
    Hnnn=zero
    do ilat=1,Nlat
       do ispin=1,Nspin
          do jspin=1,Nspin
             do iorb=1,Norb
                do jorb=1,Norb
                   is = iorb + (ispin-1)*Norb + (ilat-1)*Norb*Nspin !lattice-spin-orbit stride
                   js = jorb + (jspin-1)*Norb + (ilat-1)*Norb*Nspin !lattice-spin-orbit stride
                   Hnnn(ilat,ispin,jspin,iorb,jorb) = Hlso(is,js)
                enddo
             enddo
          enddo
       enddo
    enddo
  end function d_nlso2nnn
  function c_nlso2nnn(Hlso,Nlat,Nspin,Norb) result(Hnnn)
    integer                                               :: Nlat,Nspin,Norb
    complex(8),dimension(Nlat*Nspin*Norb,Nlat*Nspin*Norb) :: Hlso
    complex(8),dimension(Nlat,Nspin,Nspin,Norb,Norb)      :: Hnnn
    integer                                               :: iorb,ispin,ilat,is
    integer                                               :: jorb,jspin,js
    Hnnn=zero
    do ilat=1,Nlat
       do ispin=1,Nspin
          do jspin=1,Nspin
             do iorb=1,Norb
                do jorb=1,Norb
                   is = iorb + (ispin-1)*Norb + (ilat-1)*Norb*Nspin !lattice-spin-orbit stride
                   js = jorb + (jspin-1)*Norb + (ilat-1)*Norb*Nspin !lattice-spin-orbit stride
                   Hnnn(ilat,ispin,jspin,iorb,jorb) = Hlso(is,js)
                enddo
             enddo
          enddo
       enddo
    enddo
  end function c_nlso2nnn

  function d_nso2nn(Hso,Nspin,Norb) result(Hnn)
    integer                                  :: Nspin,Norb
    real(8),dimension(Nspin*Norb,Nspin*Norb) :: Hso
    real(8),dimension(Nspin,Nspin,Norb,Norb) :: Hnn
    integer                                  :: iorb,ispin,is
    integer                                  :: jorb,jspin,js
    Hnn=zero
    do ispin=1,Nspin
       do jspin=1,Nspin
          do iorb=1,Norb
             do jorb=1,Norb
                is = iorb + (ispin-1)*Norb  !spin-orbit stride
                js = jorb + (jspin-1)*Norb  !spin-orbit stride
                Hnn(ispin,jspin,iorb,jorb) = Hso(is,js)
             enddo
          enddo
       enddo
    enddo
  end function d_nso2nn
  function c_nso2nn(Hso,Nspin,Norb) result(Hnn)
    integer                                     :: Nspin,Norb
    complex(8),dimension(Nspin*Norb,Nspin*Norb) :: Hso
    complex(8),dimension(Nspin,Nspin,Norb,Norb) :: Hnn
    integer                                     :: iorb,ispin,is
    integer                                     :: jorb,jspin,js
    Hnn=zero
    do ispin=1,Nspin
       do jspin=1,Nspin
          do iorb=1,Norb
             do jorb=1,Norb
                is = iorb + (ispin-1)*Norb  !spin-orbit stride
                js = jorb + (jspin-1)*Norb  !spin-orbit stride
                Hnn(ispin,jspin,iorb,jorb) = Hso(is,js)
             enddo
          enddo
       enddo
    enddo
  end function c_nso2nn




  !+-----------------------------------------------------------------------------+!
  !PURPOSE: 
  ! reshape a matrix from the [Nlat][Nspin][Nspin][Norb][Norb] shape
  ! from/to the [Nlso][Nlso] shape.
  ! _nnn2nlso : from [Nlat][Nspin][Nspin][Norb][Norb] to [Nlso][Nlso]
  ! _nn2nso   : from [Nspin][Nspin][Norb][Norb]       to [Nso][Nso]
  !+-----------------------------------------------------------------------------+!
  function d_nnn2nlso(Hnnn,Nlat,Nspin,Norb) result(Hlso)
    integer                                            :: Nlat,Nspin,Norb
    real(8),dimension(Nlat,Nspin,Nspin,Norb,Norb)      :: Hnnn
    real(8),dimension(Nlat*Nspin*Norb,Nlat*Nspin*Norb) :: Hlso
    integer                                            :: iorb,ispin,ilat,is
    integer                                            :: jorb,jspin,js
    Hlso=zero
    do ilat=1,Nlat
       do ispin=1,Nspin
          do jspin=1,Nspin
             do iorb=1,Norb
                do jorb=1,Norb
                   is = iorb + (ispin-1)*Norb + (ilat-1)*Norb*Nspin !lattice-spin-orbit stride
                   js = jorb + (jspin-1)*Norb + (ilat-1)*Norb*Nspin !lattice-spin-orbit stride
                   Hlso(is,js) = Hnnn(ilat,ispin,jspin,iorb,jorb)
                enddo
             enddo
          enddo
       enddo
    enddo
  end function d_nnn2nlso

  function c_nnn2nlso(Hnnn,Nlat,Nspin,Norb) result(Hlso)
    integer                                            :: Nlat,Nspin,Norb
    complex(8),dimension(Nlat,Nspin,Nspin,Norb,Norb)      :: Hnnn
    complex(8),dimension(Nlat*Nspin*Norb,Nlat*Nspin*Norb) :: Hlso
    integer                                               :: iorb,ispin,ilat,is
    integer                                               :: jorb,jspin,js
    Hlso=zero
    do ilat=1,Nlat
       do ispin=1,Nspin
          do jspin=1,Nspin
             do iorb=1,Norb
                do jorb=1,Norb
                   is = iorb + (ispin-1)*Norb + (ilat-1)*Norb*Nspin !lattice-spin-orbit stride
                   js = jorb + (jspin-1)*Norb + (ilat-1)*Norb*Nspin !lattice-spin-orbit stride
                   Hlso(is,js) = Hnnn(ilat,ispin,jspin,iorb,jorb)
                enddo
             enddo
          enddo
       enddo
    enddo
  end function c_nnn2nlso

  function d_nn2nso(Hnn,Nspin,Norb) result(Hso)
    integer                                  :: Nspin,Norb
    real(8),dimension(Nspin,Nspin,Norb,Norb) :: Hnn
    real(8),dimension(Nspin*Norb,Nspin*Norb) :: Hso
    integer                                  :: iorb,ispin,is
    integer                                  :: jorb,jspin,js
    Hso=zero
    do ispin=1,Nspin
       do jspin=1,Nspin
          do iorb=1,Norb
             do jorb=1,Norb
                is = iorb + (ispin-1)*Norb  !spin-orbit stride
                js = jorb + (jspin-1)*Norb  !spin-orbit stride
                Hso(is,js) = Hnn(ispin,jspin,iorb,jorb)
             enddo
          enddo
       enddo
    enddo
  end function d_nn2nso

  function c_nn2nso(Hnn,Nspin,Norb) result(Hso)
    integer                                  :: Nspin,Norb
    complex(8),dimension(Nspin,Nspin,Norb,Norb) :: Hnn
    complex(8),dimension(Nspin*Norb,Nspin*Norb) :: Hso
    integer                                     :: iorb,ispin,is
    integer                                     :: jorb,jspin,js
    Hso=zero
    do ispin=1,Nspin
       do jspin=1,Nspin
          do iorb=1,Norb
             do jorb=1,Norb
                is = iorb + (ispin-1)*Norb  !spin-orbit stride
                js = jorb + (jspin-1)*Norb  !spin-orbit stride
                Hso(is,js) = Hnn(ispin,jspin,iorb,jorb)
             enddo
          enddo
       enddo
    enddo
  end function c_nn2nso














  !+------------------------------------------------------------------+
  !PURPOSE  : Allocate arrays and setup frequencies and times
  !+------------------------------------------------------------------+
  subroutine allocate_grids
    integer :: i
    if(.not.allocated(wm))allocate(wm(Lmats))
    if(.not.allocated(vm))allocate(vm(0:Lmats))          !bosonic frequencies
    if(.not.allocated(wr))allocate(wr(Lreal))
    if(.not.allocated(vr))allocate(vr(Lreal))
    if(.not.allocated(tau))allocate(tau(0:Ltau))
    wm     = pi/beta*real(2*arange(1,Lmats)-1,8)
    do i=0,Lmats
       vm(i) = pi/beta*2*i
    enddo
    wr     = linspace(wini,wfin,Lreal)
    vr     = linspace(0d0,wfin,Lreal)
    tau(0:)= linspace(0d0,beta,Ltau+1)
  end subroutine allocate_grids


  subroutine deallocate_grids
    if(allocated(wm))deallocate(wm)
    if(allocated(vm))deallocate(vm)
    if(allocated(tau))deallocate(tau)
    if(allocated(wr))deallocate(wr)
    if(allocated(vr))deallocate(vr)
  end subroutine deallocate_grids










  !##################################################################
  !##################################################################
  ! ROUTINES TO SEARCH CHEMICAL POTENTIAL UP TO SOME ACCURACY
  ! can be used to fix any other *var so that  *ntmp == nread
  !##################################################################
  !##################################################################

  !+------------------------------------------------------------------+
  !PURPOSE  : 
  !+------------------------------------------------------------------+
  subroutine ed_search_variable(var,ntmp,converged)
    real(8),intent(inout) :: var
    real(8),intent(in)    :: ntmp
    logical,intent(inout) :: converged
    logical               :: bool
    real(8),save          :: chich
    real(8),save          :: nold
    real(8),save          :: var_new
    real(8),save          :: var_old
    real(8)               :: var_sign
    real(8)               :: delta_n,delta_v,chi_shift
    !
    real(8)               :: ndiff
    integer,save          :: count=0,totcount=0,i
    integer               :: unit
    !
    !check actual value of the density *ntmp* with respect to goal value *nread*
    count=count+1
    totcount=totcount+1
    !  
    if(count==1)then
       chich = ndelta        !~0.2
       inquire(file="var_compressibility.restart",EXIST=bool)
       if(bool)then
          write(LOGfile,"(A)")"Reading compressibility from file"
          open(free_unit(unit),file="var_compressibility.restart")
          read(unit,*)chich
          close(unit)
       endif
       var_old = var
    endif
    !
    ndiff=ntmp-nread
    !
    open(free_unit(unit),file="var_compressibility.used")
    write(unit,*)chich
    close(unit)
    !
    ! if(abs(ndiff)>nerr)then
    !Get 'charge compressibility"
    delta_n = ntmp-nold
    delta_v = var-var_old
    if(count>1)chich = delta_v/(delta_n+1d-10) !1d-4*nerr)  !((ntmp-nold)/(var-var_old))**-1
    !
    !Add here controls on chich: not to be too small....
    if(chich>10d0)chich=2d0*chich/abs(chich) !do nothing?
    !
    chi_shift = ndiff*chich
    !
    !update chemical potential
    var_new = var - chi_shift
    !
    !
    !re-define variables:
    nold    = ntmp
    var_old = var
    var     = var_new
    !
    !Print information
    write(LOGfile,"(A11,F16.9,A,F15.9)")  "n      = ",ntmp,"| instead of",nread
    write(LOGfile,"(A11,ES16.9,A,ES16.9)")"n_diff = ",ndiff,"/",nerr
    write(LOGfile,"(A11,ES16.9,A,ES16.9)")"dv     = ",delta_v
    write(LOGfile,"(A11,ES16.9,A,ES16.9)")"dn     = ",delta_n
    write(LOGfile,"(A11,F16.9,A,F15.9)")  "dv/dn  = ",chich
    var_sign = (var-var_old)/abs(var-var_old)
    if(var_sign>0d0)then
       write(LOGfile,"(A11,ES16.9,A4)")"shift    = ",chi_shift," ==>"
    else
       write(LOGfile,"(A11,ES16.9,A4)")"shift    = ",chi_shift," <=="
    end if
    write(LOGfile,"(A11,F16.9)")"var     = ",var
    !
    ! else
    !    count=0
    ! endif


    !Save info about search variable iteration:
    open(free_unit(unit),file="search_variable_iteration_info"//reg(ed_file_suffix)//".ed",position="append")
    ! if(count==1)write(unit,*)"#var,ntmp,ndiff"
    write(unit,*)totcount,var,ntmp,ndiff
    close(unit)
    !
    !If density is not converged set convergence to .false.
    if(abs(ndiff)>nerr)converged=.false.
    !
    write(LOGfile,"(A18,I5)")"Search var count= ",count
    write(LOGfile,"(A19,L2)")"Converged       = ",converged
    print*,""
    !
    open(free_unit(unit),file="var_compressibility.restart")
    write(unit,*)chich
    close(unit)
    !
  end subroutine ed_search_variable




  subroutine search_chemical_potential(var,ntmp,converged)
    real(8),intent(inout) :: var
    real(8),intent(in)    :: ntmp
    logical,intent(inout) :: converged
    logical               :: bool
    real(8)               :: ndiff
    integer,save          :: count=0,totcount=0,i
    integer,save          :: nindex=0
    integer               :: nindex_old(3)
    real(8)               :: ndelta_old,nratio,kcompr=0d0
    integer,save          :: nth_magnitude=-2,nth_magnitude_old=-2
    real(8),save          :: nth=1.d-2,var_old,ntmp_old
    logical,save          :: ireduce=.true.
    integer               :: unit
    !
    if(count==0)then
       inquire(file="var.restart",EXIST=bool)
       if(bool)then
          open(free_unit(unit),file="var.restart")
          read(unit,*)var,ndelta
          ndelta=abs(ndelta)*ncoeff
          close(unit)
       endif
    endif
    !
    ndiff=ntmp-nread
    nratio = 0.5d0;!nratio = 1.d0/(6.d0/11.d0*pi)
    !
    !check actual value of the density *ntmp* with respect to goal value *nread*
    count=count+1
    totcount=totcount+1
    if(count>2)then
       do i=1,2
          nindex_old(i+1)=nindex_old(i)
       enddo
    endif
    nindex_old(1)=nindex
    !
    if(ndiff >= nth)then
       nindex=-1
    elseif(ndiff <= -nth)then
       nindex=1
    else
       nindex=0
    endif
    !
    !
    ndelta_old=ndelta
    bool=nindex/=0.AND.( (nindex+nindex_old(1)==0).OR.(nindex+sum(nindex_old(:))==0) )
    !if(nindex_old(1)+nindex==0.AND.nindex/=0)then !avoid loop forth and back
    if(bool)then
       ndelta=ndelta_old*nratio !decreasing the step
    else
       ndelta=ndelta_old
    endif
    !
    if(abs(ndelta_old)<1.d-9)then
       ndelta_old=0.d0
       nindex=0
    endif
    !
    !update chemical potential
    var=var+dble(nindex)*ndelta
    !
    ! if(count>0)kcompr = (ntmp - ntmp_old)/(var - var_old)       
    !
    !Print information
    write(LOGfile,"(A,f16.9,A,f15.9)")"n    = ",ntmp," /",nread
    if(nindex>0)then
       write(LOGfile,"(A,I0,A1,es16.9,A)")"shift= ",nindex,"*",ndelta," ==>"
    elseif(nindex<0)then
       write(LOGfile,"(A,I0,A1,es16.9,A)")"shift= ",nindex,"*",ndelta," <=="
    else
       write(LOGfile,"(A,I0,A1,es16.9,A)")"shift= ",nindex,"*",ndelta," ==="
    endif
    write(LOGfile,"(A,ES16.9,A,ES16.9)")"dn   = ",ndiff,"/",nth
    ! write(LOGfile,"(A,10F16.9)")"k    = ",kcompr,1d0/kcompr,ntmp,ntmp_old,ntmp-ntmp_old,var,var_old,var-var_old
    write(LOGfile,"(A,f15.9)")"var  = ",var
    !
    open(free_unit(unit),file="search_mu_iteration"//reg(ed_file_suffix)//".ed",position="append")
    write(unit,*)var,ntmp,ndiff
    close(unit)
    !
    !check convergence within actual threshold
    !if reduce is activetd
    !if density is in the actual threshold
    !if DMFT is converged
    !if threshold is larger than nerror (i.e. this is not last loop)
    bool=ireduce.AND.(abs(ndiff)<nth).AND.converged.AND.(nth>nerr)
    if(bool)then
       nth_magnitude_old=nth_magnitude        !save old threshold magnitude
       nth_magnitude=nth_magnitude_old-1      !decrease threshold magnitude || floor(log10(abs(ntmp-nread)))
       nth=max(nerr,10.d0**(nth_magnitude))   !set the new threshold 
       count=0                                !reset the counter
       converged=.false.                      !reset convergence
       ndelta=ndelta_old*nratio                  !reduce the delta step
       !
    endif
    !
    !if density is not converged set convergence to .false.
    if(abs(ntmp-nread)>nth)converged=.false.
    !
    !check convergence for this threshold
    !!---if smallest threshold-- NO MORE
    !if reduce is active (you reduced the treshold at least once)
    !if # iterations > max number
    !if not yet converged
    !set threshold back to the previous larger one.
    !bool=(nth==nerr).AND.ireduce.AND.(count>niter).AND.(.not.converged)
    bool=ireduce.AND.(count>niter).AND.(.not.converged)
    if(bool)then
       ireduce=.false.
       nth=10.d0**(nth_magnitude_old)
    endif
    !
    write(LOGfile,"(A,I5)")"count= ",count
    write(LOGfile,"(A,L2)")"Converged=",converged
    print*,""
    !
    open(free_unit(unit),file="var.restart")
    write(unit,*)var,ndelta
    close(unit)
    !
    ntmp_old = ntmp
    var_old  = var
    !
  end subroutine search_chemical_potential









END MODULE ED_AUX_FUNX
