MODULE ED_BATH
  USE SF_CONSTANTS, only: zero
  USE SF_IOTOOLS, only:free_unit,reg,file_length,txtfy
  USE SF_LINALG, only: eye,inv
  USE SF_ARRAYS, only:linspace
  USE SF_MISC, only: assert_shape
  USE ED_INPUT_VARS
  USE ED_VARS_GLOBAL
  USE ED_AUX_FUNX
  implicit none

  private



  !##################################################################
  !
  !     USER BATH ROUTINES:
  !
  !##################################################################

  !Interface for user bath I/O operations: get,set,copy
  interface get_bath_dimension
     module procedure ::  get_bath_dimension_direct
     module procedure ::  get_bath_dimension_symmetries
  end interface get_bath_dimension

  !explicit symmetries:
  interface break_symmetry_bath
     module procedure break_symmetry_bath_site
     module procedure break_symmetry_bath_lattice
  end interface break_symmetry_bath

  interface spin_symmetrize_bath
     module procedure spin_symmetrize_bath_site
     module procedure spin_symmetrize_bath_lattice
  end interface spin_symmetrize_bath

  interface orb_symmetrize_bath
     module procedure orb_symmetrize_bath_site
     module procedure orb_symmetrize_bath_lattice
  end interface orb_symmetrize_bath

  interface orb_equality_bath
     module procedure orb_equality_bath_site
     module procedure orb_equality_bath_lattice
  end interface orb_equality_bath

  interface ph_symmetrize_bath
     module procedure ph_symmetrize_bath_site
     module procedure ph_symmetrize_bath_lattice
  end interface ph_symmetrize_bath

  interface ph_trans_bath
     module procedure ph_trans_bath_site
     module procedure ph_trans_bath_lattice
  end interface ph_trans_bath


  !Aux:
  interface is_identity
     module procedure ::  is_identity_so
     module procedure ::  is_identity_nn
  end interface is_identity

  interface is_diagonal
     module procedure ::  is_diagonal_so
     module procedure ::  is_diagonal_nn
  end interface is_diagonal






  !##################################################################
  !
  !     USER BATH ROUTINES:
  !
  !##################################################################
  public  :: get_bath_dimension
  public  :: check_bath_dimension
  !
  public  :: get_bath_component_dimension
  public  :: get_bath_component
  public  :: set_bath_component
  public  :: copy_bath_component
  !
  public  :: break_symmetry_bath              
  public  :: spin_symmetrize_bath
  public  :: orb_symmetrize_bath
  public  :: orb_equality_bath
  public  :: ph_symmetrize_bath
  public  :: ph_trans_bath





  !##################################################################
  !
  !     DMFT BATH ROUTINES:
  !
  !##################################################################
  public  :: allocate_dmft_bath               !INTERNAL (for effective_bath)
  public  :: deallocate_dmft_bath             !INTERNAL (for effective_bath)
  public  :: init_dmft_bath                   !INTERNAL (for effective_bath)
  public  :: write_dmft_bath                  !INTERNAL (for effective_bath)
  public  :: save_dmft_bath                   !INTERNAL (for effective_bath)
  public  :: set_dmft_bath                    !INTERNAL (for effective_bath)
  public  :: get_dmft_bath                    !INTERNAL (for effective_bath)
  public  :: bath_from_sym                    !INTERNAL (for effective_bath)
  public  :: mask_hloc



  integer :: ibath,ilat,iorb



contains



  !+-------------------------------------------------------------------+
  !PURPOSE  : Inquire the correct bath size to allocate the 
  ! the bath array in the calling program.
  !
  ! Get size of each dimension of the component array. 
  ! The Result is an rank 1 integer array Ndim with dimension:
  ! 3 for get_component_size_bath
  ! 2 for get_spin_component_size_bath & get_orb_component_size_bath
  ! 1 for get_spin_orb_component_size_bath
  !+-------------------------------------------------------------------+
  function get_bath_dimension_direct(Hloc_nn,ispin_) result(bath_size)
    complex(8),optional,intent(in) :: Hloc_nn(:,:,:,:)
    integer,optional               :: ispin_
    integer                        :: bath_size,ndx,ispin,iorb,jspin,jorb,io,jo,counter
    real(8),allocatable            :: Hloc(:,:,:,:)
    select case(bath_type)
    case default
       !( e [Nspin][Norb][Nbath] + v [Nspin][Norb][Nbath] )
       bath_size = Norb*Nbath + Norb*Nbath
       if(.not.present(ispin_))bath_size=Nspin*bath_size
       !
       !
    case('hybrid')
       !(e [Nspin][1][Nbath] + v [Nspin][Norb][Nbath] )
       bath_size = Nbath + Norb*Nbath
       if(.not.present(ispin_))bath_size=Nspin*bath_size
       !
       !
    case('replica')
       !
       if(present(Hloc_nn))then
          allocate(Hloc(Nspin,Nspin,Norb,Norb));Hloc=dreal(Hloc_nn)
       elseif(allocated(impHloc))then
          allocate(Hloc(Nspin,Nspin,Norb,Norb));Hloc=impHloc
       else
          stop "ERROR: get_bath_dimension: bath_type=replica neither Hloc_nn present nor impHloc allocated"
       endif
       !
       !Real part of nonzero elements
       ndx=0
       do ispin=1,Nspin
          do jspin=1,Nspin
             do iorb=1,Norb
                do jorb=1,Norb
                   io=index_stride_so(ispin,iorb)
                   jo=index_stride_so(jspin,jorb)
                   if(io > jo)cycle
                   if((Hloc(ispin,jspin,iorb,jorb)/=0d0))ndx=ndx+1
                enddo
             enddo
          enddo
       enddo
       ndx = ndx * Nbath !number of non vanishing elements for each replica
       ndx = ndx + Nbath !diagonal hybridizations: Vs (different per spin)
       ndx = ndx + 1     !we also print Nbasis
       bath_size = ndx
    end select
  end function get_bath_dimension_direct

  function get_bath_dimension_symmetries(Hloc_nn) result(bath_size)
    complex(8),dimension(:,:,:,:,:),intent(in) :: Hloc_nn
    integer                                    :: bath_size,ndx,isym,Nsym
    !
    !number of symmetries
    Nsym=size(Hloc_nn,5)
    !
    ndx=Nsym
    !
    !number of replicas
    ndx = ndx * Nbath
    !diagonal hybridizations: Vs
    ndx = ndx + Nbath
    !
    !include Nbasis
    ndx=ndx+1
    !
    bath_size = ndx
    !
  end function get_bath_dimension_symmetries


  !+-------------------------------------------------------------------+
  !PURPOSE  : Check if the dimension of the bath array are consistent
  !+-------------------------------------------------------------------+
  function check_bath_dimension(bath_,Hloc_nn) result(bool)
    real(8),dimension(:)        :: bath_
    integer                     :: Ntrue,i
    logical                     :: bool
    real(8),optional,intent(in) :: Hloc_nn(:,:,:,:)
    real(8),allocatable         :: Hbasis_rebuild(:,:,:,:,:)![Nspin][:][Norb][:][Nsym]
    select case (bath_type)
    case default
       if (present(Hloc_nn))then
          Ntrue = get_bath_dimension(one*Hloc_nn)
       else
          Ntrue = get_bath_dimension()
       endif
    case ('replica')
       if(.not.allocated(H_basis))STOP "check_bath_dimension: Hbasis not allocated"
       if(.not.allocated(Hbasis_rebuild))allocate(Hbasis_rebuild(Nspin,Nspin,Norb,Norb,size(H_basis)))
       do i=1,size(H_basis)
          Hbasis_rebuild(:,:,:,:,i)=H_basis(i)%O
       enddo
       Ntrue   = get_bath_dimension(one*Hbasis_rebuild)
    end select
    bool  = ( size(bath_) == Ntrue )
  end function check_bath_dimension






  !##################################################################
  !
  !     USER BATH ROUTINES:
  !
  !##################################################################
  include 'user_aux.f90'




  !##################################################################
  !
  !     DMFT BATH ROUTINES:
  !
  !##################################################################
  include 'dmft_aux.f90'






  !##################################################################
  !
  !     AUX FUNCTIONS:
  !
  !##################################################################
  !+-------------------------------------------------------------------+
  !PURPOSE  : Check if a matrix is the identity
  !+-------------------------------------------------------------------+
  function is_identity_nn(mnnn) result(flag)
    real(8),dimension(nspin,nspin,norb,norb) :: mnnn
    real(8),dimension(nspin*norb,nspin*norb) :: mtmp
    integer                                  :: i,j
    logical                                  :: flag
    !
    flag=.true.
    !
    mtmp=nn2so_reshape(mnnn,nspin,norb)
    !
    do i=1,Nspin*Norb-1
       if((mtmp(i,i).ne.mtmp(i+1,i+1)).or.(mtmp(i,i).lt.1.d-6))flag=.false.
    enddo
    !
    do i=1,Nspin*Norb
       do j=1,Nspin*Norb
          if((i.ne.j).and.(mtmp(i,j).gt.1.d-6))flag=.false.
       enddo
    enddo
    !
  end function is_identity_nn

  function is_identity_so(mlso) result(flag)
    real(8),dimension(nspin*norb,nspin*norb) :: mlso
    real(8),dimension(nspin*norb,nspin*norb) :: mtmp
    integer                                  :: i,j
    logical                                  :: flag
    !
    flag=.true.
    !
    mtmp=mlso
    !
    do i=1,Nspin*Norb-1
       if((mtmp(i,i).ne.mtmp(i+1,i+1)).or.(mtmp(i,i).lt.1.d-6))flag=.false.
    enddo
    !
    do i=1,Nspin*Norb
       do j=1,Nspin*Norb
          if((i.ne.j).and.(mtmp(i,j).gt.1.d-6))flag=.false.
       enddo
    enddo
    !
  end function is_identity_so



  !+-------------------------------------------------------------------+
  !PURPOSE  : Check if a matrix is diagonal
  !+-------------------------------------------------------------------+
  function is_diagonal_nn(mnnn) result(flag)
    real(8),dimension(nspin,nspin,norb,norb) :: mnnn
    real(8),dimension(nspin*norb,nspin*norb) :: mtmp
    integer                                  :: i,j
    logical                                  :: flag
    !
    flag=.true.
    !
    mtmp=abs((nn2so_reshape(mnnn,nspin,norb)))
    !
    do i=1,Nspin*Norb
       do j=1,Nspin*Norb
          if((i.ne.j).and.(mtmp(i,j).gt.1.d-6))flag=.false.
       enddo
    enddo
    !
  end function is_diagonal_nn

  function is_diagonal_so(mlso) result(flag)
    real(8),dimension(Nspin*Norb,Nspin*Norb) :: mlso
    real(8),dimension(Nspin*Norb,Nspin*Norb) :: mtmp
    integer                                  :: i,j
    logical                                  :: flag
    !
    flag=.true.
    !
    mtmp=abs((mlso))
    !
    do i=1,Nspin*Norb
       do j=1,Nspin*Norb
          if((i.ne.j).and.(mtmp(i,j).gt.1.d-6))flag=.false.
       enddo
    enddo
    !
  end function is_diagonal_so



  !+-------------------------------------------------------------------+
  !PURPOSE  : Create bath mask
  !+-------------------------------------------------------------------+
  function mask_hloc(hloc,wdiag,uplo) result(Hmask)
    real(8),dimension(Nspin,Nspin,Norb,Norb) :: Hloc
    logical,optional                         :: wdiag,uplo
    logical                                  :: wdiag_,uplo_
    logical,dimension(Nspin,Nspin,Norb,Norb) :: Hmask
    integer                                  :: iorb,jorb,ispin,io,jo
    !
    wdiag_=.false.;if(present(wdiag))wdiag_=wdiag
    uplo_ =.false.;if(present(uplo))  uplo_=uplo
    !
    Hmask=.false.
    where(abs(Hloc)>1d-6)Hmask=.true.
    !
    !
    if(wdiag_)then
       do ispin=1,Nspin
          do iorb=1,Norb
             Hmask(ispin,ispin,iorb,iorb)=.true.
          enddo
       enddo
    endif
    !
    if(uplo_)then
       do ispin=1,Nspin
          do iorb=1,Norb
             do jorb=1,Norb
                io = index_stride_so(ispin,iorb)
                jo = index_stride_so(ispin,jorb)
                if(io>jo)Hmask(ispin,ispin,iorb,jorb)=.false.
             enddo
          enddo
       enddo
    endif
    !
  end function mask_hloc



END MODULE ED_BATH
