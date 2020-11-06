MODULE ED_HLOC_DECOMPOSITION
  USE ED_INPUT_VARS
  USE ED_VARS_GLOBAL
  USE ED_AUX_FUNX
  USE SF_IOTOOLS, only:free_unit,reg
  implicit none
  private


  interface print_Hloc
     module procedure print_Hloc_so
     module procedure print_Hloc_nn
  end interface print_Hloc

  interface set_Hloc
     module procedure init_Hloc_direct_so
     module procedure init_Hloc_direct_nn
     module procedure init_Hloc_symmetries
  end interface set_Hloc

  public:: set_Hloc
  public:: allocate_h_basis
  public:: deallocate_h_basis

contains

  !-------------------------------------------------------------------!
  ! PURPOSE: INITIALIZE INTERNAL HLOC STRUCTURES
  !-------------------------------------------------------------------!

  !allocate GLOBAL basis for H (used for impHloc and bath) and vectors coefficient
  subroutine allocate_h_basis(N)
    integer              :: N,isym
    !
    if(allocated(H_basis))deallocate(H_basis)
    if(allocated(lambda_impHloc))deallocate(lambda_impHloc)
    !
    allocate(H_basis(N))
    allocate(lambda_impHloc(N))
    do isym=1,N
       allocate(H_basis(isym)%O(Nspin,Nspin,Norb,Norb))
       H_basis(isym)%O=0d0
       lambda_impHloc(isym)=0d0
    enddo
  end subroutine allocate_h_basis


  !deallocate GLOBAL basis for H (used for impHloc and bath) and vectors coefficient

  subroutine deallocate_h_basis()
    integer              :: isym
    !
    do isym=1,size(H_basis)
       deallocate(H_basis(isym)%O)
    enddo
    deallocate(H_basis)
  end subroutine deallocate_h_basis

  !reconstruct [Nspin][][Norb][] hamiltonian from basis expansion given [lambda]
  function H_from_sym(lambdavec) result (H)
    real(8),dimension(:)                     :: lambdavec
    integer                                  :: isym
    real(8),dimension(Nspin,Nspin,Norb,Norb) :: H
    !
    if(size(lambdavec).ne.size(H_basis)) STOP "H_from_sym: Wrong coefficient vector size"
    H=zero
    do isym=1,size(lambdavec)
       H=H+lambdavec(isym)*H_basis(isym)%O
    enddo
  end function H_from_sym

  !initialize impHloc and the set [H_basis,lambda_impHloc]
  subroutine init_hloc_direct_so(Hloc)
    integer                                     :: ispin,jspin,iorb,jorb,counter,io,jo,Nsym
    complex(8),dimension(Nspin*Norb,Nspin*Norb) :: Hloc
    logical(8),dimension(Nspin,Nspin,Norb,Norb) :: Hmask
    !
    allocate(impHloc(Nspin,Nspin,Norb,Norb))
    impHloc=0d0
    !
    impHloc=so2nn_reshape(Hloc,Nspin,Norb)
    !
    Hmask=.false.
    !SPIN DIAGONAL
    do ispin=1,Nspin
       do jspin=1,Nspin
          do iorb=1,Norb
             do jorb=1,Norb
                io=index_stride_so(ispin,iorb)
                jo=index_stride_so(jspin,jorb)
                if((impHloc(ispin,jspin,iorb,jorb).ne.0d0).and.(io.le.jo))then
                   if(impHloc(ispin,jspin,iorb,jorb).ne.0.d0)counter=counter+1
                endif
             enddo
          enddo
       enddo
    enddo
    !
    call allocate_h_basis(counter)
    !
    counter=0
    !
    do ispin=1,Nspin
       do jspin=1,Nspin
          do iorb=1,Norb
             do jorb=1,Norb
                io=index_stride_so(ispin,iorb)
                jo=index_stride_so(jspin,jorb)
                if((impHloc(ispin,jspin,iorb,jorb).ne.0d0).and.(io.le.jo))then
                   if(impHloc(ispin,jspin,iorb,jorb).ne.0.d0)then
                      counter=counter+1
                      H_basis(counter)%O(ispin,jspin,iorb,jorb)=1d0
                      H_basis(counter)%O(ispin,jspin,jorb,iorb)=1d0
                      lambda_impHloc(counter)=impHloc(ispin,ispin,iorb,jorb)
                   endif
                   !
                endif
             enddo
          enddo
       enddo
    enddo
  end subroutine init_hloc_direct_so


  subroutine init_hloc_direct_nn(Hloc)
    integer                                     :: ispin,jspin,iorb,jorb,counter,io,jo,Nsym
    complex(8),dimension(Nspin,Nspin,Norb,Norb) :: Hloc
    logical(8),dimension(Nspin,Nspin,Norb,Norb) :: Hmask
    !
    if(.not.allocated(impHloc))allocate(impHloc(Nspin,Nspin,Norb,Norb))
    impHloc=zero
    !
    impHloc=Hloc
    !
    counter=0
    !
    Hmask=.false.
    !SPIN DIAGONAL
    do ispin=1,Nspin
       do jspin=1,Nspin
          do iorb=1,Norb
             do jorb=1,Norb
                io=index_stride_so(ispin,iorb)
                jo=index_stride_so(jspin,jorb)
                if((impHloc(ispin,jspin,iorb,jorb).ne.0d0).and.(io.le.jo))then
                   if(impHloc(ispin,jspin,iorb,jorb).ne.0.d0)counter=counter+1
                endif
             enddo
          enddo
       enddo
    enddo
    !
    call allocate_h_basis(counter)
    !
    counter=0
    !
    do ispin=1,Nspin
       do jspin=1,Nspin
          do iorb=1,Norb
             do jorb=1,Norb
                io=index_stride_so(ispin,iorb)
                jo=index_stride_so(jspin,jorb)
                if((impHloc(ispin,jspin,iorb,jorb).ne.0d0).and.(io.le.jo))then
                   if(impHloc(ispin,jspin,iorb,jorb).ne.0.d0)then
                      counter=counter+1
                      H_basis(counter)%O(ispin,jspin,iorb,jorb)=1d0
                      H_basis(counter)%O(ispin,jspin,jorb,iorb)=1d0
                      lambda_impHloc(counter)=impHloc(ispin,ispin,iorb,jorb)
                   endif
                   !
                endif
             enddo
          enddo
       enddo
    enddo
  end subroutine init_hloc_direct_nn


  subroutine init_hloc_symmetries(Hvec,lambdavec)
    integer                         :: isym,N
    complex(8),dimension(:,:,:,:,:) :: Hvec
    real(8),dimension(:)            :: lambdavec
    !
    allocate(impHloc(Nspin,Nspin,Norb,Norb))
    impHloc=0d0
    !
    N=size(lambdavec)
    !
    call allocate_h_basis(N)
    !
    do isym=1,N
       lambda_impHloc(isym)= lambdavec(isym)
       H_basis(isym)%O     = Hvec(:,:,:,:,isym)
    enddo
    !
    impHloc=H_from_sym(lambda_impHloc)
    !
  end subroutine init_hloc_symmetries








  !+------------------------------------------------------------------+
  !PURPOSE  : Print Hloc
  !+------------------------------------------------------------------+
  subroutine print_Hloc_nn(hloc,file)![Nspin][Nspin][Norb][Norb]
    real(8),dimension(Nspin,Nspin,Norb,Norb) :: hloc
    character(len=*),optional     :: file
    integer                       :: iorb,jorb,ispin,jspin
    integer                       :: unit
    unit=LOGfile
    if(present(file))then
       open(free_unit(unit),file=reg(file))
       write(LOGfile,"(A)")"print_Hloc on file :"//reg(file)
    endif
    do ispin=1,Nspin
       do iorb=1,Norb
          write(unit,"(20(F7.3,2x))")&
               ((Hloc(ispin,jspin,iorb,jorb),jorb =1,Norb),jspin=1,Nspin)
       enddo
    enddo
    write(unit,*)""
    if(present(file))close(unit)
  end subroutine print_Hloc_nn

  subroutine print_Hloc_so(hloc,file) ![Nlso][Nlso]
    real(8),dimension(Nspin*Norb,Nspin*Norb) :: hloc
    character(len=*),optional                :: file
    integer                                  :: iorb,jorb,unit,Nso
    unit=LOGfile
    if(present(file))then
       open(free_unit(unit),file=reg(file))
       write(LOGfile,"(A)")"print_Hloc on file :"//reg(file)
    endif
    !
    Nso = Nspin*Norb
    do iorb=1,Nso
       write(unit,"(20(F7.3,2x))")(Hloc(iorb,jorb),jorb =1,Nso)
    enddo
    write(unit,*)""
    if(present(file))close(unit)
  end subroutine print_Hloc_so



END MODULE ED_HLOC_DECOMPOSITION
