  subroutine ed_get_density_matrix_single(dm_,custom_rot,dm_eig_,dm_rot_)
    implicit none
    !passed
    complex(8),allocatable,intent(out)           :: dm_(:,:)
    complex(8),allocatable,intent(in) ,optional  :: custom_rot(:,:)
    real(8),allocatable,intent(out)   ,optional  :: dm_eig_(:)
    complex(8),allocatable,intent(out),optional  :: dm_rot_(:,:)
    !internal
    integer                                      :: unit
    integer                                      :: iorb,jorb,ispin,jspin,io,jo
    complex(8)                                   :: Tr
    complex(8),allocatable                       :: dm_custom_rot(:,:)
    real(8)                                      :: soc
    !
    if (.not.allocated(imp_density_matrix)) then
       write(LOGfile,"(A)") "imp_density_matrix is not allocated"
       stop
    endif
    !
    if(allocated(dm_))                         deallocate(dm_)          ;allocate(dm_(Nspin*Norb,Nspin*Norb))          ;dm_ = zero
    if(allocated(dm_custom_rot))               deallocate(dm_custom_rot);allocate(dm_custom_rot(Nspin*Norb,Nspin*Norb));dm_custom_rot = zero
    if(present(dm_eig_).and.allocated(dm_eig_))deallocate(dm_eig_)      ;allocate(dm_eig_(Nspin*Norb))                 ;dm_eig_ = 0.0d0
    if(present(dm_rot_).and.allocated(dm_rot_))deallocate(dm_rot_)      ;allocate(dm_rot_(Nspin*Norb,Nspin*Norb))      ;dm_rot_ = zero
    !
    ! dm in the impurity problem basis
    dm_ = nn2so_reshape(imp_density_matrix,Nspin,Norb)
    !
    if(bath_type=="replica")then
       !
       ! dm in her diagonal basis
       if(present(dm_eig_).and.present(dm_rot_))then
          dm_rot_=dm_
          call eigh(dm_rot_,dm_eig_,jobz='V',uplo='U')
       endif
       !
       ! dm in the basis defined by custom_rot
       dm_custom_rot=matmul(transpose(conjg(custom_rot)),matmul(dm_,custom_rot))
       !
    elseif(bath_type=="normal")then !.and.SOC/=0.d0
       !
       ! here I assume that custom_rot is: {J,jz}-->{t2g,Sz} / {Lz,Sz}
       dm_custom_rot = nn2so_reshape(imp_density_matrix,Nspin,Norb)
       dm_=matmul(custom_rot,matmul(dm_custom_rot,transpose(conjg(custom_rot))))
       !
       ! dm in her diagonal basis
       if(present(dm_eig_).and.present(dm_rot_))then
          dm_rot_=dm_
          call eigh(dm_rot_,dm_eig_,jobz='V',uplo='U')
       endif
       !
    endif
    !
    call print_dm(dm_,dm_rot_,dm_eig_,dm_custom_rot,1)
    !
  end subroutine ed_get_density_matrix_single



  subroutine ed_get_density_matrix_lattice(dm_,custom_rot,dm_eig_,dm_rot_)
    implicit none
    !passed
    complex(8),allocatable,intent(out)           :: dm_(:,:,:)
    complex(8),allocatable,intent(in) ,optional  :: custom_rot(:,:)
    real(8),allocatable,intent(out)   ,optional  :: dm_eig_(:,:)
    complex(8),allocatable,intent(out),optional  :: dm_rot_(:,:,:)
    !internal
    integer                                      :: unit
    integer                                      :: iorb,jorb,ispin,jspin,io,jo,ilat,Nlat
    complex(8)                                   :: Tr
    complex(8),allocatable                       :: dm_custom_rot(:,:,:)
    real(8)                                      :: soc
    complex(8),allocatable                       :: dm_tmp(:,:)
    real(8),allocatable                          :: dm_eig_tmp(:)
    complex(8),allocatable                       :: dm_rot_tmp(:,:)
    complex(8),allocatable                       :: dm_custom_rot_tmp(:,:)
    !
    Nlat=size(imp_density_matrix_ii,1)
    !
    if (.not.allocated(imp_density_matrix)) then
       write(LOGfile,"(A)") "imp_density_matrix is not allocated"
       stop
    endif
    !
    if(allocated(dm_))deallocate(dm_)                             ;allocate(dm_(Nlat,Nspin*Norb,Nspin*Norb))          ;dm_ = zero
    if(allocated(dm_custom_rot))deallocate(dm_custom_rot)         ;allocate(dm_custom_rot(Nlat,Nspin*Norb,Nspin*Norb));dm_custom_rot = zero
    if(present(dm_eig_).and.allocated(dm_eig_))deallocate(dm_eig_);allocate(dm_eig_(Nlat,Nspin*Norb))                 ;dm_eig_ = 0.0d0
    if(present(dm_rot_).and.allocated(dm_rot_))deallocate(dm_rot_);allocate(dm_rot_(Nlat,Nspin*Norb,Nspin*Norb))      ;dm_rot_ = zero
    !
    if(allocated(dm_tmp))deallocate(dm_)                          ;allocate(dm_tmp(Nspin*Norb,Nspin*Norb))            ;dm_tmp = zero
    if(allocated(dm_custom_rot_tmp))deallocate(dm_custom_rot_tmp) ;allocate(dm_custom_rot_tmp(Nspin*Norb,Nspin*Norb)) ;dm_custom_rot_tmp = zero
    if(allocated(dm_eig_tmp))deallocate(dm_eig_tmp)               ;allocate(dm_eig_tmp(Nspin*Norb))                   ;dm_eig_tmp = 0.0d0
    if(allocated(dm_rot_tmp))deallocate(dm_rot_tmp)               ;allocate(dm_rot_tmp(Nspin*Norb,Nspin*Norb))        ;dm_rot_tmp = zero
    !
    do ilat=1,Nlat
       !
       ! dm in the impurity problem basis
       dm_(ilat,:,:) = nn2so_reshape(imp_density_matrix_ii(ilat,:,:,:,:),Nspin,Norb)
       !
       if(bath_type=="replica")then
          !
          ! dm in her diagonal basis
          if(present(dm_eig_).and.present(dm_rot_))then
             dm_rot_(ilat,:,:)=dm_(ilat,:,:)
             call eigh(dm_rot_(ilat,:,:),dm_eig_(ilat,:),jobz='V',uplo='U')
          endif
          !
          ! dm in the basis defined by custom_rot
          dm_custom_rot(ilat,:,:)=matmul(transpose(conjg(custom_rot)),matmul(dm_(ilat,:,:),custom_rot))
          !
       elseif(bath_type=="normal")then !.and.SOC/=0.d0
          !
          ! here I assume that custom_rot is: {J,jz}-->{t2g,Sz} / {Lz,Sz}
          dm_custom_rot(ilat,:,:)=nn2so_reshape(imp_density_matrix_ii(ilat,:,:,:,:),Nspin,Norb)
          dm_(ilat,:,:)=matmul(custom_rot,matmul(dm_custom_rot(ilat,:,:),transpose(conjg(custom_rot))))
          !
          ! dm in her diagonal basis
          if(present(dm_eig_).and.present(dm_rot_))then
             dm_rot_(ilat,:,:)=dm_(ilat,:,:)
             call eigh(dm_rot_(ilat,:,:),dm_eig_(ilat,:),jobz='V',uplo='U')
          endif
          !
       endif
       !
       dm_tmp            = zero ; dm_tmp            = dm_(ilat,:,:)
       dm_rot_tmp        = zero ; dm_rot_tmp        = dm_rot_(ilat,:,:)
       dm_eig_tmp        = zero ; dm_eig_tmp        = dm_eig_(ilat,:)
       dm_custom_rot_tmp = zero ; dm_custom_rot_tmp = dm_custom_rot(ilat,:,:)
       call print_dm(dm_tmp,dm_rot_tmp,dm_eig_tmp,dm_custom_rot_tmp,ilat)
       !
    enddo
  end subroutine ed_get_density_matrix_lattice


  subroutine print_dm(dm_,dm_rot_,dm_eig_,dm_custom_rot,ndx)
    implicit none
    integer               ,intent(in)            :: ndx
    complex(8),allocatable,intent(in)            :: dm_(:,:)
    complex(8),allocatable,intent(in)            :: dm_custom_rot(:,:)
    real(8),allocatable   ,intent(in),optional   :: dm_eig_(:)
    complex(8),allocatable,intent(in),optional   :: dm_rot_(:,:)
    !internal
    integer                                      :: unit
    character(len=24)                            :: suffix
    integer                                      :: iorb,jorb,ispin,jspin,io,jo
    !
    suffix="imp_density_matrix_"//reg(str(ndx))//".dat"
    !
    unit = free_unit()
    open(unit,file=suffix,action="write",position="rewind",status='unknown')
    !
    write(unit,"(A90)")"# density matrix in the impurity problem basis REAL part:"
    do io=1,Nspin*Norb
       write(unit,"(90(F15.9,1X))") (real(dm_(io,jo)),jo=1,Nspin*Norb)
    enddo
    write(unit,*)
    !
    write(unit,"(A90)")"# density matrix in the impurity problem basis IMAGINARY part:"
    do io=1,Nspin*Norb
       write(unit,"(90(F15.9,1X))") (aimag(dm_(io,jo)),jo=1,Nspin*Norb)
    enddo
    write(unit,*)
    !
    if(present(dm_eig_).and.present(dm_rot_))then
       write(unit,"(A90)")"# eigenvalues of density matrix"
       write(unit,'(10F22.12)') dm_eig_
       write(unit,*)
       !
       write(unit,"(A90)")"# density matrix eigenvector matrix REAL part:"
       do io=1,Nspin*Norb
          write(unit,"(90(F15.9,1X))") (real(dm_rot_(io,jo)),jo=1,Nspin*Norb)
       enddo
       write(unit,*)
       !
       write(unit,"(A90)")"# density matrix eigenvector matrix IMAGINARY part:"
       do io=1,Nspin*Norb
          write(unit,"(90(F15.9,1X))") (aimag(dm_rot_(io,jo)),jo=1,Nspin*Norb)
       enddo
       write(unit,*)
    endif
    !
    write(unit,"(A90)")"# density matrix in the basis defined by custom_rot REAL part:"
    do io=1,Nspin*Norb
       write(unit,"(90(F15.9,1X))") (real(dm_custom_rot(io,jo)),jo=1,Nspin*Norb)
    enddo
    write(unit,*)
    !
    write(unit,"(A90)")"# density matrix in the basis defined by custom_rot IMAGINARY part:"
    do io=1,Nspin*Norb
       write(unit,"(90(F15.9,1X))") (aimag(dm_custom_rot(io,jo)),jo=1,Nspin*Norb)
    enddo
    write(unit,*)
    write(unit,"(A30)")"# J basis densities"
    write(unit,"(90(F15.9,1X))") (real(dm_custom_rot(io,io)),io=1,Nspin*Norb)
    !
    close(unit)
    !
  end subroutine print_dm
