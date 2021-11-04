MODULE ED_IO
  USE ED_INPUT_VARS
  USE ED_VARS_GLOBAL
  USE ED_AUX_FUNX
  USE ED_BATH
  USE ED_BATH_FUNCTIONS
  USE SF_LINALG
  USE SF_ARRAYS, only: linspace,arange
  USE SF_IOTOOLS, only: str,reg,free_unit,splot,sread
  implicit none
  private


  !Retrieve imp GF through routines.
  interface ed_get_gimp_matsubara
     module procedure ed_get_gimp_matsubara_main
     module procedure ed_get_gimp_matsubara_lattice
  end interface ed_get_gimp_matsubara


  interface ed_get_gimp_realaxis
     module procedure ed_get_gimp_realaxis_main
     module procedure ed_get_gimp_realaxis_lattice
  end interface ed_get_gimp_realaxis


  !Retrieve self-energy through routines:
  interface ed_get_sigma_matsubara
     module procedure ed_get_sigma_matsubara_main
     module procedure ed_get_sigma_matsubara_lattice
  end interface ed_get_sigma_matsubara

  interface ed_get_sigma_realaxis
     module procedure ed_get_sigma_realaxis_main
     module procedure ed_get_sigma_realaxis_lattice
  end interface ed_get_sigma_realaxis


  !Retrieve static common observables  
  interface ed_get_dens
     module procedure ed_get_dens_main
     module procedure ed_get_dens_lattice
  end interface ed_get_dens

  interface ed_get_mag
     module procedure ed_get_mag_main
     module procedure ed_get_mag_lattice
  end interface ed_get_mag

  interface ed_get_docc
     module procedure ed_get_docc_main
     module procedure ed_get_docc_lattice
  end interface ed_get_docc

  interface ed_get_eimp
     module procedure :: ed_get_eimp_main
     module procedure :: ed_get_eimp_lattice
  end interface ed_get_eimp

  interface ed_get_doubles
     module procedure :: ed_get_doubles_main
     module procedure :: ed_get_doubles_lattice
  end interface ed_get_doubles


  public :: ed_get_sigma_matsubara
  public :: ed_get_sigma_realaxis

  public :: ed_get_gimp_matsubara
  public :: ed_get_gimp_realaxis

  public :: ed_get_dens
  public :: ed_get_mag
  public :: ed_get_docc
  public :: ed_get_eimp
  public :: ed_get_doubles


  !****************************************************************************************!
  !****************************************************************************************!

  public :: ed_print_impSigma
  public :: ed_print_impG
  public :: ed_print_impG0
  public :: ed_print_impD
  public :: ed_print_impChi


  !****************************************************************************************!
  !****************************************************************************************!



  !Frequency and time arrays:
  !=========================================================
  ! real(8),dimension(:),allocatable :: wm,tau,wr,vm
  character(len=64)                :: suffix





contains



  !+--------------------------------------------------------------------------+!
  ! PURPOSE: Retrieve measured values of the impurity GF and self-energy 
  !+--------------------------------------------------------------------------+!
  include "get_gimp.f90"




  !+--------------------------------------------------------------------------+!
  ! PURPOSE: Retrieve measured values of the local observables
  !+--------------------------------------------------------------------------+!
  include "get_obs.f90"







  !+------------------------------------------------------------------+
  !                         PRINT SIGMA:
  !+------------------------------------------------------------------+  
  subroutine ed_print_impSigma
    integer                                           :: i,ispin,isign,unit(2),iorb,jorb
    character(len=20)                                 :: suffix
    integer,dimension(:),allocatable                  :: getIorb,getJorb
    integer                                           :: totNorb,l
    !
    ! if(.not.allocated(wm))allocate(wm(Lmats))
    ! if(.not.allocated(wr))allocate(wr(Lreal))
    ! wm     = pi/beta*real(2*arange(1,Lmats)-1,8)
    ! wr     = linspace(wini,wfin,Lreal)
    call allocate_grids()
    !
    select case(bath_type)
    case default                !Diagonal in both spin and orbital
       totNorb=Norb
       allocate(getIorb(totNorb),getJorb(totNorb))
       l=0
       do iorb=1,Norb
          L=l+1
          getIorb(l)=iorb
          getJorb(l)=iorb
       enddo
       totNorb=l
    case ('hybrid')             !Diagonal in spin only. Full Orbital structure
       totNorb=Norb*(Norb+1)/2
       allocate(getIorb(totNorb),getJorb(totNorb))
       l=0
       do iorb=1,Norb
          do jorb=iorb,Norb
             l=l+1
             getIorb(l)=iorb
             getJorb(l)=jorb
          enddo
       enddo
    end select
    if(l/=totNorb)stop "print_gf_normal error counting the orbitals"
    !!
    !Print the impurity functions:
    do ispin=1,Nspin
       do l=1,totNorb
          iorb=getIorb(l)
          jorb=getJorb(l)
          suffix="_l"//str(iorb)//str(jorb)//"_s"//str(ispin)
          call splot("impSigma"//reg(suffix)//"_iw"//reg(ed_file_suffix)//".ed"   ,wm,impSmats(ispin,ispin,iorb,jorb,:))
          call splot("impSigma"//reg(suffix)//"_realw"//reg(ed_file_suffix)//".ed",wr,impSreal(ispin,ispin,iorb,jorb,:))
       enddo
    enddo
    !
    ! if(allocated(wm))deallocate(wm)
    ! if(allocated(wr))deallocate(wr)
    call deallocate_grids()
    !
  end subroutine ed_print_impSigma




  !+------------------------------------------------------------------+
  !                         PRINT G
  !+------------------------------------------------------------------+  
  subroutine ed_print_impG
    integer                                           :: i,ispin,isign,unit(2),iorb,jorb
    character(len=20)                                 :: suffix
    integer,dimension(:),allocatable                  :: getIorb,getJorb
    integer                                           :: totNorb,l
    !
    call allocate_grids()
    !
    select case(bath_type)
    case default                !Diagonal in both spin and orbital
       totNorb=Norb
       allocate(getIorb(totNorb),getJorb(totNorb))
       l=0
       do iorb=1,Norb
          L=l+1
          getIorb(l)=iorb
          getJorb(l)=iorb
       enddo
       totNorb=l
    case ('hybrid')             !Diagonal in spin only. Full Orbital structure
       totNorb=Norb*(Norb+1)/2
       allocate(getIorb(totNorb),getJorb(totNorb))
       l=0
       do iorb=1,Norb
          do jorb=iorb,Norb
             l=l+1
             getIorb(l)=iorb
             getJorb(l)=jorb
          enddo
       enddo
    end select
    if(l/=totNorb)stop "print_gf_normal error counting the orbitals"
    !!
    !Print the impurity functions:
    do ispin=1,Nspin
       do l=1,totNorb
          iorb=getIorb(l)
          jorb=getJorb(l)
          suffix="_l"//str(iorb)//str(jorb)//"_s"//str(ispin)
          call splot("impG"//reg(suffix)//"_iw"//reg(ed_file_suffix)//".ed"   ,wm,impGmats(ispin,ispin,iorb,jorb,:))
          call splot("impG"//reg(suffix)//"_realw"//reg(ed_file_suffix)//".ed",wr,impGreal(ispin,ispin,iorb,jorb,:))
       enddo
    enddo
    !
    call deallocate_grids()
    !
  end subroutine ed_print_impG



  !+------------------------------------------------------------------+
  !                         PRINT G0
  !+------------------------------------------------------------------+  
  subroutine ed_print_impG0
    integer                                           :: i,ispin,isign,unit(2),iorb,jorb
    character(len=20)                                 :: suffix
    integer,dimension(:),allocatable                  :: getIorb,getJorb
    integer                                           :: totNorb,l
    !
    call allocate_grids()
    !
    select case(bath_type)
    case default                !Diagonal in both spin and orbital
       totNorb=Norb
       allocate(getIorb(totNorb),getJorb(totNorb))
       l=0
       do iorb=1,Norb
          L=l+1
          getIorb(l)=iorb
          getJorb(l)=iorb
       enddo
       totNorb=l
    case ('hybrid')             !Diagonal in spin only. Full Orbital structure
       totNorb=Norb*(Norb+1)/2
       allocate(getIorb(totNorb),getJorb(totNorb))
       l=0
       do iorb=1,Norb
          do jorb=iorb,Norb
             l=l+1
             getIorb(l)=iorb
             getJorb(l)=jorb
          enddo
       enddo
    end select
    if(l/=totNorb)stop "print_gf_normal error counting the orbitals"
    !!
    !Print the impurity functions:
    do ispin=1,Nspin
       do l=1,totNorb
          iorb=getIorb(l)
          jorb=getJorb(l)
          suffix="_l"//str(iorb)//str(jorb)//"_s"//str(ispin)
          call splot("impG0"//reg(suffix)//"_iw"//reg(ed_file_suffix)//".ed"   ,wm,impG0mats(ispin,ispin,iorb,jorb,:))
          call splot("impG0"//reg(suffix)//"_realw"//reg(ed_file_suffix)//".ed",wr,impG0real(ispin,ispin,iorb,jorb,:))
       enddo
    enddo
    !
    call deallocate_grids()
    !
  end subroutine ed_print_impG0



  !+------------------------------------------------------------------+
  !                         PRINT D (phonon Green's function)
  !+------------------------------------------------------------------+  
  subroutine ed_print_impD
    !
    call allocate_grids()
    !
    !Print the impurity functions:
    call splot("impDph_iw.ed"   ,vm,impDmats(:))
    call splot("impDph_realw.ed",vr,impDreal(:))
    !
    call deallocate_grids()
    !
  end subroutine ed_print_impD



  !+------------------------------------------------------------------+
  !                         PRINT CHI:
  !+------------------------------------------------------------------+  
  subroutine ed_print_impChi
    if(chispin_flag)call print_chi_spin
    if(chidens_flag)call print_chi_dens
    if(chipair_flag)call print_chi_pair
    if(chiexct_flag)call print_chi_exct
  end subroutine ed_print_impChi

  !                         SPIN-SPIN
  subroutine print_chi_spin
    integer                               :: i,j,iorb,jorb
    call allocate_grids()
    do iorb=1,Norb
       do jorb=1,Norb
          call splot("spinChi_l"//str(iorb)//str(jorb)//"_tau"//reg(ed_file_suffix)//".ed",tau,spinChi_tau(iorb,jorb,0:))
          call splot("spinChi_l"//str(iorb)//str(jorb)//"_realw"//reg(ed_file_suffix)//".ed",vr,spinChi_w(iorb,jorb,:))
          call splot("spinChi_l"//str(iorb)//str(jorb)//"_iw"//reg(ed_file_suffix)//".ed",vm,spinChi_iv(iorb,jorb,:))
       enddo
    enddo
    call deallocate_grids()
  end subroutine print_chi_spin
  !                     DENSITY-DENSITY
  subroutine print_chi_dens
    integer                               :: i,j,iorb,jorb
    call allocate_grids()
    do iorb=1,Norb
       do jorb=1,Norb
          call splot("densChi_l"//str(iorb)//str(jorb)//"_tau"//reg(ed_file_suffix)//".ed",tau,densChi_tau(iorb,jorb,0:))
          call splot("densChi_l"//str(iorb)//str(jorb)//"_realw"//reg(ed_file_suffix)//".ed",vr,densChi_w(iorb,jorb,:))
          call splot("densChi_l"//str(iorb)//str(jorb)//"_iw"//reg(ed_file_suffix)//".ed",vm,densChi_iv(iorb,jorb,:))
       enddo
    enddo
    call deallocate_grids()
  end subroutine print_chi_dens
  !                     PAIR-PAIR
  subroutine print_chi_pair
    integer                               :: i,j,iorb,jorb
    call allocate_grids()
    do iorb=1,Norb
       do jorb=1,Norb
          call splot("pairChi_l"//str(iorb)//str(jorb)//"_tau"//reg(ed_file_suffix)//".ed",tau,pairChi_tau(iorb,jorb,0:))
          call splot("pairChi_l"//str(iorb)//str(jorb)//"_realw"//reg(ed_file_suffix)//".ed",vr,pairChi_w(iorb,jorb,:))
          call splot("pairChi_l"//str(iorb)//str(jorb)//"_iw"//reg(ed_file_suffix)//".ed",vm,pairChi_iv(iorb,jorb,:))
       enddo
    enddo
    call deallocate_grids()
  end subroutine print_chi_pair
  !                     EXCITON
  subroutine print_chi_exct
    integer                               :: i,j,iorb,jorb
    call allocate_grids()
    do iorb=1,Norb
       do jorb=iorb+1,Norb
          call splot("exctChi_singlet_l"//str(iorb)//str(jorb)//"_tau"//reg(ed_file_suffix)//".ed",tau,exctChi_tau(0,iorb,jorb,0:))
          call splot("exctChi_singlet_l"//str(iorb)//str(jorb)//"_realw"//reg(ed_file_suffix)//".ed",vr,exctChi_w(0,iorb,jorb,:))
          call splot("exctChi_singlet_l"//str(iorb)//str(jorb)//"_iw"//reg(ed_file_suffix)//".ed",vm,exctChi_iv(0,iorb,jorb,:))
       enddo
    enddo
    !
    do iorb=1,Norb
       do jorb=iorb+1,Norb
          call splot("exctChi_tripletXY_l"//str(iorb)//str(jorb)//"_tau"//reg(ed_file_suffix)//".ed",tau,exctChi_tau(1,iorb,jorb,0:))
          call splot("exctChi_tripletXY_l"//str(iorb)//str(jorb)//"_realw"//reg(ed_file_suffix)//".ed",vr,exctChi_w(1,iorb,jorb,:))
          call splot("exctChi_tripletXY_l"//str(iorb)//str(jorb)//"_iw"//reg(ed_file_suffix)//".ed",vm,exctChi_iv(1,iorb,jorb,:))
       enddo
    enddo
    !
    do iorb=1,Norb
       do jorb=iorb+1,Norb
          call splot("exctChi_tripletZ_l"//str(iorb)//str(jorb)//"_tau"//reg(ed_file_suffix)//".ed",tau,exctChi_tau(2,iorb,jorb,0:))
          call splot("exctChi_tripletZ_l"//str(iorb)//str(jorb)//"_realw"//reg(ed_file_suffix)//".ed",vr,exctChi_w(2,iorb,jorb,:))
          call splot("exctChi_tripletZ_l"//str(iorb)//str(jorb)//"_iw"//reg(ed_file_suffix)//".ed",vm,exctChi_iv(2,iorb,jorb,:))
       enddo
    enddo
    !
    call deallocate_grids()
  end subroutine print_chi_exct









  ! !+--------------------------------------------------------------------------+!
  ! ! PURPOSE: Retrieve Anderson non-interacting green's functions 
  ! !+--------------------------------------------------------------------------+!
  ! function delta_bath_user(x,bath_) result(Delta)
  !   complex(8),dimension(:),intent(in)                  :: x
  !   type(effective_bath)                                :: dmft_bath_
  !   complex(8),dimension(Nspin,Nspin,Norb,Norb,size(x)) :: Delta
  !   real(8),dimension(:)                                :: bath_
  !   logical                                             :: check
  !   check= check_bath_dimension(bath_)
  !   if(.not.check)stop "delta_bath_mats_main_ error: wrong bath dimensions"
  !   call allocate_dmft_bath(dmft_bath_)
  !   call set_dmft_bath(bath_,dmft_bath_)
  !   Delta = delta_bath_function(x,dmft_bath_)
  !   call deallocate_dmft_bath(dmft_bath_)
  ! end function delta_bath_user


  ! function g0and_bath_user(x,bath_) result(G0and)
  !   complex(8),dimension(:),intent(in)                  :: x
  !   type(effective_bath)                                :: dmft_bath_
  !   complex(8),dimension(Nspin,Nspin,Norb,Norb,size(x)) :: G0and
  !   real(8),dimension(:)                                :: bath_
  !   logical                                             :: check
  !   check= check_bath_dimension(bath_)
  !   if(.not.check)stop "g0and_bath_mats_main_ error: wrong bath dimensions"
  !   call allocate_dmft_bath(dmft_bath_)
  !   call set_dmft_bath(bath_,dmft_bath_)
  !   G0and = g0and_bath_function(x,dmft_bath_)
  !   call deallocate_dmft_bath(dmft_bath_)
  ! end function g0and_bath_user



END MODULE ED_IO







