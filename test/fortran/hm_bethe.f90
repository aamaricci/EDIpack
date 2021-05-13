program hm_bethe
  USE EDIPACK
  USE SCIFOR
  USE MPI
  implicit none
  !INPUT:
  character(len=16)                           :: finput
  integer                                     :: Le
  real(8),dimension(5)                        :: Wbethe,Dbethe
  logical                                     :: mixG0
  real(8)                                     :: wmixing
  !
  !Bath:
  integer                                     :: Nb
  real(8),dimension(:),allocatable            :: Bath,Bath_prev
  !Dimension:
  integer                                     :: Nso
  !Local Variables:
  real(8),dimension(:,:,:,:),allocatable      :: Hloc
  real(8),dimension(:),allocatable            :: H0
  real(8),dimension(:,:),allocatable          :: Dbands
  real(8),dimension(:,:),allocatable          :: Ebands
  real(8),dimension(:),allocatable            :: Wband
  real(8),dimension(:),allocatable            :: de
  !Local Functions:  
  complex(8),allocatable,dimension(:,:,:,:,:) :: Smats,Sreal
  complex(8),allocatable,dimension(:,:,:,:,:) :: Gmats,Greal
  complex(8),allocatable,dimension(:,:,:,:,:) :: Delta
  !Other:
  integer                                     :: iloop,iorb,i
  logical                                     :: converged
  complex(8)                                  :: zeta
  complex(8),allocatable,dimension(:)         :: Gtest
  real(8),dimension(:),allocatable            :: wm,wr
  integer                                     :: comm,rank
  logical                                     :: master


  !MPI Wrapper from SciFortran:
  call init_MPI()
  comm = MPI_COMM_WORLD
  call StartMsg_MPI(comm)
  rank = get_Rank_MPI(comm)
  master = get_Master_MPI(comm)

  !Input Parser from SciFortran
  call parse_cmd_variable(finput,"FINPUT",default='inputED.conf')
  call parse_input_variable(Le,"LE",finput,default=1000)
  call parse_input_variable(Wbethe,"WBETHE",finput,default=[1d0,1d0,1d0,1d0,1d0])
  call parse_input_variable(Dbethe,"DBETHE",finput,default=[0d0,0d0,0d0,0d0,0d0])
  call parse_input_variable(wmixing,"WMIXING",finput,default=0.5d0)
  call parse_input_variable(mixG0,"mixG0",finput,default=.false.)
  !
  !READ ED INPUT:
  call ed_read_input(trim(finput),comm)

  Nso=Nspin*Norb

  !Build Density of States:
  !linspace, dens_bethe, arange, pi in SciFortran
  allocate(Ebands(Nso,Le))
  allocate(Dbands(Nso,Le))
  allocate(Wband(Nso))
  allocate(H0(Nso))
  allocate(de(Nso))
  Wband = Wbethe(:Norb)
  H0    = Dbethe(:Norb)
  do iorb=1,Norb
     Ebands(iorb,:) = linspace(-Wband(iorb),Wband(iorb),Le,mesh=de(iorb))
     Dbands(iorb,:) = dens_bethe(Ebands(iorb,:),Wband(iorb))*de(iorb)
  enddo


  !Frequency array:
  allocate(wm(Lmats),wr(Lreal))
  wm = pi/beta*(2*arange(1,Lmats)-1)
  wr = linspace(wini,wfin,Lreal)


  !All functions have shape [Nspin,Nspin,Norb,Norb(,L)]:
  allocate(Delta(Nspin,Nspin,Norb,Norb,Lmats))
  allocate(Smats(Nspin,Nspin,Norb,Norb,Lmats))
  allocate(Gmats(Nspin,Nspin,Norb,Norb,Lmats))
  allocate(Sreal(Nspin,Nspin,Norb,Norb,Lreal))
  allocate(Greal(Nspin,Nspin,Norb,Norb,Lreal))
  allocate(Hloc(Nspin,Nspin,Norb,Norb))
  allocate(Gtest(Lmats))
  Hloc(1,1,:,:)=diag(H0)


  !Setup solver
  Nb=ed_get_bath_dimension()
  allocate(bath(Nb),bath_prev(Nb))
  call ed_init_solver(comm,bath,Hloc)

  !DMFT loop
  iloop=0;converged=.false.
  do while(.not.converged.AND.iloop<nloop)
     iloop=iloop+1
     write(*,*)iloop,nloop,"DMFT-loop"
     !
     !Solve the EFFECTIVE IMPURITY PROBLEM (first w/ a guess for the bath)
     call ed_solve(comm,bath)
     !
     !Get Self-energies
     call ed_get_sigma_matsubara(Smats)
     call ed_get_sigma_realaxis(Sreal)
     !
     !Compute the local gf:
     do iorb=1,Norb
        do i=1,Lmats
           zeta    = xi*wm(i) + xmu - H0(iorb) - Smats(1,1,iorb,iorb,i)
           Gmats(1,1,iorb,iorb,i) = sum(Dbands(iorb,:)/(zeta - Ebands(iorb,:)))
        enddo
        do i=1,Lreal
           zeta     = cmplx(wr(i),eps) + xmu - H0(iorb) - Sreal(1,1,iorb,iorb,i)
           Greal(1,1,iorb,iorb,i) = sum(Dbands(iorb,:)/(zeta - Ebands(iorb,:)))
        enddo
        call splot("Gloc_l"//str(iorb)//"_iw.dat",wm,Gmats(1,1,iorb,iorb,:))    !in SciFortran
        call splot("Gloc_l"//str(iorb)//"_realw.dat",wr,Greal(1,1,iorb,iorb,:)) !in SciFortran
     enddo
     !
     !Get the Delta function and FIT:
     do iorb=1,Norb
        Delta(1,1,iorb,iorb,:) = 0.25d0*Wband(iorb)*Gmats(1,1,iorb,iorb,:)        
     enddo
     call ed_chi2_fitgf(comm,Delta,Bath,Hloc,ispin=1)
     !
     !MIXING:
     if(iloop>1)Bath = wmixing*Bath + (1.d0-wmixing)*Bath_prev
     Bath_prev=Bath
     !
     !Check convergence (if required change chemical potential)
     Gtest=zero
     do iorb=1,Norb
        Gtest=Gtest+Delta(1,1,iorb,iorb,:)/Norb
     enddo
     converged = check_convergence(Gtest,dmft_error,nsuccess,nloop)
     write(*,*)"---------"
  enddo

  call finalize_MPI()


contains

  function check_convergence(Xnew,eps,N1,N2) result(convergence)
    complex(8),intent(in)           :: Xnew(:)
    real(8),intent(in)              :: eps
    integer,intent(in)              :: N1,N2
    logical                         :: convergence  
    integer                         :: i,j,Msum
    real(8)                         :: err
    real(8)                         :: M,S
    complex(8),save,allocatable     :: Xold(:)
    integer,save                    :: success=0,check=1
    character(len=100)              :: file_
    file_='error.err'
    Msum=size(Xnew)
    if(.not.allocated(Xold))then
       allocate(Xold(Msum))
       Xold=0.d0
    endif
    S=0.d0 ; M=0.d0
    do i=1,Msum
       M=M + abs(Xnew(i)-Xold(i))
       S=S + abs(Xnew(i))
    enddo
    err= M/S
    Xold=Xnew
    open(10,file=reg(file_),position="append")
    write(10,"(I5,ES15.7)")check,err
    close(10)
    if(err < eps)then
       success=success+1
    else
       success=0
    endif
    convergence=.false.
    if(success > N1)convergence=.true.
    if(check>=N2)then
       open(10,file="ERROR.README")
       write(10,*)""
       close(10)
       write(*,"(A,I4,A)")"Not converged after",N2," iterations."
    endif
    if(convergence)then
       write(*,"(A,ES15.7)")bold_green("error="),err
    else
       if(err < eps)then
          write(*,"(A,ES15.7)")bold_yellow("error="),err
       else
          write(*,"(A,ES15.7)")bold_red("error="),err
       endif
    endif
    check=check+1
  end function check_convergence


end program hm_bethe

  






  ! !+----------------------------------------+
  ! subroutine get_delta_bethe
  !   integer                     :: i,j,iorb
  !   complex(8)                  :: iw,zita,g0loc
  !   complex(8),dimension(Lmats) :: gloc,sigma,Tiw
  !   complex(8),dimension(Lreal) :: grloc
  !   real(8)                     :: wm(Lmats),wr(Lreal),tau(0:Lmats),C0,C1,n0
  !   real(8),dimension(0:Lmats)  :: sigt,gtau,Ttau
  !   real(8),dimension(3)  :: Scoeff
  !   wm = pi/beta*(2*arange(1,Lmats)-1)
  !   wr = linspace(wini,wfin,Lreal)
  !      do i=1,Lmats
  !         iw = xi*wm(i)
  !         zita    = iw + xmu - impSmats(1,1,1,1,i)
  !         gloc(i) = gfbethe(wm(i),zita,Wband)
  !         if(cg_scheme=='weiss')then
  !            delta(i)= one/(one/gloc(i) + impSmats(1,1,1,1,i))
  !         else
  !            delta(i)= iw + xmu - impSmats(1,1,1,1,i) - one/gloc(i)
  !         endif
  !      enddo
  !      do i=1,Lreal
  !         iw=cmplx(wr(i),eps)
  !         zita     = iw + xmu - impSreal(1,1,1,1,i)
  !         grloc(i) = gfbether(wr(i),zita,Wband)
  !      enddo
  !      if(ED_MPI_ID==0)then
  !         call splot("Gloc_iw.ed",wm,gloc)
  !         call splot("Delta_iw.ed",wm,delta)
  !         call splot("Gloc_realw.ed",wr,-dimag(grloc)/pi,dreal(grloc))
  !      endif
  !      ! tau(0:) = linspace(0.d0,beta,Lmats+1)
  !      ! C0=Uloc(1)*(ed_dens_up(1)-0.5d0)
  !      ! C1=Uloc(1)**2*ed_dens_up(1)*(1.d0-ed_dens_dw(1))
  !      ! Tiw=dcmplx(C0,-C1/wm)
  !      ! call splot("Tail_iw.ed",wm,Tiw)
  !      ! Ttau = -C1/2.d0
  !      ! Sigma = impSmats(1,1,1,1,:)  - Tiw
  !      ! call fftgf_iw2tau(Sigma,Sigt(0:),beta,notail=.true.)
  !      ! Sigt=Sigt + Ttau
  !      ! call splot("Sigma_tau.ed",tau,sigt)
  !      ! Sigt=Sigt !+ Ttau
  !      ! call fftgf_tau2iw(sigt(0:),sigma,beta)
  !      ! Sigma=Sigma !+ Tiw
  !      ! call splot("Sigma_iw.ed",wm,sigma)
  !      ! call fftgf_iw2tau(gloc,gtau(0:),beta)
  !      ! call splot("Gloc_tau.ed",tau(0:),gtau(0:))
  ! end subroutine get_delta_bethe
  ! !+----------------------------------------+






  ! !<DEBUG:
  ! ! subroutine get_ed_energy(Lk) 
  ! !   integer               :: Lk
  ! !   real(8),dimension(Lk) :: ek
  ! !   real(8)               :: de
  ! !   real(8),dimension(Lk) :: Wtk
  ! !   ek  = linspace(-Wband,Wband,Lk,mesh=de)
  ! !   Wtk = dens_bethe(ek,wband)*de
  ! !   call ed_kinetic_energy(impSmats(1,1,1,1,:),ek,wtk)
  ! ! end subroutine get_ed_energy


  ! function get_energy(Lk) result(H0)
  !   integer                     :: Lk
  !   complex(8),dimension(Lk)    :: Hk
  !   complex(8),dimension(Lmats) :: Sigma
  !   real(8)                     :: H0
  !   real(8),dimension(Lk)       :: Wtk
  !   real(8)                     :: Tail0,Tail1
  !   real(8)                     :: Sigma_HF,Ak,Bk
  !   complex(8)                  :: Ck,Dk,Zk
  !   complex(8)                  :: Zeta,Gk,Tk
  !   integer                     :: i,ik,iorb
  !   real(8),dimension(Lmats)    :: wm
  !   real(8)                     :: de
  !   !
  !   wm = pi/beta*dble(2*arange(1,Lmats)-1)
  !   !
  !   Hk  = one*linspace(-Wband,Wband,Lk,mesh=de)
  !   Wtk = dens_bethe(dreal(Hk(:)),wband)*de
  !   Sigma = impSmats(1,1,1,1,:) 
  !   Sigma_HF = dreal(Sigma(Lmats))
  !   !
  !   H0=0.d0
  !   do ik=1,Lk
  !      Ak = Hk(ik)
  !      Bk =-Hk(ik) - Sigma_hf
  !      do i=1,Lmats
  !         Gk = one/(xi*wm(i) + xmu - Hk(ik) - Sigma(i) )
  !         Tk = one/(xi*wm(i)) - Bk/(xi*wm(i))**2
  !         Ck = Ak*(Gk - Tk)
  !         H0 = H0 + Ck*Wtk(ik)
  !      enddo
  !   enddo
  !   H0=H0/beta*4d0
  !   !
  !   Tail0=zero
  !   Tail1=zero
  !   do ik=1,Lk
  !      Ak= Hk(ik)
  !      Bk =-Hk(ik) - Sigma_hf
  !      Ck= Ak*Bk
  !      Tail0 = Tail0 + 0.5d0*Ak*Wtk(ik)
  !      Tail1 = Tail1 + 0.25d0*Ck*Wtk(ik)*beta
  !   enddo
  !   Tail0=2d0*Tail0
  !   Tail1=2d0*Tail1
  !   H0 = H0 + Tail0 + Tail1
  ! end function get_energy
  ! !>DEBUG




