program hm_bethe
  USE EDIPACK
  USE SCIFOR
  USE MPI
  implicit none
  !INPUT:
  integer                                     :: Le
  real(8)                                     :: wmixing
  real(8)                                     :: Wband
  !
  !Bath:
  integer                                     :: Nb
  real(8),dimension(:),allocatable            :: Bath,Bath_prev
  !Local Variables:
  real(8),dimension(:,:,:,:),allocatable      :: Hloc
  real(8),dimension(:),allocatable            :: Ebands
  real(8),dimension(:),allocatable            :: Dbands

  real(8)                                     :: de
  !Local Functions:  
  complex(8),allocatable,dimension(:,:,:,:,:) :: Smats,Sreal
  complex(8),allocatable,dimension(:,:,:,:,:) :: Gmats,Greal
  complex(8),allocatable,dimension(:,:,:,:,:) :: Delta
  !Other:
  integer                                     :: iloop,iorb,i
  logical                                     :: converged
  complex(8)                                  :: zeta
  real(8),dimension(:),allocatable            :: wm,wr
  integer                                     :: comm,rank
  logical                                     :: master


  !MPI Wrapper from SciFortran:
  call init_MPI()
  comm = MPI_COMM_WORLD
  call StartMsg_MPI(comm)
  rank = get_Rank_MPI(comm)
  master = get_Master_MPI(comm)


  !READ ED INPUT:
  call ed_read_input('inputED.conf')
  if(Nspin/=1.OR.Norb/=1)stop "This test code is for Nspin=1 + Norb=1."
  Le     = 1000
  wmixing= 0.5d0
  Wband  = 1d0

  !BUILD Density of States:
  !linspace, dens_bethe, arange, pi in SciFortran
  allocate(Ebands(Le))
  allocate(Dbands(Le))
  Ebands = linspace(-Wband,Wband,Le,mesh=de)
  Dbands = dens_bethe(Ebands,Wband)*de

  !BUILD frequency array:
  allocate(wm(Lmats),wr(Lreal))
  wm = pi/beta*(2*arange(1,Lmats)-1)
  wr = linspace(wini,wfin,Lreal)


  !ALL functions must have shape [Nspin,Nspin,Norb,Norb(,L)]:
  allocate(Delta(Nspin,Nspin,Norb,Norb,Lmats))
  allocate(Smats(Nspin,Nspin,Norb,Norb,Lmats))
  allocate(Gmats(Nspin,Nspin,Norb,Norb,Lmats))
  allocate(Sreal(Nspin,Nspin,Norb,Norb,Lreal))
  allocate(Greal(Nspin,Nspin,Norb,Norb,Lreal))
  allocate(Hloc(Nspin,Nspin,Norb,Norb))
  Hloc = 0d0


  !SETUP SOLVER
  Nb=ed_get_bath_dimension()
  allocate(bath(Nb),bath_prev(Nb))
  call ed_init_solver(bath)

  !DMFT CYCLE
  iloop=0;converged=.false.
  do while(.not.converged.AND.iloop<nloop)
     iloop=iloop+1
     if(master)write(*,*)"DMFT-loop:",iloop,"/",nloop
     !
     !Solve the EFFECTIVE IMPURITY PROBLEM (first w/ a guess for the bath)
     call ed_solve(bath,Hloc)
     !
     !Get Self-energies
     call ed_get_sigma_matsubara(Smats)
     call ed_get_sigma_realaxis(Sreal)
     !
     !Compute the local gf:
     do i=1,Lmats
        zeta    = xi*wm(i) + xmu  - Smats(1,1,1,1,i)
        Gmats(1,1,1,1,i) = sum(Dbands(:)/(zeta - Ebands(:)))
     enddo
     do i=1,Lreal
        zeta     = cmplx(wr(i),eps) + xmu - Sreal(1,1,1,1,i)
        Greal(1,1,1,1,i) = sum(Dbands(:)/(zeta - Ebands(:)))
     enddo
     if(master)call splot("Gloc_iw.dat",wm,Gmats(1,1,1,1,:))    !splot is in SciFortran
     if(master)call splot("Gloc_realw.dat",wr,Greal(1,1,1,1,:)) !splot is in SciFortran
     !
     !Get the Delta function and FIT:
     Delta(1,1,1,1,:) = 0.25d0*Wband*Gmats(1,1,1,1,:)        
     call ed_chi2_fitgf(Delta,Bath,ispin=1)
     !
     !MIXING:
     if(iloop>1)Bath = wmixing*Bath + (1.d0-wmixing)*Bath_prev
     Bath_prev=Bath
     !
     !Check convergence
     converged = check_convergence(Delta(1,1,1,1,:),dmft_error,nsuccess,nloop)
     if(master)write(*,*)"---------"
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
    if(master)then
       open(10,file=reg(file_),position="append")
       write(10,"(I5,ES15.7)")check,err
       close(10)
    endif
    if(err < eps)then
       success=success+1
    else
       success=0
    endif
    convergence=.false.
    if(success > N1)convergence=.true.
    if(check>=N2)then
       if(master)then
          open(10,file="ERROR.README")
          write(10,*)""
          close(10)
          write(*,"(A,I4,A)")"Not converged after",N2," iterations."
       endif
    endif
    if(convergence)then
       if(master)write(*,"(A,ES15.7)")bold_green("error="),err
    else
       if(err < eps)then
          if(master)write(*,"(A,ES15.7)")bold_yellow("error="),err
       else
          if(master)write(*,"(A,ES15.7)")bold_red("error="),err
       endif
    endif
    check=check+1
  end function check_convergence


end program hm_bethe












