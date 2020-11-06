!##################################################################
! THE CALCULATION OF THE \chi^2 FUNCTIONS USE PROCEDURES FURTHER 
! BELOW TO EVALUATE INDEPENDENTLY THE ANDERSON MODEL:
!  - DELTA, 
!  -\GRAD DELTA
!  - G0
! THE LATTER ARE ADAPTED FROM THE PROCEDURES:
! DELTA_BATH_MATS
! GRAD_DELTA_BATH_MATS
! G0 BATH_MATS
! FOR, YOU NEED TO DECOMPOSE THE a INPUT ARRAY INTO ELEMENTS.
!##################################################################


!+-------------------------------------------------------------+
!PURPOSE  : Chi^2 interface for Irreducible bath normal phase
!+-------------------------------------------------------------+
subroutine chi2_fitgf_normal_normal(fg,bath_,ispin,iorb)
  complex(8),dimension(:,:,:)                 :: fg ![Norb][Norb][Lmats]
  real(8),dimension(:),intent(inout)          :: bath_
  integer                                     :: ispin
  integer,optional                            :: iorb
  real(8),dimension(:),allocatable            :: array_bath
  integer                                     :: iter,stride,jorb,myorb,i,io,j,Asize
  real(8)                                     :: chi
  logical                                     :: check
  type(effective_bath)                        :: dmft_bath
  character(len=256)                          :: suffix
  integer                                     :: unit
  complex(8),dimension(:,:,:,:,:),allocatable :: fgand ![Nspin][][Norb][][Ldelta]
  !
  if(size(fg,1)/=Norb)stop "chi2_fitgf_normal_normal error: size[fg,1]!=Norb"
  if(size(fg,2)/=Norb)stop "chi2_fitgf_normal_normal error: size[fg,2]!=Norb"
  !
  check= check_bath_dimension(bath_)
  if(.not.check)stop "chi2_fitgf_normal_normal error: wrong bath dimensions"
  !
  Ldelta = Lfit ; if(Ldelta>size(fg,3))Ldelta=size(fg,3)
  !
  allocate(Gdelta(1,Ldelta))
  allocate(Xdelta(Ldelta))
  allocate(Wdelta(Ldelta))
  !
  Xdelta = pi/beta*(2*arange(1,Ldelta)-1)
  !
  select case(cg_weight)
  case default
     Wdelta=1d0
  case(2)
     Wdelta=1d0*arange(1,Ldelta)
  case(3)
     Wdelta=Xdelta
  end select
  !
  call allocate_dmft_bath(dmft_bath)
  call set_dmft_bath(bath_,dmft_bath)
  !
  !Asize = get_chi2_bath_size()
  !E_{\s,\a}(:)  [ 1 ][ 1 ][Nbath]
  !V_{\s,\a}(:)  [ 1 ][ 1 ][Nbath]
  Asize = Nbath + Nbath
  allocate(array_bath(Asize))
  !
  do jorb=1,Norb
     if(present(iorb))then
        if(jorb/=iorb)cycle
     endif
     Orb_indx=jorb
     Spin_indx=ispin
     !
     Gdelta(1,1:Ldelta) = fg(jorb,jorb,1:Ldelta)
     !
     !Nbath + Nbath
     stride = 0
     do i=1,Nbath
        io = stride + i
        array_bath(io) = dmft_bath%e(ispin,jorb,i)
     enddo
     stride = Nbath
     do i=1,Nbath
        io = stride + i
        array_bath(io) = dmft_bath%v(ispin,jorb,i)
     enddo
     !
     select case(cg_method)     !0=NR-CG[default]; 1=CG-MINIMIZE; 2=CG+
     case default
        if(cg_grad==0)then
           select case (cg_scheme)
           case ("weiss")
              call fmin_cg(array_bath,chi2_weiss_normal_normal,grad_chi2_weiss_normal_normal,iter,chi,&
                   itmax=cg_niter,&
                   ftol=cg_Ftol,  &
                   istop=cg_stop, &
                   iverbose=(ed_verbose>3))
           case ("delta")
              call fmin_cg(array_bath,chi2_delta_normal_normal,grad_chi2_delta_normal_normal,iter,chi,&
                   itmax=cg_niter,&
                   ftol=cg_Ftol,  &
                   istop=cg_stop, &
                   iverbose=(ed_verbose>3))
           case default
              stop "chi2_fitgf_normal_normal error: cg_scheme != [weiss,delta]"
           end select
        else
           select case (cg_scheme)
           case ("weiss")
              call fmin_cg(array_bath,chi2_weiss_normal_normal,iter,chi,&
                   itmax=cg_niter,&
                   ftol=cg_Ftol,  &
                   istop=cg_stop, &
                   iverbose=(ed_verbose>3))
           case ("delta")
              call fmin_cg(array_bath,chi2_delta_normal_normal,iter,chi,&
                   itmax=cg_niter,&
                   ftol=cg_Ftol,  &
                   istop=cg_stop, &
                   iverbose=(ed_verbose>3))
           case default
              stop "chi2_fitgf_normal_normal error: cg_scheme != [weiss,delta]"
           end select
        endif
        !
        !
     case (1)
        select case (cg_scheme)
        case ("weiss")
           call fmin_cgminimize(array_bath,chi2_weiss_normal_normal,&
                iter,chi,itmax=cg_niter,ftol=cg_Ftol,&
                new_version=cg_minimize_ver,&
                hh_par=cg_minimize_hh,&
                iverbose=(ed_verbose>3))
        case ("delta")
           call fmin_cgminimize(array_bath,chi2_delta_normal_normal,&
                iter,chi,itmax=cg_niter,ftol=cg_Ftol,&
                new_version=cg_minimize_ver,&
                hh_par=cg_minimize_hh,&
                iverbose=(ed_verbose>3))
        case default
           stop "chi2_fitgf_normal_normal error: cg_scheme != [weiss,delta]"
        end select
        !
     end select
     !
     !
     write(LOGfile,"(A,ES18.9,A,I5,A)")&
          "chi^2|iter"//reg(ed_file_suffix)//'= ',chi," | ",iter,&
          "  <--  Orb"//reg(txtfy(jorb))//" Spin"//reg(txtfy(ispin))
     !
     suffix="_orb"//reg(txtfy(jorb))//"_s"//reg(txtfy(ispin))//reg(ed_file_suffix)
     unit=free_unit()
     open(unit,file="chi2fit_results"//reg(suffix)//".ed",position="append")
     write(unit,"(ES18.9,1x,I5)") chi,iter
     close(unit)
     !
     !Nbath + Nbath
     stride = 0
     do i=1,Nbath
        io = stride + i
        dmft_bath%e(ispin,jorb,i) = array_bath(io) 
     enddo
     stride = Nbath
     do i=1,Nbath
        io = stride + i
        dmft_bath%v(ispin,jorb,i) = array_bath(io)
     enddo
     !
  enddo
  !
  call write_dmft_bath(dmft_bath,LOGfile)
  !
  call save_dmft_bath(dmft_bath)
  !
  allocate(fgand(Nspin,Nspin,Norb,Norb,Ldelta))
  if(cg_scheme=='weiss')then
     fgand = g0and_bath_function(xi*Xdelta(:),dmft_bath)
  else
     fgand = delta_bath_function(xi*Xdelta(:),dmft_bath)
  endif
  call write_fit_result(ispin)
  deallocate(fgand)
  !
  call get_dmft_bath(dmft_bath,bath_)
  call deallocate_dmft_bath(dmft_bath)
  deallocate(Gdelta,Xdelta,Wdelta)
  !
contains
  !
  subroutine write_fit_result(ispin)
    integer           :: i,jorb,ispin
    do jorb=1,Norb
       if(present(iorb))then
          if(jorb/=iorb)cycle
       endif
       suffix="_orb"//reg(txtfy(jorb))//"_s"//reg(txtfy(ispin))//reg(ed_file_suffix)
       unit=free_unit()
       if(cg_scheme=='weiss')then
          open(unit,file="fit_weiss"//reg(suffix)//".ed")
       else
          open(unit,file="fit_delta"//reg(suffix)//".ed")
       endif
       do i=1,Ldelta
          write(unit,"(5F24.15)")Xdelta(i),&
               dimag(fg(jorb,jorb,i)),dimag(fgand(ispin,ispin,jorb,jorb,i)),&
               dreal(fg(jorb,jorb,i)),dreal(fgand(ispin,ispin,jorb,jorb,i))
       enddo
       close(unit)
    enddo
  end subroutine write_fit_result
end subroutine chi2_fitgf_normal_normal





! subroutine chi2_fitgf_normal_normal_OneOrb(fg,bath_,ispin,iorb)
!   complex(8),dimension(:,:,:)                 :: fg ![Norb][Norb][Lmats]
!   real(8),dimension(:),intent(inout)          :: bath_
!   integer                                     :: ispin,iorb
!   real(8),dimension(:),allocatable            :: array_bath
!   integer                                     :: iter,stride,i,io,j,Asize
!   real(8)                                     :: chi
!   logical                                     :: check
!   type(effective_bath)                        :: dmft_bath
!   character(len=256)                          :: suffix
!   integer                                     :: unit
!   complex(8),dimension(:,:,:,:,:),allocatable :: fgand ![Nspin][][Norb][][Ldelta]
!   !
!   if(size(fg,1)/=Norb)stop "chi2_fitgf_normal_normal error: size[fg,1]!=Norb"
!   if(size(fg,2)/=Norb)stop "chi2_fitgf_normal_normal error: size[fg,2]!=Norb"
!   !
!   check= check_bath_dimension(bath_)
!   if(.not.check)stop "chi2_fitgf_normal_normal error: wrong bath dimensions"
!   !
!   Ldelta = Lfit ; if(Ldelta>size(fg,3))Ldelta=size(fg,3)
!   !
!   allocate(Gdelta(1,Ldelta))
!   allocate(Xdelta(Ldelta))
!   allocate(Wdelta(Ldelta))
!   !
!   Xdelta = pi/beta*(2*arange(1,Ldelta)-1)
!   !
!   select case(cg_weight)
!   case default
!      Wdelta=1d0
!   case(2)
!      Wdelta=arange(1,Ldelta)
!   case(3)
!      Wdelta=Xdelta
!   end select
!   !
!   call allocate_dmft_bath(dmft_bath)
!   call set_dmft_bath(bath_,dmft_bath)
!   !
!   !Asize = get_chi2_bath_size()
!   !E_{\s,\a}(:)  [ 1 ][ 1 ][Nbath]
!   !V_{\s,\a}(:)  [ 1 ][ 1 ][Nbath]
!   Asize = Nbath + Nbath
!   allocate(array_bath(Asize))
!   !
!   Orb_indx=iorb
!   Spin_indx=ispin
!   !
!   Gdelta(1,1:Ldelta) = fg(iorb,iorb,1:Ldelta)
!   !
!   !Nbath + Nbath
!   stride = 0
!   do i=1,Nbath
!      io = stride + i
!      array_bath(io) = dmft_bath%e(ispin,iorb,i)
!   enddo
!   stride = Nbath
!   do i=1,Nbath
!      io = stride + i
!      array_bath(io) = dmft_bath%v(ispin,iorb,i)
!   enddo
!   !
!   select case(cg_method)     !0=NR-CG[default]; 1=CG-MINIMIZE; 2=CG+
!   case default
!      if(cg_grad==0)then
!         select case (cg_scheme)
!         case ("weiss")
!            call fmin_cg(array_bath,chi2_weiss_normal_normal,grad_chi2_weiss_normal_normal,iter,chi,&
!                 itmax=cg_niter,&
!                 ftol=cg_Ftol,  &
!                 istop=cg_stop, &
!                 iverbose=(ed_verbose>3))
!         case ("delta")
!            call fmin_cg(array_bath,chi2_delta_normal_normal,grad_chi2_delta_normal_normal,iter,chi,&
!                 itmax=cg_niter,&
!                 ftol=cg_Ftol,  &
!                 istop=cg_stop, &
!                 iverbose=(ed_verbose>3))
!         case default
!            stop "chi2_fitgf_normal_normal error: cg_scheme != [weiss,delta]"
!         end select
!      else
!         select case (cg_scheme)
!         case ("weiss")
!            call fmin_cg(array_bath,chi2_weiss_normal_normal,iter,chi,&
!                 itmax=cg_niter,&
!                 ftol=cg_Ftol,  &
!                 istop=cg_stop, &
!                 iverbose=(ed_verbose>3))
!         case ("delta")
!            call fmin_cg(array_bath,chi2_delta_normal_normal,iter,chi,&
!                 itmax=cg_niter,&
!                 ftol=cg_Ftol,  &
!                 istop=cg_stop, &
!                 iverbose=(ed_verbose>3))
!         case default
!            stop "chi2_fitgf_normal_normal error: cg_scheme != [weiss,delta]"
!         end select
!      endif
!      !
!      !
!   case (1)
!      select case (cg_scheme)
!      case ("weiss")
!         call fmin_cgminimize(array_bath,chi2_weiss_normal_normal,iter,chi,&
!              itmax=cg_niter, &
!              ftol=cg_Ftol,&
!              new_version=cg_minimize_ver,&
!              hh_par=cg_minimize_hh,&
!              iverbose=(ed_verbose>3))
!      case ("delta")
!         call fmin_cgminimize(array_bath,chi2_delta_normal_normal,iter,chi,&
!              itmax=cg_niter, &
!              ftol=cg_Ftol,&
!              new_version=cg_minimize_ver,&
!              hh_par=cg_minimize_hh,&
!              iverbose=(ed_verbose>3))
!      case default
!         stop "chi2_fitgf_normal_normal error: cg_scheme != [weiss,delta]"
!      end select
!      !
!   end select
!   !
!   !
!   write(LOGfile,"(A,ES18.9,A,I5,A)")&
!        "chi^2|iter"//reg(ed_file_suffix)//'= ',chi," | ",iter,&
!        "  <--  Orb"//reg(txtfy(iorb))//" Spin"//reg(txtfy(ispin))
!   !
!   suffix="_orb"//reg(txtfy(iorb))//"_s"//reg(txtfy(ispin))//reg(ed_file_suffix)
!   unit=free_unit()
!   open(unit,file="chi2fit_results"//reg(suffix)//".ed",position="append")
!   write(unit,"(ES18.9,1x,I5)") chi,iter
!   close(unit)
!   !
!   !Nbath + Nbath
!   stride = 0
!   do i=1,Nbath
!      io = stride + i
!      dmft_bath%e(ispin,iorb,i) = array_bath(io) 
!   enddo
!   stride = Nbath
!   do i=1,Nbath
!      io = stride + i
!      dmft_bath%v(ispin,iorb,i) = array_bath(io)
!   enddo
!   !
!   call write_dmft_bath(dmft_bath,LOGfile)
!   !
!   call save_dmft_bath(dmft_bath)
!   !
!   allocate(fgand(Nspin,Nspin,Norb,Norb,Ldelta))
!   if(cg_scheme=='weiss')then
!      fgand = g0and_bath_function(xi*Xdelta(:),dmft_bath)
!   else
!      fgand = delta_bath_function(xi*Xdelta(:),dmft_bath)
!   endif
!   call write_fit_result(ispin,iorb)
!   deallocate(fgand)
!   !
!   call get_dmft_bath(dmft_bath,bath_)
!   call deallocate_dmft_bath(dmft_bath)
!   deallocate(Gdelta,Xdelta,Wdelta)
!   !
! contains
!   !
!   subroutine write_fit_result(ispin,iorb)
!     integer           :: i,j,iorb,ispin
!     suffix="_orb"//reg(txtfy(iorb))//"_s"//reg(txtfy(ispin))//reg(ed_file_suffix)
!     unit=free_unit()
!     if(cg_scheme=='weiss')then
!        open(unit,file="fit_weiss"//reg(suffix)//".ed")
!     else
!        open(unit,file="fit_delta"//reg(suffix)//".ed")
!     endif
!     do i=1,Ldelta
!        write(unit,"(5F24.15)")Xdelta(i),&
!             dimag(fg(iorb,iorb,i)),dimag(fgand(ispin,ispin,iorb,iorb,i)),&
!             dreal(fg(iorb,iorb,i)),dreal(fgand(ispin,ispin,iorb,iorb,i))
!     enddo
!     close(unit)
!   end subroutine write_fit_result
! end subroutine chi2_fitgf_normal_normal_OneOrb







!##################################################################
! THESE PROCEDURES EVALUATES THE \chi^2 FUNCTIONS TO MINIMIZE. 
!##################################################################
!+-------------------------------------------------------------+
!PURPOSE: Evaluate the \chi^2 distance of \Delta_Anderson function.
!+-------------------------------------------------------------+
function chi2_delta_normal_normal(a) result(chi2)
  real(8),dimension(:)         ::  a
  real(8)                      ::  chi2
  complex(8),dimension(Ldelta) ::  Delta
  real(8),dimension(Ldelta)    ::  Ctmp
  !
  Delta = delta_normal_normal(a)
  !
  Ctmp = abs(Gdelta(1,:)-Delta(:))
  chi2=sum( Ctmp**cg_pow/Wdelta )
  chi2=chi2/Ldelta
  !
end function chi2_delta_normal_normal

!+-------------------------------------------------------------+
!PURPOSE: Evaluate the gradient \Grad\chi^2 of 
! \Delta_Anderson function.
!+-------------------------------------------------------------+
function grad_chi2_delta_normal_normal(a) result(dchi2)
  real(8),dimension(:)                 :: a
  real(8),dimension(size(a))           :: dchi2
  real(8),dimension(size(a))           :: df
  complex(8),dimension(Ldelta)         :: Delta
  complex(8),dimension(Ldelta)         :: Ftmp
  real(8),dimension(Ldelta)            :: Ctmp
  complex(8),dimension(Ldelta,size(a)) :: dDelta
  integer                              :: j
  !
  Delta   = delta_normal_normal(a)
  dDelta  = grad_delta_normal_normal(a)
  !
  Ftmp = Gdelta(1,:)-Delta(:)
  Ctmp = abs(Ftmp)**(cg_pow-2)
  do j=1,size(a)
     df(j)=sum( dreal(Ftmp)*dreal(dDelta(:,j))*Ctmp/Wdelta ) + &
          sum(  dimag(Ftmp)*dimag(dDelta(:,j))*Ctmp/Wdelta )
  enddo
  !
  dchi2 = -cg_pow*df/Ldelta
  !
end function grad_chi2_delta_normal_normal

!+-------------------------------------------------------------+
!PURPOSE: Evaluate the \chi^2 distance of G_0_Anderson function 
! The Gradient is not evaluated, so the minimization requires 
! a numerical estimate of the gradient. 
!+-------------------------------------------------------------+
function chi2_weiss_normal_normal(a) result(chi2)
  real(8),dimension(:)         ::  a
  complex(8),dimension(Ldelta) ::  g0and
  real(8),dimension(Ldelta)    ::  Ctmp
  real(8)                      ::  chi2,w
  !
  g0and  = g0and_normal_normal(a)
  !
  Ctmp = abs(Gdelta(1,:)-g0and(:))
  chi2 = sum( Ctmp**cg_pow/Wdelta )
  chi2 = chi2/Ldelta
  !
end function chi2_weiss_normal_normal

!+-------------------------------------------------------------+
!PURPOSE: Evaluate the gradient \Grad\chi^2 of 
! \Delta_Anderson function.
!+-------------------------------------------------------------+
function grad_chi2_weiss_normal_normal(a) result(dchi2)
  real(8),dimension(:)                 :: a
  real(8),dimension(size(a))           :: dchi2
  real(8),dimension(size(a))           :: df
  complex(8),dimension(Ldelta)         :: g0and,Ftmp
  real(8),dimension(Ldelta)            :: Ctmp
  complex(8),dimension(Ldelta,size(a)) :: dg0and
  integer                              :: j
  !
  g0and  = g0and_normal_normal(a)
  dg0and = grad_g0and_normal_normal(a)
  !
  Ftmp = Gdelta(1,:)-g0and(:)
  Ctmp = abs(Ftmp)**(cg_pow-2)
  do j=1,size(a)
     df(j)=sum( dreal(Ftmp)*dreal(dg0and(:,j))*Ctmp/Wdelta ) + &
          sum(  dimag(Ftmp)*dimag(dg0and(:,j))*Ctmp/Wdelta )
  enddo
  !
  dchi2 = -cg_pow*df/Ldelta
  !
end function grad_chi2_weiss_normal_normal





!##################################################################
! THESE PROCEDURES EVALUATES THE 
! - \delta
! - \grad \delta
! - g0
! FUNCTIONS. 
!##################################################################

function delta_normal_normal(a) result(Delta)
  real(8),dimension(:)         :: a
  complex(8),dimension(Ldelta) :: Delta
  integer                      :: i,io,stride
  real(8),dimension(Nbath)     :: eps,vps
  !
  !\Delta_{aa} = \sum_k [ V_{a}(k) * V_{a}(k)/(iw_n - E_{a}(k)) ]
  !
  stride = 0
  do i=1,Nbath
     io = stride + i
     eps(i) = a(io) 
  enddo
  stride = Nbath
  do i=1,Nbath
     io = stride + i
     vps(i) = a(io)
  enddo
  !
  do i=1,Ldelta
     Delta(i) = sum( vps(:)*vps(:)/(xi*Xdelta(i) - eps(:)) )
  enddo
  !
end function delta_normal_normal

function grad_delta_normal_normal(a) result(dDelta)
  real(8),dimension(:)                 :: a
  complex(8),dimension(Ldelta,size(a)) :: dDelta
  integer                              :: i,k,ik,io,stride
  real(8),dimension(Nbath)             :: eps,vps
  complex(8)                           :: iw
  !
  !
  !\grad_{E_{a}(k)} \Delta_{bb}^{rr} = [ V_{a}(k)*V_{a}(k) / ( iw_n - E_{a}(k) )**2 ]
  !
  !\grad_{V_{a}(k)} \Delta_{bb}^{rr} = [ 2*V_{a}(k) / ( iw_n - E_{a}(k) ) ]
  !
  stride = 0
  do i=1,Nbath
     io = stride + i
     eps(i) = a(io) 
  enddo
  stride = Nbath
  do i=1,Nbath
     io = stride + i
     vps(i) = a(io)
  enddo
  !
  stride = 0
  do k=1,Nbath
     ik = stride + k
     dDelta(:,ik) = vps(k)*vps(k)/(xi*Xdelta(:) - eps(k))**2
  enddo
  stride = Nbath
  do k=1,Nbath
     ik = stride + k
     dDelta(:,ik) = 2d0*vps(k)/(xi*Xdelta(:) - eps(k))
  enddo
  !
end function grad_delta_normal_normal


function g0and_normal_normal(a) result(G0and)
  real(8),dimension(:)         :: a
  complex(8),dimension(Ldelta) :: G0and,Delta
  integer                      :: i,io,iorb,ispin
  !
  iorb   = Orb_indx
  ispin  = Spin_indx
  !
  Delta(:) = delta_normal_normal(a)
  G0and(:) = xi*Xdelta(:) + xmu - impHloc(ispin,ispin,iorb,iorb) - Delta(:)
  G0and(:) = one/G0and(:)
  !
end function g0and_normal_normal

function grad_g0and_normal_normal(a) result(dG0and)
  real(8),dimension(:)                 :: a
  complex(8),dimension(Ldelta)         :: G0and,Delta
  complex(8),dimension(Ldelta,size(a)) :: dG0and,dDelta
  integer                              :: k,iorb,ispin
  !
  iorb   = Orb_indx
  ispin  = Spin_indx
  !
  Delta  = delta_normal_normal(a)
  dDelta = grad_delta_normal_normal(a)
  G0and  = xi*Xdelta + xmu - impHloc(ispin,ispin,iorb,iorb) - Delta
  do k=1,size(a)
     dG0and(:,k) = one/G0and/G0and*dDelta(:,k)
  enddo
  !
end function grad_g0and_normal_normal























