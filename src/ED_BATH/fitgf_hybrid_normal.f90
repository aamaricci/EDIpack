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
!PURPOSE  : Chi^2 interface for Hybrid bath normal phase
!+-------------------------------------------------------------+
subroutine chi2_fitgf_hybrid_normal(fg,bath_,ispin)
  complex(8),dimension(:,:,:)                 :: fg ![Norb][Norb][Lmats]
  real(8),dimension(:),intent(inout)          :: bath_
  integer                                     :: ispin
  real(8),dimension(:),allocatable            :: array_bath
  integer                                     :: iter,stride,i,io,j,corb,l,Asize
  integer                                     :: iorb,jorb
  real(8)                                     :: chi
  logical                                     :: check
  type(effective_bath)                        :: dmft_bath
  character(len=256)                          :: suffix
  integer                                     :: unit
  complex(8),dimension(:,:,:,:,:),allocatable :: fgand ![Nspin][][Norb][][Ldelta]
  !
  if(size(fg,1)/=Norb)stop "chi2_fitgf_hybrid_normal error: size[fg,1]!=Norb"
  if(size(fg,2)/=Norb)stop "chi2_fitgf_hybrid_normal error: size[fg,2]!=Norb"
  !
  check= check_bath_dimension(bath_)
  if(.not.check)stop "chi2_fitgf_hybrid_normal error: wrong bath dimensions"
  !
  Ldelta = Lfit ; if(Ldelta>size(fg,3))Ldelta=size(fg,3)
  !
  allocate(getIorb(Norb*(Norb+1)/2),getJorb(Norb*(Norb+1)/2))
  corb=0
  do iorb=1,Norb
     do jorb=iorb,Norb
        corb=corb+1
        getIorb(corb)=iorb
        getJorb(corb)=jorb
     enddo
  enddo
  totNorb=corb
  if(totNorb/=(Norb*(Norb+1)/2))stop "chi2_fitgf_hybrid_normal error counting the orbitals"
  !
  allocate(Gdelta(totNorb,Ldelta))
  allocate(Xdelta(Ldelta))
  allocate(Wdelta(Ldelta))
  !
  do i=1,totNorb
     Gdelta(i,1:Ldelta) = fg(getIorb(i),getJorb(i),1:Ldelta)
  enddo
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
  Spin_indx=ispin
  !
  call allocate_dmft_bath(dmft_bath)
  call set_dmft_bath(bath_,dmft_bath)
  !
  !E_{\s,1}(:)  [ 1 ][ 1 ][Nbath]
  !V_{\s,:}(:)  [ 1 ][ Norb][Nbath]
  Asize = Nbath + Norb*Nbath
  allocate(array_bath(Asize))
  !
  !Nbath + Norb*Nbath
  stride = 0
  do i=1,Nbath
     io = stride + i
     array_bath(io) = dmft_bath%e(ispin,1,i)
  enddo
  stride = Nbath
  do iorb=1,Norb
     do i=1,Nbath
        io = stride + i + (iorb-1)*Nbath
        array_bath(io) = dmft_bath%v(ispin,iorb,i)
     enddo
  enddo
  !
  select case(cg_method)     !0=NR-CG[default]; 1=CG-MINIMIZE; 2=CG+
  case default
     if(cg_grad==0)then
        select case (cg_scheme)
        case ("weiss")
           call fmin_cg(array_bath,chi2_weiss_hybrid_normal,grad_chi2_weiss_hybrid_normal,iter,chi,&
                itmax=cg_niter,&
                ftol=cg_Ftol,  &
                istop=cg_stop, &
                iverbose=(ed_verbose>3))
        case ("delta")
           call fmin_cg(array_bath,chi2_delta_hybrid_normal,grad_chi2_delta_hybrid_normal,iter,chi,&
                itmax=cg_niter,&
                ftol=cg_Ftol,  &
                istop=cg_stop, &
                iverbose=(ed_verbose>3))
        case default
           stop "chi2_fitgf_hybrid_normal error: cg_scheme != [weiss,delta]"
        end select
     else
        select case (cg_scheme)
        case ("weiss")
           call fmin_cg(array_bath,chi2_weiss_hybrid_normal,iter,chi,&
                itmax=cg_niter,&
                ftol=cg_Ftol, &
                istop=cg_stop,&
                iverbose=(ed_verbose>3))
        case ("delta")
           call fmin_cg(array_bath,chi2_delta_hybrid_normal,iter,chi,&
                itmax=cg_niter,&
                ftol=cg_Ftol, &
                istop=cg_stop,&
                iverbose=(ed_verbose>3))
        case default
           stop "chi2_fitgf_hybrid_normal error: cg_scheme != [weiss,delta]"
        end select
     endif
     !
     !
  case (1)
     select case (cg_scheme)
     case ("weiss")
        call fmin_cgminimize(array_bath,chi2_weiss_hybrid_normal,&
             iter,chi,itmax=cg_niter,ftol=cg_Ftol,&
             new_version=cg_minimize_ver,&
             hh_par=cg_minimize_hh,&
             iverbose=(ed_verbose>3))
     case ("delta")
        call fmin_cgminimize(array_bath,chi2_delta_hybrid_normal,&
             iter,chi,itmax=cg_niter,ftol=cg_Ftol,&
             new_version=cg_minimize_ver,&
             hh_par=cg_minimize_hh,&
             iverbose=(ed_verbose>3))
     case default
        stop "chi2_fitgf_hybrid_normal error: cg_scheme != [weiss,delta]"
     end select
     !
  end select
  !
  write(LOGfile,"(A,ES18.9,A,I5)")&
       'chi^2|iter'//reg(ed_file_suffix)//'=',chi," | ",iter,&
       "  <--  all Orbs, Spin"//reg(txtfy(ispin))
  !
  suffix="_ALLorb_s"//reg(txtfy(ispin))//reg(ed_file_suffix)
  unit=free_unit()
  open(unit,file="chi2fit_results"//reg(suffix)//".ed",position="append")
  write(unit,"(ES18.9,1x,I5)") chi,iter
  close(unit)
  !
  stride = 0
  do i=1,Nbath
     io = stride + i
     dmft_bath%e(ispin,1,i) = array_bath(io)
  enddo
  stride = Nbath
  do iorb=1,Norb
     do i=1,Nbath
        io = stride + i + (iorb-1)*Nbath
        dmft_bath%v(ispin,iorb,i) = array_bath(io)
     enddo
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
  deallocate(getIorb,getJorb)
  !
contains
  !
  subroutine write_fit_result(ispin)
    integer :: i,j,l,m,iorb,jorb,ispin,jspin
    !
    do l=1,totNorb
       iorb=getIorb(l)
       jorb=getJorb(l)
       suffix="_l"//reg(txtfy(iorb))//"_m"//reg(txtfy(jorb))//reg(ed_file_suffix)
       unit=free_unit()
       if(cg_scheme=='weiss')then
          open(unit,file="fit_weiss"//reg(suffix)//".ed")
       else
          open(unit,file="fit_delta"//reg(suffix)//".ed")
       endif
       do i=1,Ldelta
          write(unit,"(5F24.15)")Xdelta(i),&
               dimag(Gdelta(l,i)),dimag(fgand(ispin,ispin,iorb,jorb,i)),&
               dreal(Gdelta(l,i)),dreal(fgand(ispin,ispin,iorb,jorb,i))
       enddo
       close(unit)
    enddo
  end subroutine write_fit_result
end subroutine chi2_fitgf_hybrid_normal








!##################################################################
! THESE PROCEDURES EVALUATES THE \chi^2 FUNCTIONS TO MINIMIZE. 
!##################################################################
!+-------------------------------------------------------------+
!PURPOSE: Evaluate the \chi^2 distance of \Delta_Anderson function.
!+-------------------------------------------------------------+
function chi2_delta_hybrid_normal(a) result(chi2)
  real(8),dimension(:)                   :: a
  real(8)                                :: chi2
  real(8),dimension(totNorb)             :: chi2_orb
  complex(8),dimension(Norb,Norb,Ldelta) :: Delta
  real(8),dimension(Ldelta)              :: Ctmp
  integer                                :: i,l,iorb,jorb
  !
  Delta = delta_hybrid_normal(a)
  !
  do l=1,totNorb
     iorb=getIorb(l)
     jorb=getJorb(l)
     Ctmp=abs(Gdelta(l,:)-Delta(iorb,jorb,:))
     chi2_orb(l) = sum( Ctmp**cg_pow/Wdelta )
  enddo
  !
  chi2=sum(chi2_orb)
  chi2=chi2/Ldelta
  !
end function chi2_delta_hybrid_normal

!+-------------------------------------------------------------+
!PURPOSE: Evaluate the gradient \Grad\chi^2 of 
! \Delta_Anderson function.
!+-------------------------------------------------------------+
function grad_chi2_delta_hybrid_normal(a) result(dchi2)
  real(8),dimension(:)                           :: a
  real(8),dimension(size(a))                     :: dchi2
  real(8),dimension(totNorb,size(a))             :: df
  complex(8),dimension(Norb,Norb,Ldelta)         :: Delta
  complex(8),dimension(Norb,Norb,Ldelta,size(a)) :: dDelta
  complex(8),dimension(Ldelta)                   :: Ftmp
  real(8),dimension(Ldelta)                      :: Ctmp
  integer                                        :: i,j,l,iorb,jorb
  !
  Delta  = delta_hybrid_normal(a)
  dDelta = grad_delta_hybrid_normal(a)
  !
  do l=1,totNorb
     iorb=getIorb(l)
     jorb=getJorb(l)
     !
     Ftmp = Gdelta(l,:)-Delta(iorb,jorb,:)
     Ctmp = abs(Ftmp)**(cg_pow-2)
     do j=1,size(a)
        df(l,j)=&
             sum( dreal(Ftmp)*dreal(dDelta(iorb,jorb,:,j))*Ctmp/Wdelta ) + &
             sum( dimag(Ftmp)*dimag(dDelta(iorb,jorb,:,j))*Ctmp/Wdelta )
     enddo
  enddo
  !
  dchi2 = -cg_pow*sum(df,1)/Ldelta     !sum over all orbital indices
  !
end function grad_chi2_delta_hybrid_normal




!+-------------------------------------------------------------+
!PURPOSE: Evaluate the \chi^2 distance of G_0_Anderson function 
! The Gradient is not evaluated, so the minimization requires 
! a numerical estimate of the gradient. 
!+-------------------------------------------------------------+
function chi2_weiss_hybrid_normal(a) result(chi2)
  real(8),dimension(:)                   :: a
  real(8),dimension(totNorb)             :: chi2_orb
  complex(8),dimension(Norb,Norb,Ldelta) :: g0and
  real(8),dimension(Ldelta)              :: Ctmp
  real(8)                                :: chi2
  integer                                :: i,l,iorb,jorb
  !
  g0and(:,:,:) = g0and_hybrid_normal(a)
  !

  do l=1,totNorb
     iorb = getIorb(l)
     jorb = getJorb(l)
     Ctmp = abs(Gdelta(l,:)-g0and(iorb,jorb,:))
     chi2_orb(l) = sum( Ctmp**cg_pow/Wdelta )
  enddo
  !
  chi2=sum(chi2_orb)/Norb
  chi2=chi2/Ldelta
  !
end function chi2_weiss_hybrid_normal

!+-------------------------------------------------------------+
!PURPOSE: Evaluate the gradient \Grad\chi^2 of 
! \Delta_Anderson function.
!+-------------------------------------------------------------+
function grad_chi2_weiss_hybrid_normal(a) result(dchi2)
  real(8),dimension(:)                           :: a
  real(8),dimension(size(a))                     :: dchi2
  real(8),dimension(totNorb,size(a))             :: df
  complex(8),dimension(Norb,Norb,Ldelta)         :: g0and
  complex(8),dimension(Norb,Norb,Ldelta,size(a)) :: dg0and
  complex(8),dimension(Ldelta)                   :: Ftmp
  real(8),dimension(Ldelta)                      :: Ctmp
  integer                                        :: i,j,l,iorb,jorb
  !
  g0and  = g0and_hybrid_normal(a)
  dg0and = grad_g0and_hybrid_normal(a)
  !
  do l=1,totNorb
     iorb=getIorb(l)
     jorb=getJorb(l)
     !
     Ftmp = Gdelta(l,:)-g0and(iorb,jorb,:)
     Ctmp = abs(Ftmp)**(cg_pow-2)
     do j=1,size(a)
        df(l,j)=&
             sum( dreal(Ftmp)*dreal(dg0and(iorb,jorb,:,j))*Ctmp/Wdelta ) + &
             sum( dimag(Ftmp)*dimag(dg0and(iorb,jorb,:,j))*Ctmp/Wdelta )
     enddo
  enddo
  !
  dchi2 = -cg_pow*sum(df,1)/Ldelta     !sum over all orbital indices
  !
end function grad_chi2_weiss_hybrid_normal





!##################################################################
! THESE PROCEDURES EVALUATES THE 
! - \delta
! - \grad \delta
! - g0
! FUNCTIONS. 
!##################################################################
function delta_hybrid_normal(a) result(Delta)
  real(8),dimension(:)                   :: a
  complex(8),dimension(Norb,Norb,Ldelta) :: Delta
  integer                                :: iorb,jorb
  integer                                :: i,io,l,stride
  real(8),dimension(Nbath)               :: eps
  real(8),dimension(Norb,Nbath)          :: vops
  !
  !\Delta_{ab} = \sum_k [ V_{a}(k) * V_{b}(k)/(iw_n - E(k)) ]
  !
  stride = 0
  do i=1,Nbath
     io = stride + i
     eps(i)    = a(io)
  enddo
  stride = Nbath
  do l=1,Norb
     do i=1,Nbath
        io = stride + i + (l-1)*Nbath
        vops(l,i) = a(io)
     enddo
  enddo
  !
  do i=1,Ldelta
     do iorb=1,Norb
        Delta(iorb,iorb,i) = sum( vops(iorb,:)*vops(iorb,:)/(xi*Xdelta(i) - eps(:)) )
        do jorb=iorb+1,Norb
           Delta(iorb,jorb,i) = sum( vops(iorb,:)*vops(jorb,:)/(xi*Xdelta(i) - eps(:)) )
           Delta(jorb,iorb,i) = sum( vops(jorb,:)*vops(iorb,:)/(xi*Xdelta(i) - eps(:)) )
        enddo
     enddo
  enddo
  !
end function delta_hybrid_normal

function grad_delta_hybrid_normal(a) result(dDelta)
  real(8),dimension(:)                           :: a
  complex(8),dimension(Norb,Norb,Ldelta,size(a)) :: dDelta
  integer                                        :: iorb,jorb
  integer                                        :: i,k,ik,l,io,stride
  real(8),dimension(Nbath)                       :: eps
  real(8),dimension(Norb,Nbath)                  :: vops
  real(8),dimension(Norb,Norb)                   :: delta_orb
  !
  !\grad_{E_{1}(k)} \Delta_{ab} = [ V_{a}(k)*V_{b}(k) / ( iw_n - E_{1}(k) )**2 ]
  !
  !\grad_{V_{g}(k)} \Delta_{ab} = [ \d(g,a)*V_{b}(k)+\d(g,b)*V_{a}(k) ] / ( iw_n - E_{1}(k) )
  !
  stride = 0
  do i=1,Nbath
     io = stride + i
     eps(i)    = a(io)
  enddo
  stride = Nbath
  do l=1,Norb
     do i=1,Nbath
        io = stride + i + (l-1)*Nbath
        vops(l,i) = a(io)
     enddo
  enddo
  !
  delta_orb = eye(Norb)
  !
  do iorb=1,Norb
     do jorb=1,Norb
        stride = 0
        do k=1,Nbath
           ik = stride + k
           dDelta(iorb,jorb,:,ik) = vops(iorb,k)*vops(jorb,k)/(xi*Xdelta(:) - eps(k))**2
        enddo
        stride = Nbath
        do l=1,Norb
           do k=1,Nbath
              ik = stride + k + (l-1)*Nbath
              dDelta(iorb,jorb,:,ik) = (delta_orb(l,iorb)*vops(jorb,k) + delta_orb(l,jorb)*vops(iorb,k))/(xi*Xdelta(:) - eps(k))
           enddo
        enddo
     enddo
  enddo
  !
end function grad_delta_hybrid_normal




function g0and_hybrid_normal(a) result(G0and)
  real(8),dimension(:)                   :: a
  complex(8),dimension(Norb,Norb,Ldelta) :: G0and,Delta
  complex(8),dimension(Norb,Norb)        :: zeta,fgorb
  integer                                :: i,ispin
  !
  ispin  = Spin_indx
  !
  Delta = delta_hybrid_normal(a)
  !
  do i=1,Ldelta
     fgorb=zero
     zeta = (xi*Xdelta(i)+xmu)*eye(Norb)
     fgorb(:,:)   = zeta(:,:) - impHloc(ispin,ispin,:,:) - Delta(:,:,i)
     call inv(fgorb)
     G0and(:,:,i) = fgorb
  enddo
  !
end function g0and_hybrid_normal

function grad_g0and_hybrid_normal(a) result(dG0and)
  real(8),dimension(:)                           :: a
  complex(8),dimension(Norb,Norb,Ldelta,size(a)) :: dG0and
  complex(8),dimension(Norb,Norb,Ldelta)         :: G0and
  complex(8),dimension(Norb,Norb,Ldelta,size(a)) :: dDelta
  integer                                        :: iorb,jorb
  integer                                        :: ik
  !
  G0and  = g0and_hybrid_normal(a)
  dDelta = grad_delta_hybrid_normal(a)
  do iorb=1,Norb
     do jorb=1,Norb
        do ik=1,size(a)
           dG0and(iorb,jorb,:,ik) = G0and(iorb,jorb,:)*G0and(iorb,jorb,:)*dDelta(iorb,jorb,:,ik)
        enddo
     enddo
  enddo
  !
end function grad_g0and_hybrid_normal
