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
!PURPOSE  : Chi^2 interface for
!+-------------------------------------------------------------+
subroutine chi2_fitgf_replica(fg,bath_)
  complex(8),dimension(:,:,:,:,:)             :: fg ![Nspin][Nspin][Norb][Norb][Lmats]
  logical(8),dimension(Nspin,Nspin,Norb,Norb) :: Hmask
  real(8),dimension(:),intent(inout)          :: bath_
  real(8),dimension(:),allocatable            :: array_bath
  integer                                     :: i,j,iorb,jorb,ispin,jspin,io,jo,ibath
  integer                                     :: iter,stride,counter,Asize
  real(8)                                     :: chi
  logical                                     :: check
  type(effective_bath)                        :: dmft_bath
  character(len=256)                          :: suffix
  integer                                     :: unit
  complex(8),dimension(:,:,:,:,:),allocatable :: fgand ![Nspin][][Norb][][Ldelta]  
  !
  if(size(fg,1)/=Nspin)stop "chi2_fitgf_replica error: size[fg,1]!=Nspin"
  if(size(fg,2)/=Nspin)stop "chi2_fitgf_replica error: size[fg,2]!=Nspin"
  if(size(fg,3)/=Norb)stop "chi2_fitgf_replica error: size[fg,3]!=Norb"
  if(size(fg,4)/=Norb)stop "chi2_fitgf_replica error: size[fg,4]!=Norb"
  !
  check= check_bath_dimension(bath_,impHloc)
  if(.not.check)stop "chi2_fitgf_replica error: wrong bath dimensions"
  !
  call allocate_dmft_bath(dmft_bath)
  call set_dmft_bath(bath_,dmft_bath)
  allocate(array_bath(size(bath_)-Nbath))
  allocate(Nlambdas(Nbath))
  counter=0
  do ibath=1,Nbath
     counter=counter+1
     Nlambdas(ibath)=bath_(counter) !N_dec, one per bath copy
  enddo
  array_bath=bath_(Nbath+1:size(bath_))
  !
  Hmask=mask_hloc(impHloc,wdiag=.true.,uplo=.true.)
  totNso=count(Hmask)
  !
  allocate(getIspin(totNso),getJspin(totNso))
  allocate(getIorb(totNso) ,getJorb(totNso))
  counter=0
  do ispin=1,Nspin
     do jspin=1,Nspin
        do iorb=1,Norb
           do jorb=1,Norb
              if (Hmask(ispin,jspin,iorb,jorb))then
                 counter=counter+1
                 getIspin(counter) = ispin
                 getIorb(counter)  = iorb
                 getJspin(counter) = jspin
                 getJorb(counter)  = jorb
              endif
           enddo
        enddo
     enddo
  enddo
  !
  Ldelta = Lfit ; if(Ldelta>size(fg,5))Ldelta=size(fg,5)
  !
  allocate(Gdelta(totNso,Ldelta))
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
  write(LOGfile,*)"  fitted functions",totNso
  do i=1,totNso
     Gdelta(i,1:Ldelta) = fg(getIspin(i),getJspin(i),getIorb(i),getJorb(i),1:Ldelta)
  enddo
  !
  select case(cg_method)     !0=NR-CG[default]; 1=CG-MINIMIZE; 2=CG+
  case default
     if(cg_grad==0)then
        select case (cg_scheme)
        case ("weiss")
           call fmin_cg(array_bath,chi2_weiss_replica,grad_chi2_weiss_replica,iter,chi,&
                itmax=cg_niter,&
                ftol=cg_Ftol,  &
                istop=cg_stop, &
                iverbose=(ed_verbose>3))
        case ("delta")
           call fmin_cg(array_bath,chi2_delta_replica,grad_chi2_delta_replica,iter,chi,&
                itmax=cg_niter,&
                ftol=cg_Ftol,  &
                istop=cg_stop, &
                iverbose=(ed_verbose>3))
        case default
           stop "chi2_fitgf_replica error: cg_scheme != [weiss,delta]"
        end select
     else
        select case (cg_scheme)
        case ("weiss")
           call fmin_cg(array_bath,chi2_weiss_replica,iter,chi,&
                itmax=cg_niter,&
                ftol=cg_Ftol,  &
                istop=cg_stop, &
                iverbose=(ed_verbose>3))
        case ("delta")
           call fmin_cg(array_bath,chi2_delta_replica,iter,chi,&
                itmax=cg_niter,&
                ftol=cg_Ftol,  &
                istop=cg_stop, &
                iverbose=(ed_verbose>3))
        case default
           stop "chi2_fitgf_replica error: cg_scheme != [weiss,delta]"
        end select
     endif
     !
     !
  case (1)
     select case (cg_scheme)
     case ("weiss")
        call fmin_cgminimize(array_bath,chi2_weiss_replica,&
             iter,chi,itmax=cg_niter,ftol=cg_Ftol,new_version=cg_minimize_ver,hh_par=cg_minimize_hh,&
             iverbose=(ed_verbose>3))
     case ("delta")
        call fmin_cgminimize(array_bath,chi2_delta_replica,&
             iter,chi,itmax=cg_niter,ftol=cg_Ftol,new_version=cg_minimize_ver,hh_par=cg_minimize_hh,&
             iverbose=(ed_verbose>3))
     case default
        stop "chi2_fitgf_replica error: cg_scheme != [weiss,delta]"
     end select
     !
  end select
  !
  write(LOGfile,"(A,ES18.9,A,I5,A)")"chi^2|iter"//reg(ed_file_suffix)//'= ',chi," | ",iter,"  <--  All Orbs, All Spins"
  !
  suffix="_ALLorb_ALLspins"//reg(ed_file_suffix)
  unit=free_unit()
  open(unit,file="chi2fit_results"//reg(suffix)//".ed",position="append")
  write(unit,"(ES18.9,1x,I5)") chi,iter
  close(unit)
  !
  bath_(Nbath+1:size(bath_))=array_bath
  call set_dmft_bath(bath_,dmft_bath) ! *** array_bath --> dmft_bath *** (per write fit result)
  call write_dmft_bath(dmft_bath,LOGfile)
  call save_dmft_bath(dmft_bath)
  !
  allocate(fgand(Nspin,Nspin,Norb,Norb,Ldelta))
  if(cg_scheme=='weiss')then
     fgand = g0and_bath_function(xi*Xdelta(:),dmft_bath)
  else
     fgand = delta_bath_function(xi*Xdelta(:),dmft_bath)
  endif
  call write_fit_result()
  deallocate(fgand)
  !
  call get_dmft_bath(dmft_bath,bath_)                ! ***  dmft_bath --> bath_ ***    (bath in output)
  call deallocate_dmft_bath(dmft_bath)
  deallocate(Gdelta,Xdelta,Wdelta)
  deallocate(getIspin,getJspin)
  deallocate(getIorb,getJorb)
  deallocate(array_bath)
  deallocate(Nlambdas)
  !
contains
  !
  subroutine write_fit_result()
    integer   :: i,j,s,l,iorb,jorb,ispin,jspin
    !
    do l=1,totNso
       iorb = getIorb(l)
       jorb = getJorb(l)
       ispin = getIspin(l)
       jspin = getJspin(l)
       suffix="_l"//reg(txtfy(iorb))//&
            "_m"//reg(txtfy(jorb))//&
            "_s"//reg(txtfy(ispin))//&
            "_r"//reg(txtfy(jspin))//reg(ed_file_suffix)
       unit=free_unit()
       if(cg_scheme=='weiss')then
          open(unit,file="fit_weiss"//reg(suffix)//".ed")
       else
          open(unit,file="fit_delta"//reg(suffix)//".ed")
       endif
       do i=1,Ldelta
          write(unit,"(5F24.15)")Xdelta(i),&
               dimag(fg(ispin,jspin,iorb,jorb,i)),dimag(fgand(ispin,jspin,iorb,jorb,i)),&
               dreal(fg(ispin,jspin,iorb,jorb,i)),dreal(fgand(ispin,jspin,iorb,jorb,i))
       enddo
       close(unit)
    enddo
  end subroutine write_fit_result
  !
end subroutine chi2_fitgf_replica


!##################################################################
! THESE PROCEDURES EVALUATES THE \chi^2 FUNCTIONS TO MINIMIZE. 
!##################################################################
!+-------------------------------------------------------------+
!PURPOSE: Evaluate the \chi^2 distance of \Delta_Anderson function.
!+-------------------------------------------------------------+
function chi2_delta_replica(a) result(chi2)
  real(8),dimension(:)                               :: a
  real(8)                                            :: chi2
  real(8),dimension(totNso)                          :: chi2_so
  complex(8),dimension(Nspin,Nspin,Norb,Norb,Ldelta) :: Delta
  real(8),dimension(Ldelta)                          :: Ctmp
  integer                                            :: i,l,iorb,jorb,ispin,jspin
  !
  Delta = delta_replica(a)
  !
  do l=1,totNso
     iorb = getIorb(l)
     jorb = getJorb(l)
     ispin = getIspin(l)
     jspin = getJspin(l)
     !
     Ctmp =  abs(Gdelta(l,:)-Delta(ispin,jspin,iorb,jorb,:))
     chi2_so(l) = sum( Ctmp**cg_pow/Wdelta )
  enddo
  !
  chi2=sum(chi2_so)
  chi2=chi2/Ldelta
  !
end function chi2_delta_replica


!+-------------------------------------------------------------+
!PURPOSE: Evaluate the gradient \Grad\chi^2 of 
! \Delta_Anderson function.
!+-------------------------------------------------------------+
function grad_chi2_delta_replica(a) result(dchi2)
  real(8),dimension(:)                                       :: a
  real(8),dimension(size(a))                                 :: dchi2
  real(8),dimension(totNso,size(a))                          :: df
  complex(8),dimension(Nspin,Nspin,Norb,Norb,Ldelta)         :: Delta
  complex(8),dimension(Nspin,Nspin,Norb,Norb,Ldelta,size(a)) :: dDelta
  complex(8),dimension(Ldelta)                               :: Ftmp
  real(8),dimension(Ldelta)                                  :: Ctmp
  integer                                                    :: i,j,l,iorb,jorb,ispin,jspin
  !
  Delta  = delta_replica(a)
  dDelta = grad_delta_replica(a)
  !
  do l=1,totNso
     iorb = getIorb(l)
     jorb = getJorb(l)
     ispin = getIspin(l)
     jspin = getJspin(l)
     !
     Ftmp = Gdelta(l,:)-Delta(ispin,jspin,iorb,jorb,:)
     Ctmp = abs(Ftmp)**(cg_pow-2)
     do j=1,size(a)
        df(l,j)=&
             sum( dreal(Ftmp)*dreal(dDelta(ispin,jspin,iorb,jorb,:,j))*Ctmp/Wdelta ) + &
             sum( dimag(Ftmp)*dimag(dDelta(ispin,jspin,iorb,jorb,:,j))*Ctmp/Wdelta )
     enddo
  enddo
  !
  dchi2 = -cg_pow*sum(df,1)/Ldelta     !sum over all orbital indices
  !
end function grad_chi2_delta_replica



!+-------------------------------------------------------------+
!PURPOSE: Evaluate the \chi^2 distance of G_0_Anderson function 
! The Gradient is not evaluated, so the minimization requires 
! a numerical estimate of the gradient. 
!+-------------------------------------------------------------+
function chi2_weiss_replica(a) result(chi2)
  real(8),dimension(:)                               :: a
  real(8)                                            :: chi2
  real(8),dimension(totNso)                          :: chi2_so
  complex(8),dimension(Nspin,Nspin,Norb,Norb,Ldelta) :: g0and
  real(8),dimension(Ldelta)                          :: Ctmp
  integer                                            :: i,l,iorb,jorb,ispin,jspin
  !
  g0and = g0and_replica(a)
  !
  do l=1,totNso
     iorb = getIorb(l)
     jorb = getJorb(l)
     ispin = getIspin(l)
     jspin = getJspin(l)
     !
     Ctmp = abs(Gdelta(l,:)-g0and(ispin,jspin,iorb,jorb,:))
     chi2_so(l) = sum( Ctmp**cg_pow/Wdelta )
  enddo
  !
  chi2=sum(chi2_so)
  chi2=chi2/Ldelta
  !
end function chi2_weiss_replica

!+-------------------------------------------------------------+
!PURPOSE: Evaluate the gradient \Grad\chi^2 of 
! \Delta_Anderson function.
!+-------------------------------------------------------------+
function grad_chi2_weiss_replica(a) result(dchi2)
  real(8),dimension(:)                                       :: a
  real(8),dimension(size(a))                                 :: dchi2
  real(8),dimension(totNso,size(a))                          :: df
  complex(8),dimension(Nspin,Nspin,Norb,Norb,Ldelta)         :: g0and
  complex(8),dimension(Nspin,Nspin,Norb,Norb,Ldelta,size(a)) :: dg0and
  complex(8),dimension(Ldelta)                               :: Ftmp
  real(8),dimension(Ldelta)                                  :: Ctmp
  integer                                                    :: i,j,l,iorb,jorb,ispin,jspin
  !
  g0and  = g0and_replica(a)
  dg0and = grad_g0and_replica(a)
  !
  do l=1,totNso
     iorb = getIorb(l)
     jorb = getJorb(l)
     ispin = getIspin(l)
     jspin = getJspin(l)
     !
     Ftmp = Gdelta(l,:)-g0and(ispin,jspin,iorb,jorb,:)
     Ctmp = abs(Ftmp)**(cg_pow-2)
     do j=1,size(a)
        df(l,j)=&
             sum( dreal(Ftmp)*dreal(dg0and(ispin,jspin,iorb,jorb,:,j))*Ctmp/Wdelta ) + &
             sum( dimag(Ftmp)*dimag(dg0and(ispin,jspin,iorb,jorb,:,j))*Ctmp/Wdelta )
     enddo
  enddo
  !
  dchi2 = -cg_pow*sum(df,1)/Ldelta     !sum over all orbital indices
  !
end function grad_chi2_weiss_replica





!##################################################################
! THESE PROCEDURES EVALUATES THE 
! - \delta
! - g0
! FUNCTIONS. 
!##################################################################
!ACHTUNG! We use a direct dump of the array into the necessary element of the bath.
! rather than using aux functions in ED_BATH. This improves execution speed. 
function delta_replica(a) result(Delta)
  real(8),dimension(:)                               :: a
  complex(8),dimension(Nspin,Nspin,Norb,Norb,Ldelta) :: Delta
  integer                                            :: ispin,jspin,iorb,jorb,ibath
  integer                                            :: i,stride
  complex(8),dimension(Nspin*Norb,Nspin*Norb)        :: Haux,Htmp
  complex(8),dimension(Nspin,Nspin,Norb,Norb)        :: invH_knn
  real(8),dimension(Nspin,Nbath)                     :: dummy_Vbath
  type(nsymm_vector),dimension(Nbath)                :: dummy_lambda
  !
  !Get Hs
  stride = 0
  do ibath=1,Nbath
     allocate(dummy_lambda(ibath)%element(Nlambdas(ibath)))
     !
     do ispin=1,Nspin
        stride = stride + 1
        dummy_vbath(ispin,ibath) = a(stride)
     enddo
     dummy_lambda(ibath)%element=a(stride+1:stride+Nlambdas(ibath))
     stride=stride+Nlambdas(ibath)
  enddo
  !
  Delta=zero
  do ibath=1,Nbath
     invH_knn = bath_from_sym(dummy_lambda(ibath)%element)
     Htmp     = nn2so_reshape( invH_knn,Nspin,Norb)
     do i=1,Ldelta
        Haux     = zeye(Nspin*Norb)*xi*Xdelta(i) - Htmp
        call inv(Haux)
        invH_knn = so2nn_reshape(Haux,Nspin,Norb)
        !
        do ispin=1,Nspin
           Delta(ispin,ispin,:,:,i)=Delta(ispin,ispin,:,:,i) + &
                dummy_Vbath(ispin,ibath)*invH_knn(ispin,ispin,:,:)*dummy_Vbath(ispin,ibath)
        enddo
     enddo
     deallocate(dummy_lambda(ibath)%element)
  enddo
  !
end function delta_replica

function g0and_replica(a) result(G0and)
  real(8),dimension(:)                               :: a
  complex(8),dimension(Nspin,Nspin,Norb,Norb,Ldelta) :: G0and,Delta
  complex(8),dimension(Norb*Nspin,Norb*Nspin)        :: FGorb
  integer                                            :: i,Nso
  !
  Nso   = Norb*Nspin
  Delta = delta_replica(a)
  do i=1,Ldelta
     FGorb = (xi*Xdelta(i)+xmu)*zeye(Nso) - nn2so_reshape(impHloc + Delta(:,:,:,:,i), Nspin,Norb)
     call inv(FGorb)
     G0and(:,:,:,:,i) = so2nn_reshape(FGorb,Nspin,Norb)
  enddo
end function g0and_replica



function grad_delta_replica(a) result(dDelta)
  real(8),dimension(:)                                      :: a
  complex(8),dimension(Nspin,Nspin,Norb,Norb,Ldelta,size(a)) :: dDelta
  integer                                                   :: ispin,iorb,jorb,ibath
  integer                                                   :: i,k,ik,l,io,counter
  complex(8),dimension(Nspin*Norb,Nspin*Norb)               :: H_reconstructed,Htmp,Hbasis_so
  complex(8),dimension(Nspin*Norb,Nspin*Norb,Ldelta)        :: Haux
  complex(8),dimension(Nspin,Nspin,Norb,Norb,Ldelta)        :: invH_knn
  real(8),dimension(Nspin,Nbath)                            :: dummy_Vbath
  type(nsymm_vector),dimension(Nbath)                       :: dummy_lambda
  !
  !
  !Get Hs
  counter = 0
  do ibath=1,Nbath
     allocate(dummy_lambda(ibath)%element(Nlambdas(ibath)))
     !
     do ispin=1,Nspin
        counter = counter + 1
        dummy_vbath(ispin,ibath) = a(counter)
     enddo
     dummy_lambda(ibath)%element=a(counter+1:counter+Nlambdas(ibath))
     counter=counter+Nlambdas(ibath)
  enddo
  !
  dDelta=zero
  counter=0
  do ibath=1,Nbath
     H_reconstructed     = nn2so_reshape( bath_from_sym(dummy_lambda(ibath)%element) ,Nspin,Norb)
     do i=1,Ldelta
        Haux(:,:,i) = zeye(Nspin*Norb)*xi*Xdelta(i) - H_reconstructed
        call inv(Haux(:,:,i))
        invH_knn(:,:,:,:,i) = so2nn_reshape(Haux(:,:,i),Nspin,Norb)
     enddo
     !Derivate_Vp
     do ispin=1,Nspin
        counter = counter + 1
        do iorb=1,Norb
           do jorb=1,Norb
              dDelta(ispin,ispin,iorb,jorb,:,counter)=2d0*dummy_Vbath(ispin,ibath)*invH_knn(ispin,ispin,iorb,jorb,:)
           enddo
        enddo
     enddo
     !Derivate_lambda_p
     do ispin=1,Nspin
        do k=1,size(dummy_lambda(ibath)%element)
           counter = counter + 1
           Hbasis_so=nn2so_reshape(H_basis(k)%O,Nspin,Norb)
           do l=1,Ldelta
              ! Hbasis_lso=nnn2lso_reshape(H_basis(k)%O,Nlat,Nspin,Norb)
              ! Htmp=matmul(Haux(:,:,l),Hbasis_lso)
              ! Htmp=matmul(Htmp,Haux(:,:,l))
              Htmp = ((Haux(:,:,l) .x. Hbasis_so)) .x. Haux(:,:,l)
              dDelta(:,:,:,:,l,counter)=so2nn_reshape((dummy_Vbath(ispin,ibath)**2)*Htmp,Nspin,Norb)
           enddo
        enddo
     enddo
     !
  enddo
end function grad_delta_replica


function grad_g0and_replica(a) result(dG0and)
  real(8),dimension(:)                                       :: a
  complex(8),dimension(Nspin,Nspin,Norb,Norb,Ldelta,size(a)) :: dG0and
  complex(8),dimension(Nspin,Nspin,Norb,Norb,Ldelta)         :: G0and
  complex(8),dimension(Nspin,Nspin,Norb,Norb,Ldelta,size(a)) :: dDelta
  complex(8),dimension(Nspin*Norb,Nspin*Norb)                :: dDelta_so,dG0and_so,G0and_so
  integer                                                    :: ispin,iorb,jorb
  integer                                                    :: ik,l
  !
  G0and  = g0and_replica(a)
  dDelta = grad_delta_replica(a)
  !
  dG0and = zero
  do l=1,Ldelta
     G0and_so=nn2so_reshape(g0and(:,:,:,:,l),Nspin,Norb)
     do ik=1,size(a)
        dDelta_so=nn2so_reshape(dDelta(:,:,:,:,l,ik),Nspin,Norb)
        ! dG0and_lso=matmul(-G0and_lso,dDelta_lso)
        ! dG0and_lso=matmul(dG0and_lso,G0and_lso)
        dG0and_so = (G0and_so .x. dDelta_so) .x. G0and_so
        dG0and(:,:,:,:,l,ik)=-so2nn_reshape(dG0and_so,Nspin,Norb)
     enddo
  enddo
  !
end function grad_g0and_replica











! function delta_replica(a) result(Delta)
!   real(8),dimension(:)                                :: a
!   complex(8),dimension(Nspin,Nspin,Norb,Norb,Ldelta)  :: Delta
!   integer                                             :: ispin,jspin,iorb,jorb,ibath
!   integer                                             :: i,io,jo,ndx
!   complex(8),dimension(Nspin*Norb,Nspin*Norb)         :: V_k
!   complex(8),dimension(Nspin*Norb,Nspin*Norb,Ldelta)  :: invH_k
!   complex(8),dimension(Nspin,Nspin,Norb,Norb,Nbath)   :: invH_knn
!   type(effective_bath)                                :: dmft_bath_tmp
!   !
!   call allocate_dmft_bath(dmft_bath_tmp)
!   call init_dmft_bath_mask(dmft_bath_tmp)
!   call set_dmft_bath(a,dmft_bath_tmp)
!   !
!   Delta=zero
!   invH_k=zero
!   do i=1,Ldelta
!      invH_knn=zero
!      do ibath=1,Nbath
!         !
!         do ispin=1,Nspin
!            do jspin=1,Nspin
!               do iorb=1,Norb
!                  do jorb=1,Norb
!                     io = iorb + (ispin-1) * Norb
!                     jo = jorb + (jspin-1) * Norb
!                     invH_k(io,jo,i)=dmft_bath_tmp%h(ispin,jspin,iorb,jorb,ibath)
!                  enddo
!               enddo
!            enddo
!         enddo
!         !
!         invH_k(:,:,i) = zeye(Nspin*Norb) * xi * Xdelta(i) - invH_k(:,:,i)
!         call inv(invH_k(:,:,i))
!         invH_knn(:,:,:,:,ibath)=so2nn_reshape(invH_k(:,:,i),Nspin,Norb)
!         !
!      enddo
!      !
!      do ibath=1,Nbath
!         do ispin=1,Nspin
!            do jspin=1,Nspin
!               do iorb=1,Norb
!                  do jorb=1,Norb
!                     Delta(ispin,jspin,iorb,jorb,i)=Delta(ispin,jspin,iorb,jorb,i)+ &
!                          (dmft_bath_tmp%vr(ibath))*invH_knn(ispin,jspin,iorb,jorb,ibath)*dmft_bath_tmp%vr(ibath)
!                  enddo
!               enddo
!            enddo
!         enddo
!      enddo
!   enddo
!   !
!   call deallocate_dmft_bath(dmft_bath_tmp)
!   !
! end function delta_replica

! function g0and_replica(a) result(G0and)
!   real(8),dimension(:)                               :: a
!   complex(8),dimension(Nspin,Nspin,Norb,Norb,Ldelta) :: G0and,Delta
!   complex(8),dimension(Norb*Nspin,Norb*Nspin)        :: zeta,fgorb
!   integer                                            :: i,Nso
!   integer                                            :: iorb,jorb,ispin,jspin,io,jo
!   !
!   Nso = Norb*Nspin
!   !
!   Delta = delta_replica(a)
!   !
!   do i=1,Ldelta
!      fgorb=zero
!      zeta = (xi*Xdelta(i)+xmu)*zeye(Nso)
!      do ispin=1,Nspin
!         do jspin=1,Nspin
!            do iorb=1,Norb
!               do jorb=1,Norb
!                  io = iorb + (ispin-1)*Norb
!                  jo = jorb + (jspin-1)*Norb
!                  fgorb(io,jo)   = zeta(io,jo) - impHloc(ispin,jspin,iorb,jorb) - Delta(ispin,jspin,iorb,jorb,i)
!               enddo
!            enddo
!         enddo
!      enddo
!      call inv(fgorb)
!      do ispin=1,Nspin
!         do jspin=1,Nspin
!            do iorb=1,Norb
!               do jorb=1,Norb
!                  io = iorb + (ispin-1)*Norb
!                  jo = jorb + (jspin-1)*Norb
!                  G0and(ispin,jspin,iorb,jorb,i) = fgorb(io,jo)
!               enddo
!            enddo
!         enddo
!      enddo
!   enddo
!   !
! end function g0and_replica


