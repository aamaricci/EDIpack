MODULE ED_CHI_EXCT
  USE SF_CONSTANTS, only:one,xi,zero,pi
  USE SF_TIMER  
  USE SF_IOTOOLS, only: str,free_unit,reg,free_units,txtfy
  USE SF_LINALG,  only: inv,eigh,eye
  USE SF_SP_LINALG, only: sp_lanc_tridiag
  USE ED_INPUT_VARS
  USE ED_VARS_GLOBAL
  USE ED_IO                     !< this contains the routine to print GF,Sigma and G0
  USE ED_EIGENSPACE
  USE ED_BATH
  USE ED_BATH_FUNCTIONS
  USE ED_SETUP
  USE ED_SECTOR
  USE ED_HAMILTONIAN
  USE ED_AUX_FUNX

  implicit none
  private


  public :: build_chi_exct

  integer                      :: istate,iorb,jorb,ispin,jspin
  integer                      :: isector,jsector,ksector
  real(8),allocatable          :: vvinit(:),vvinit_tmp(:)
  real(8),allocatable          :: alfa_(:),beta_(:)
  integer                      :: ialfa
  integer                      :: jalfa
  integer                      :: ipos,jpos
  integer                      :: i,j,k
  real(8)                      :: sgn,norm2
  real(8),dimension(:),pointer :: state_cvec
  real(8)                      :: state_e

contains


  !+------------------------------------------------------------------+
  !                            EXCITON
  !PURPOSE  : Evaluate the Exciton susceptibility \Chi_exct for a 
  ! \chi_ab = <O*_a(\tau)O_b(0)>
  ! a/=b
  ! Singlet: \sum_\sigma <C^+_{a\sigma}C_{b\sigma} 
  ! Triplet: \sum_{\sigma\rho} C^+_{a\sigma} \tau_{\sigma\rho} C_{b\rho}
  !+------------------------------------------------------------------+
  subroutine build_chi_exct()
    if(Norb>1)then       
       write(LOGfile,"(A)")"Get impurity exciton Chi:"
       do iorb=1,Norb
          do jorb=iorb+1,Norb
             write(LOGfile,"(A)")"Get singlet Chi_exct_l"//reg(txtfy(iorb))//reg(txtfy(jorb))
             if(MPIMASTER)call start_timer()
             call lanc_ed_build_exctChi_singlet(iorb,jorb)
             if(MPIMASTER)call stop_timer(unit=LOGfile)
             !
             write(LOGfile,"(A)")"Get triplet Chi_exct_l"//reg(txtfy(iorb))//reg(txtfy(jorb))
             if(MPIMASTER)call start_timer()
             call lanc_ed_build_exctChi_tripletXY(iorb,jorb)
             call lanc_ed_build_exctChi_tripletZ(iorb,jorb)
             if(MPIMASTER)call stop_timer(unit=LOGfile)
             !
             exctChi_w(0:,jorb,iorb,:)   = exctChi_w(0:,iorb,jorb,:)
             exctChi_tau(0:,jorb,iorb,:) = exctChi_tau(0:,iorb,jorb,:)
             exctChi_iv(0:,jorb,iorb,:)  = exctChi_iv(0:,iorb,jorb,:)
          end do
       end do
    endif
  end subroutine build_chi_exct




  ! \chi_ab  = <Delta*_ab(\tau)Delta_ab(0)>
  !\Delta_ab = \sum_\sigma C^+_{a\sigma}C_{b\sigma}
  subroutine lanc_ed_build_exctChi_singlet(iorb,jorb)
    integer      :: iorb,jorb
    type(sector) :: sectorI,sectorK
    !
    if(ed_total_ud)then
       ialfa = 1
       jalfa = 1
       ipos  = iorb
       jpos  = jorb
    else
       write(LOGfile,"(A)")"ED_CHI_EXCITION warning: can not evaluate \Chi_exct_singlet with ed_total_ud=F"
       return
    endif
    !
    do istate=1,state_list%size
       isector    =  es_return_sector(state_list,istate)
       state_e    =  es_return_energy(state_list,istate)
#ifdef _MPI
       if(MpiStatus)then
          state_cvec => es_return_cvector(MpiComm,state_list,istate)
       else
          state_cvec => es_return_cvector(state_list,istate)
       endif
#else
       state_cvec => es_return_cvector(state_list,istate)
#endif
       !
       !C^+_as C_bs => jsector == isector
       if(MpiMaster)then
          call build_sector(isector,sectorI)
          if(ed_verbose>=3)write(LOGfile,"(A,I6,20I4)")&
               'Apply \sum_s C^+_as.C_bs :',isector,sectorI%Nups,sectorI%Ndws
          allocate(vvinit(sectorI%Dim))     ;  vvinit=0d0
       else
          allocate(vvinit(1));vvinit=0.d0
       endif
       !
       ksector = getCsector(jalfa,2,isector)       
       if(ksector/=0)then
          if(MpiMaster)then
             call build_sector(ksector,sectorK)
             allocate(vvinit_tmp(sectorK%Dim)) ;  vvinit_tmp=0d0
             !C_b,up|gs>=|tmp>
             do i=1,sectorI%Dim
                call apply_op_C(i,k,sgn,jpos,jalfa,2,sectorI,sectorK)
                if(sgn==0.OR.k==0)cycle
                vvinit_tmp(k) = sgn*state_cvec(i)
             enddo
             !C^+_a,up|tmp>=|vvinit>
             do k=1,sectorK%Dim
                call apply_op_CDG(k,i,sgn,ipos,ialfa,2,sectorK,sectorI)
                if(sgn==0.OR.k==0)cycle
                vvinit(i) = sgn*vvinit_tmp(k)
             enddo
             deallocate(vvinit_tmp)
             call delete_sector(sectorK)
          endif
       endif
       ksector = getCsector(jalfa,1,isector)
       if(ksector/=0)then
          if(MpiMaster)then
             call build_sector(ksector,sectorK)
             allocate(vvinit_tmp(sectorK%Dim)) ;  vvinit_tmp=0d0
             !C_b,dw|gs>=|tmp>
             do i=1,sectorI%Dim
                call apply_op_C(i,k,sgn,jpos,jalfa,1,sectorI,sectorK)
                if(sgn==0.OR.k==0)cycle
                vvinit_tmp(k) = sgn*state_cvec(i)
             enddo
             !C^+_a,dw|tmp>=|vvinit>
             do k=1,sectorK%Dim
                call apply_op_CDG(k,i,sgn,ipos,ialfa,1,sectorK,sectorI)
                if(sgn==0.OR.k==0)cycle
                vvinit(i) = vvinit(i) + sgn*vvinit_tmp(k)
             enddo
             deallocate(vvinit_tmp)
             call delete_sector(sectorK)
          endif
       endif
       !
       if(MpiMaster)call delete_sector(sectorI)
       !
       call tridiag_Hv_sector(isector,vvinit,alfa_,beta_,norm2)
       call add_to_lanczos_exctChi(norm2,state_e,alfa_,beta_,iorb,jorb,0)
       deallocate(alfa_,beta_)
       if(allocated(vvinit))deallocate(vvinit)
       !
#ifdef _MPI
       if(MpiStatus)then
          if(associated(state_cvec))deallocate(state_cvec)
       else
          if(associated(state_cvec))nullify(state_cvec)
       endif
#else
       if(associated(state_cvec))nullify(state_cvec)
#endif
       !
    enddo
    return
  end subroutine lanc_ed_build_exctChi_singlet







  ! \chi_ab  = <Z_ab(\tau)Z_ab(0)>
  !Z_ab = \sum_sp C^+_{as}.tau^z_{sp}.C_{bp}
  subroutine lanc_ed_build_exctChi_tripletZ(iorb,jorb)
    integer      :: iorb,jorb
    type(sector) :: sectorI,sectorK,sectorG
    !
    if(ed_total_ud)then
       ialfa = 1
       jalfa = 1
       ipos  = iorb
       jpos  = jorb
    else
       write(LOGfile,"(A)")"ED_CHI_EXCITION warning: can not evaluate \Chi_exc=t_triplet with ed_total_ud=F"
       return
    endif
    !
    do istate=1,state_list%size
       isector    =  es_return_sector(state_list,istate)
       state_e    =  es_return_energy(state_list,istate)
#ifdef _MPI
       if(MpiStatus)then
          state_cvec => es_return_cvector(MpiComm,state_list,istate)
       else
          state_cvec => es_return_cvector(state_list,istate)
       endif
#else
       state_cvec => es_return_cvector(state_list,istate)
#endif
       !
       !Z - Component:
       !Z_{ab}= C^+_{a,up}C_{b,up} - C^+_{a,dw}C_{b,dw}
       if(MpiMaster)then
          call build_sector(isector,sectorI)
          allocate(vvinit(sectorI%Dim))     ;  vvinit=0d0
          if(ed_verbose>=3)write(LOGfile,"(A,I6,20I4)")&
               'Apply \sum_s C^+_as.C_bs :',isector,sectorI%Nups,sectorI%Ndws
       else
          allocate(vvinit(1));vvinit=0.d0
       endif
       !Intermediate sectors:
       ksector = getCsector(jalfa,2,isector)
       if(ksector/=0)then
          if(MpiMaster)then
             call build_sector(ksector,sectorK)
             allocate(vvinit_tmp(sectorK%Dim)) ;  vvinit_tmp=0d0
             !C_b,dw|gs>=|tmp>
             do i=1,sectorI%Dim
                call apply_op_C(i,k,sgn,jpos,jalfa,2,sectorI,sectorK)
                if(sgn==0.OR.k==0)cycle
                vvinit_tmp(k) = sgn*state_cvec(i)
             enddo
             !C^+_a,dw|tmp>=|vvinit>
             do k=1,sectorK%Dim
                call apply_op_CDG(k,i,sgn,ipos,ialfa,2,sectorK,sectorI)
                if(sgn==0.OR.k==0)cycle
                vvinit(i) = sgn*vvinit_tmp(k)
             enddo
             deallocate(vvinit_tmp)
             call delete_sector(sectorK)
          endif
       endif
       ksector = getCsector(jalfa,1,isector)
       if(ksector/=0)then
          if(MpiMaster)then
             call build_sector(ksector,sectorK)
             allocate(vvinit_tmp(sectorK%Dim)) ;  vvinit_tmp=0d0
             !C_b,up|gs>=|tmp>
             do i=1,sectorI%Dim
                call apply_op_C(i,k,sgn,jpos,jalfa,1,sectorI,sectorK)
                if(sgn==0.OR.k==0)cycle
                vvinit_tmp(k) = sgn*state_cvec(i)
             enddo
             !C^+_a,up|tmp>=|vvinit>
             do k=1,sectorK%Dim
                call apply_op_CDG(k,i,sgn,ipos,ialfa,1,sectorK,sectorI)
                if(sgn==0.OR.k==0)cycle
                vvinit(i) =  sgn*vvinit_tmp(k) - vvinit(i)
             enddo
             deallocate(vvinit_tmp)
             call delete_sector(sectorK)
          endif
       endif
       if(MpiMaster)call delete_sector(sectorI)
       !
       call tridiag_Hv_sector(isector,vvinit,alfa_,beta_,norm2)
       call add_to_lanczos_exctChi(norm2,state_e,alfa_,beta_,iorb,jorb,2)
       deallocate(alfa_,beta_)
       if(allocated(vvinit))deallocate(vvinit)
       !
#ifdef _MPI
       if(MpiStatus)then
          if(associated(state_cvec))deallocate(state_cvec)
       else
          if(associated(state_cvec))nullify(state_cvec)
       endif
#else
       if(associated(state_cvec))nullify(state_cvec)
#endif
       !
    enddo
    return
  end subroutine lanc_ed_build_exctChi_tripletZ





  ! \chi_ab  = <O_ab(\tau)O_ab(0)>
  ! O_ab = \sum_sp C^+_{as}.tau^o_{sp}.C_{bp} with o=X,Y
  ! O_ab|0> X:=   [C^+_{a,up}C_{b,dw} + C^+_{a,dw}C_{b,up}]|0>
  !         Y:= -i[C^+_{a,up}C_{b,dw} - C^+_{a,dw}C_{b,up}]|0>
  !         X:=   [P_{up,dw} +  P_{dw,up}]|0> = |v> + |w> 
  !         X:= -i[P_{up,dw} -  P_{dw,up}]|0> = |v> + |w>
  ! If |0>\in\SS(N_up,N_dw) => |v>\in\SS(N_up+1,N_dw-1), |w>\in\SS(N_up-1,N_dw+1)
  ! so that the sum |v> + |w> can not be accumulated onto a single vector |vvinit>
  ! Yet, we can recast the \Chi_ab expression in:
  ! \chi_ab = (<w| + <v|)[z-H]^{-1}(|v> + |w>)
  ! the direct terms: <v|[z-H]^{-1}|v> and <w|[z-H]^{-1}|w> are evaluated as usual.
  ! the mixed terms: <v|[z-H]^{-1}|w> and <w|[z-H]^{-1}|v> are indeed null.
  ! Proof:
  ! |v> and |w> belong to different sectors. H has a sector-block structure and so
  ! does its inverse H^{-1} made of the inverse of each block. 
  ! The expected values <v|H^{-1}|w> are taken across different sectors, but unless
  ! spin-conservation is broken these sectors are not connected and, as such, these
  ! numbers have to be zero.
  ! Note that while H can have a sparse nature, its inverse in general does not.
  ! For this reason the same argument as above DOES NOT apply to the Z case as
  ! in that case |v> and |w> belong to the same sector (the same as |0>) and the
  ! mixed term is in general non null.
  subroutine lanc_ed_build_exctChi_tripletXY(iorb,jorb)
    integer      :: iorb,jorb
    type(sector) :: sectorI,sectorK,sectorJ
    !
    if(ed_total_ud)then
       ialfa = 1
       jalfa = 1
       ipos  = iorb
       jpos  = jorb
    else
       write(LOGfile,"(A)")"ED_CHI_EXCITION warning: can not evaluate \Chi_exc=t_triplet with ed_total_ud=F"
       return
    endif
    !
    do istate=1,state_list%size
       isector    =  es_return_sector(state_list,istate)
       state_e    =  es_return_energy(state_list,istate)
#ifdef _MPI
       if(MpiStatus)then
          state_cvec => es_return_cvector(MpiComm,state_list,istate)
       else
          state_cvec => es_return_cvector(state_list,istate)
       endif
#else
       state_cvec => es_return_cvector(state_list,istate)
#endif
       !
       !X - Component == Y -Component 
       !X_{ab}= C^+_{a,up}C_{b,dw} + C^+_{a,dw}C_{b,up}
       !
       !C^+_{a,dw}C_{b,up}:
       ksector = getCsector(jalfa,1,isector)
       if(ksector/=0)then
          jsector = getCDGsector(ialfa,2,ksector)
          if(jsector/=0)then
             if(MpiMaster)then
                call build_sector(isector,sectorI)
                call build_sector(ksector,sectorK)
                call build_sector(jsector,sectorJ)
                if(ed_verbose>=3)write(LOGfile,"(A,I6,20I4)")&
                     'Apply C^+_{a,dw}C_{b,up}: :',isector,sectorI%Nups,sectorI%Ndws
                allocate(vvinit_tmp(sectorK%Dim)) ;  vvinit_tmp=0d0
                allocate(vvinit(sectorJ%Dim))     ;  vvinit=0d0
                !C_{b,up}|0>=|tmp>
                do i=1,sectorI%Dim
                   call apply_op_C(i,k,sgn,jpos,jalfa,1,sectorI,sectorK)
                   if(sgn==0.OR.k==0)cycle
                   vvinit_tmp(k) = sgn*state_cvec(i)
                enddo
                !C^+_{a,dw}|tmp>=|vvinit>
                do k=1,sectorK%Dim
                   call apply_op_CDG(k,i,sgn,ipos,ialfa,2,sectorK,sectorJ)
                   if(sgn==0.OR.k==0)cycle
                   vvinit(i) = sgn*vvinit_tmp(k)
                enddo
                deallocate(vvinit_tmp)
                call delete_sector(sectorI)
                call delete_sector(sectorK)
                call delete_sector(sectorJ)
             else
                allocate(vvinit(1));vvinit=0.d0
             endif
             !
             call tridiag_Hv_sector(jsector,vvinit,alfa_,beta_,norm2)
             call add_to_lanczos_exctChi(norm2,state_e,alfa_,beta_,iorb,jorb,1)
             deallocate(alfa_,beta_)
             if(allocated(vvinit))deallocate(vvinit)
          endif
       endif
       !
       !C^+_{a,up}C_{b,dw}:
       ksector = getCsector(jalfa,2,isector)
       if(ksector/=0)then
          jsector = getCDGsector(ialfa,1,ksector)
          if(jsector/=0)then
             if(MpiMaster)then
                call build_sector(isector,sectorI)
                call build_sector(ksector,sectorK)
                call build_sector(jsector,sectorJ)
                if(ed_verbose>=3)write(LOGfile,"(A,I6,20I4)")&
                     'Apply C^+_{a,dw}C_{b,up}: :',isector,sectorI%Nups,sectorI%Ndws
                allocate(vvinit_tmp(sectorK%Dim)) ;  vvinit_tmp=0d0
                allocate(vvinit(sectorJ%Dim))     ;  vvinit=0d0
                !C_{b,dw}|0>=|tmp>
                do i=1,sectorI%Dim
                   call apply_op_C(i,k,sgn,jpos,jalfa,2,sectorI,sectorK)
                   if(sgn==0.OR.k==0)cycle
                   vvinit_tmp(k) = sgn*state_cvec(i)
                enddo
                !C^+_{a,up}|tmp>=|vvinit>
                do k=1,sectorK%Dim
                   call apply_op_CDG(k,i,sgn,ipos,ialfa,1,sectorK,sectorJ)
                   if(sgn==0.OR.k==0)cycle
                   vvinit(i) = sgn*vvinit_tmp(k)
                enddo
                deallocate(vvinit_tmp)
                call delete_sector(sectorI)
                call delete_sector(sectorK)
                call delete_sector(sectorJ)
             else
                allocate(vvinit(1));vvinit=0.d0
             endif
             !
             call tridiag_Hv_sector(jsector,vvinit,alfa_,beta_,norm2)
             call add_to_lanczos_exctChi(norm2,state_e,alfa_,beta_,iorb,jorb,1)
             deallocate(alfa_,beta_)
             if(allocated(vvinit))deallocate(vvinit)
          endif
       endif
       !
       !
#ifdef _MPI
       if(MpiStatus)then
          if(associated(state_cvec))deallocate(state_cvec)
       else
          if(associated(state_cvec))nullify(state_cvec)
       endif
#else
       if(associated(state_cvec))nullify(state_cvec)
#endif
       !
    enddo
    return
  end subroutine lanc_ed_build_exctChi_tripletXY







  subroutine add_to_lanczos_exctChi(vnorm2,Ei,alanc,blanc,iorb,jorb,indx)
    integer                                    :: iorb,jorb,indx
    real(8)                                    :: pesoF,pesoAB,pesoBZ,peso,vnorm2  
    real(8)                                    :: Ei,Ej,Egs,de
    integer                                    :: nlanc
    real(8),dimension(:)                       :: alanc
    real(8),dimension(size(alanc))             :: blanc 
    real(8),dimension(size(alanc),size(alanc)) :: Z
    real(8),dimension(size(alanc))             :: diag,subdiag
    integer                                    :: i,j,ierr
    complex(8)                                 :: iw,chisp
    !
    if(vnorm2==0)return
    !
    Egs = state_list%emin       !get the gs energy
    !
    Nlanc = size(alanc)
    !
    if(.not.any([0,1,2]==indx))stop "add_to_lanczos_exctChi ERROR: indx/=any[0,1,2]"
    !
    pesoF  = vnorm2/zeta_function 
    pesoBZ = 1d0
    if(finiteT)pesoBZ = exp(-beta*(Ei-Egs))
    !
#ifdef _MPI
    if(MpiStatus)then
       call Bcast_MPI(MpiComm,alanc)
       call Bcast_MPI(MpiComm,blanc)
    endif
#endif
    diag(1:Nlanc)    = alanc(1:Nlanc)
    subdiag(2:Nlanc) = blanc(2:Nlanc)
    call eigh(diag(1:Nlanc),subdiag(2:Nlanc),Ev=Z(:Nlanc,:Nlanc))
    !
    do j=1,nlanc
       Ej     = diag(j)
       dE     = Ej-Ei
       pesoAB = Z(1,j)*Z(1,j)
       peso   = pesoF*pesoAB*pesoBZ
       ! the correct behavior for beta*dE << 1 is recovered only by assuming that v_n is still finite
       ! beta*dE << v_n for v_n--> 0 slower. First limit beta*dE--> 0 and only then v_n -->0.
       ! This ensures that the correct null contribution is obtained.
       ! So we impose that: if (beta*dE is larger than a small qty) we sum up the contribution, else
       ! we do not include the contribution (because we are in the situation described above).
       ! For the real-axis case this problem is circumvented by the usual i*0+ = xi*eps
       if(beta*dE > 1d-3)exctChi_iv(indx,iorb,jorb,0)=exctChi_iv(indx,iorb,jorb,0) + peso*2*(1d0-exp(-beta*dE))/dE 
       do i=1,Lmats
          exctChi_iv(indx,iorb,jorb,i)=exctChi_iv(indx,iorb,jorb,i) + &
               peso*(1d0-exp(-beta*dE))*2d0*dE/(vm(i)**2+dE**2)
       enddo
       do i=0,Ltau
          exctChi_tau(indx,iorb,jorb,i)=exctChi_tau(indx,iorb,jorb,i) + exp(-tau(i)*dE)*peso
       enddo
       do i=1,Lreal
          exctChi_w(indx,iorb,jorb,i)=exctChi_w(indx,iorb,jorb,i) - &
               peso*(1d0-exp(-beta*dE))*(1d0/(dcmplx(vr(i),eps) - dE) - 1d0/(dcmplx(vr(i),eps) + dE))
       enddo
    enddo
    !
  end subroutine add_to_lanczos_exctChi



END MODULE ED_CHI_EXCT





























! I think this is wrong because it neglects mixed terms in the
! trace, e.g. < O^+_\up(t) O_\dw(0)>
! with O_\sigma = C^+_{a\sigma}C_{b\sigma} which corresponds to
! off-diagonal terms of the resolvant (z-H)^-1. 
!   ! \chi_ab  = <S_ab(\tau)S_ab(0)>
!   !S_ab = \sum_\sigma C^+_{a\sigma}C_{b\sigma}
!   subroutine lanc_ed_build_exctChi_singlet_(iorb,jorb)
!     integer      :: iorb,jorb
!     type(sector) :: sectorI,sectorK,sectorG
!     !
!     if(ed_total_ud)then
!        ialfa = 1
!        jalfa = 1
!        ipos  = iorb
!        jpos  = jorb
!     else
!        write(LOGfile,"(A)")"ED_CHI_EXCITION warning: can not evaluate \Chi_exct_singlet with ed_total_ud=F"
!        return
!     endif
!     !
!     do istate=1,state_list%size
!        isector    =  es_return_sector(state_list,istate)
!        state_e    =  es_return_energy(state_list,istate)
! #ifdef _MPI
!        if(MpiStatus)then
!           state_cvec => es_return_cvector(MpiComm,state_list,istate)
!        else
!           state_cvec => es_return_cvector(state_list,istate)
!        endif
! #else
!        state_cvec => es_return_cvector(state_list,istate)
! #endif
!        !
!        !C^+_as C_bs => jsector == isector
!        do ispin=1,2
!           ksector = getCsector(jalfa,ispin,isector)
!           if(ksector/=0)then
!              if(MpiMaster)then
!                 call build_sector(isector,sectorI)
!                 call build_sector(ksector,sectorK)
!                 if(ed_verbose>=3)write(LOGfile,"(A,I6,20I4)")&
!                      'Apply \sum_s C^+_as.C_bs :',isector,sectorI%Nups,sectorI%Ndws
!                 allocate(vvinit(sectorI%Dim))     ;  vvinit=0d0
!                 allocate(vvinit_tmp(sectorK%Dim)) ;  vvinit_tmp=0d0
!                 !C_b,s|gs>=|tmp>
!                 do i=1,sectorI%Dim
!                    call apply_op_C(i,k,sgn,jpos,jalfa,ispin,sectorI,sectorK)
!                    if(sgn==0.OR.k==0)cycle
!                    vvinit_tmp(k) = sgn*state_cvec(i)
!                 enddo
!                 !C^+_a,s|tmp>=|vvinit>
!                 do k=1,sectorK%Dim
!                    call apply_op_CDG(k,i,sgn,ipos,ialfa,ispin,sectorK,sectorI)
!                    if(sgn==0.OR.k==0)cycle
!                    vvinit(i) = sgn*vvinit_tmp(k)
!                 enddo
!                 deallocate(vvinit_tmp)
!                 call delete_sector(sectorK)
!                 call delete_sector(sectorI)
!              else
!                 allocate(vvinit(1));vvinit=0.d0
!              endif
!              !
!              call tridiag_Hv_sector(isector,vvinit,alfa_,beta_,norm2)
!              call add_to_lanczos_exctChi(norm2,state_e,alfa_,beta_,iorb,jorb,1)
!              deallocate(alfa_,beta_)
!              if(allocated(vvinit))deallocate(vvinit)
!           endif
!        enddo
!        !
! #ifdef _MPI
!        if(MpiStatus)then
!           if(associated(state_cvec))deallocate(state_cvec)
!        else
!           if(associated(state_cvec))nullify(state_cvec)
!        endif
! #else
!        if(associated(state_cvec))nullify(state_cvec)
! #endif
!        !
!     enddo
!     return
!   end subroutine lanc_ed_build_exctChi_singlet_
