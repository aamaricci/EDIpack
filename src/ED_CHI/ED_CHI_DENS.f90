MODULE ED_CHI_DENS
  USE SF_CONSTANTS, only:one,xi,zero,pi
  USE SF_TIMER  
  USE SF_IOTOOLS, only: str,reg,txtfy
  USE SF_LINALG,  only: inv,eigh,eye
  USE ED_INPUT_VARS
  USE ED_VARS_GLOBAL
  USE ED_EIGENSPACE
  USE ED_BATH
  USE ED_BATH_FUNCTIONS
  USE ED_SETUP
  USE ED_HAMILTONIAN
  USE ED_AUX_FUNX

  implicit none
  private


  public :: build_chi_dens

  integer                      :: istate,iorb,jorb,ispin,jspin
  integer                      :: isector
  real(8),allocatable          :: vvinit(:)
  real(8),allocatable          :: alfa_(:),beta_(:)
  integer                      :: ialfa
  integer                      :: jalfa
  integer                      :: ipos,jpos
  integer                      :: i,j
  integer                      :: iph,i_el
  real(8)                      :: sgn,norm2
  real(8),dimension(:),pointer :: state_cvec
  real(8)                      :: state_e


contains


  !+------------------------------------------------------------------+
  !                            DENS
  !PURPOSE  : Evaluate the Dens susceptibility \Chi_dens for a 
  ! \chi_ab = <n_a(\tau)n_b(0)>
  !+------------------------------------------------------------------+
  subroutine build_chi_dens()
    write(LOGfile,"(A)")"Get impurity dens Chi:"
    do iorb=1,Norb
       write(LOGfile,"(A)")"Get Chi_dens_l"//reg(txtfy(iorb))
       if(MPIMASTER)call start_timer()
       select case(ed_diag_type)
       case default
          call lanc_ed_build_densChi_diag(iorb)
       case ("full")
          call full_ed_build_densChi_main(iorb,iorb)
       end select
       if(MPIMASTER)call stop_timer(unit=LOGfile)
    enddo
    !
    if(Norb>1)then
       do iorb=1,Norb
          do jorb=iorb+1,Norb
             write(LOGfile,"(A)")"Get Chi_dens_mix_l"//reg(txtfy(iorb))//reg(txtfy(jorb))
             if(MPIMASTER)call start_timer()
             select case(ed_diag_type)
             case default
                call lanc_ed_build_densChi_mix(iorb,jorb)
             case ("full")
                call full_ed_build_densChi_main(iorb,jorb)
             end select
             if(MPIMASTER)call stop_timer(unit=LOGfile)
          end do
       end do
       !
       !
       do iorb=1,Norb
          do jorb=iorb+1,Norb
             select case(ed_diag_type)
             case default
                densChi_w(iorb,jorb,:)   = 0.5d0*(densChi_w(iorb,jorb,:) - densChi_w(iorb,iorb,:) - densChi_w(jorb,jorb,:))
                densChi_tau(iorb,jorb,:) = 0.5d0*(densChi_tau(iorb,jorb,:) - densChi_tau(iorb,iorb,:) - densChi_tau(jorb,jorb,:))
                densChi_iv(iorb,jorb,:)  = 0.5d0*(densChi_iv(iorb,jorb,:) - densChi_iv(iorb,iorb,:) - densChi_iv(jorb,jorb,:))
                !
             case ("full")
                ! The previous calculation is not needed in the FULL ED case
             end select
             !
             densChi_w(jorb,iorb,:)   = densChi_w(iorb,jorb,:)
             densChi_tau(jorb,iorb,:) = densChi_tau(iorb,jorb,:)
             densChi_iv(jorb,iorb,:)  = densChi_iv(iorb,jorb,:)
          enddo
       enddo
    endif
    !
  end subroutine build_chi_dens






  !################################################################
  !################################################################
  !################################################################
  !################################################################






  subroutine lanc_ed_build_densChi_diag(iorb)
    integer                     :: iorb
    type(sector)                :: sectorI,sectorJ
    !
    if(ed_total_ud)then
       ialfa = 1
       ipos  = iorb
    else
       ialfa = iorb
       ipos  = 1
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
       if(MpiMaster)then
          call build_sector(isector,sectorI)
          if(ed_verbose>=3)write(LOGfile,"(A,I6,20I4)")&
               'From sector  :',isector,sectorI%Nups,sectorI%Ndws
          if(ed_verbose==3)write(LOGfile,"(A,I12)")'Apply N:',isector
          allocate(vvinit(sectorI%Dim)) ; vvinit=zero
          do i=1,sectorI%Dim
             call apply_op_N(i,sgn,ipos,ialfa,sectorI)
             vvinit(i) = sgn*state_cvec(i)
          enddo
          call delete_sector(sectorI)
       else
          allocate(vvinit(1));vvinit=0.d0
       endif
       !
       call tridiag_Hv_sector(isector,vvinit,alfa_,beta_,norm2)
       call add_to_lanczos_densChi(norm2,state_e,alfa_,beta_,iorb,iorb)
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
    enddo
    return
  end subroutine lanc_ed_build_densChi_diag



  !################################################################



  subroutine lanc_ed_build_densChi_mix(iorb,jorb)
    integer                     :: iorb,jorb
    type(sector)                :: sectorI,sectorJ
    real(8)                     :: Niorb,Njorb
    !
    if(ed_total_ud)then
       ialfa = 1
       jalfa = 1
       ipos  = iorb
       jpos  = jorb
    else
       ialfa = iorb
       jalfa = jorb
       ipos  = 1
       jpos  = 1
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
       !EVALUATE (N_jorb + N_iorb)|gs> = N_jorb|gs> + N_iorb|gs>
       if(MpiMaster)then
          call build_sector(isector,sectorI)
          if(ed_verbose>=3)write(LOGfile,"(A,I6,20I4)")&
               'From sector  :',isector,sectorI%Nups,sectorI%Ndws
          if(ed_verbose==3)write(LOGfile,"(A,I15)")'Apply Na+Nb:',isector
          allocate(vvinit(sectorI%Dim)) ; vvinit=zero
          do i=1,sectorI%Dim
             call apply_op_N(i,Niorb,ipos,ialfa,sectorI)
             call apply_op_N(i,Njorb,jpos,jalfa,sectorI)
             sgn       = Niorb + Njorb
             vvinit(i) = sgn*state_cvec(i)
          enddo
          call delete_sector(sectorI)
       else
          allocate(vvinit(1));vvinit=0.d0
       endif
       !
       call tridiag_Hv_sector(isector,vvinit,alfa_,beta_,norm2)
       call add_to_lanczos_densChi(norm2,state_e,alfa_,beta_,iorb,jorb)
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
    enddo
    return
  end subroutine lanc_ed_build_densChi_mix




  !################################################################




  subroutine add_to_lanczos_densChi(vnorm2,Ei,alanc,blanc,iorb,jorb)
    integer                                    :: iorb,jorb
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
    Egs = state_list%emin       !get the gs energy
    !
    Nlanc = size(alanc)
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
       if(beta*dE > 1d-3)densChi_iv(iorb,jorb,0)=densChi_iv(iorb,jorb,0) + peso*2*(1d0-exp(-beta*dE))/dE 
       do i=1,Lmats
          densChi_iv(iorb,jorb,i)=densChi_iv(iorb,jorb,i) + peso*(1d0-exp(-beta*dE))*2d0*dE/(vm(i)**2+dE**2)
       enddo
       do i=0,Ltau
          densChi_tau(iorb,jorb,i)=densChi_tau(iorb,jorb,i) + exp(-tau(i)*dE)*peso
       enddo
       do i=1,Lreal
          densChi_w(iorb,jorb,i)=densChi_w(iorb,jorb,i) - &
               peso*(1d0-exp(-beta*dE))*(1d0/(dcmplx(vr(i),eps) - dE) - 1d0/(dcmplx(vr(i),eps) + dE))
       enddo
    enddo
    !
  end subroutine add_to_lanczos_densChi





  !################################################################
  !################################################################
  !################################################################
  !################################################################



  subroutine full_ed_build_densChi_main(iorb,jorb)
    integer                     :: iorb,jorb
    type(sector)                :: sectorI,sectorJ
    real(8)                     :: Chiorb,Chjorb,Niorb,Njorb
    integer                     :: i,j,ll,m,isector
    real(8)                     :: Ei,Ej,cc,peso,pesotot
    real(8)                     :: expterm,de,w0,it
    complex(8)                  :: iw 
    !
    !
    !Dens susceptibility \X(tau).
    !
    if(ed_total_ud)then
       ialfa = 1
       jalfa = 1
       ipos  = iorb
       jpos  = jorb
    else
       ialfa = iorb
       jalfa = jorb
       ipos  = 1
       jpos  = 1
    endif
    !
    do isector=1,Nsectors !loop over <i| total particle number
       call eta(isector,Nsectors,LOGfile)
       call build_sector(isector,sectorI)
       !
       do i=1,sectorI%Dim 
          do j=1,sectorI%Dim
             Chiorb=0d0
             Chjorb=0d0
             expterm=exp(-beta*espace(isector)%e(i))+exp(-beta*espace(isector)%e(j))
             if(expterm<cutoff)cycle
             do ll=1,sectorI%Dim
                call apply_op_N(i,Niorb,ipos,ialfa,sectorI)
                Chiorb   = Chiorb + espace(isector)%M(ll,i)*Niorb*espace(isector)%M(ll,j)
                call apply_op_N(i,Njorb,jpos,jalfa,sectorI)
                Chjorb   = Chjorb + espace(isector)%M(ll,i)*Njorb*espace(isector)%M(ll,j)
             enddo
             Ei=espace(isector)%e(i)
             Ej=espace(isector)%e(j)
             de=Ei-Ej
             peso = Chiorb*Chjorb/zeta_function
             !
             !Matsubara (bosonic) frequency
             if(beta*dE > 1d-3)densChi_iv(iorb,jorb,0)=densChi_iv(iorb,jorb,0) + peso*2*exp(-beta*Ej)*(1d0-exp(-beta*dE))/dE
             do m=1,Lmats
                densChi_iv(iorb,jorb,m)=densChi_iv(iorb,jorb,m)+ peso*exp(-beta*Ej)*2*dE/(vm(m)**2 + de**2)
             enddo
             !
             !Imaginary time: V
             do m=0,Ltau 
                it=tau(m)
                densChi_tau(iorb,jorb,m)=densChi_tau(iorb,jorb,m) + exp(-it*Ei)*exp(-(beta-it)*Ej)*peso
             enddo
             !
             !Real-frequency: Retarded = Commutator = response function
             do m=1,Lreal
                iw=dcmplx(vr(m),eps)
                densChi_w(iorb,jorb,m)=densChi_w(iorb,jorb,m)-peso*(exp(-beta*Ei) - exp(-beta*Ej))/(iw+de)
             enddo
             !
          enddo
       enddo
    enddo
  end subroutine full_ed_build_densChi_main


END MODULE ED_CHI_DENS
























