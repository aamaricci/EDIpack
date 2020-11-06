MODULE ED_GF_PHONON
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
  USE ED_HAMILTONIAN
  USE ED_AUX_FUNX
  implicit none
  private

  public :: build_gf_phonon

  integer                               :: istate
  integer                               :: isector
  integer                               :: idim,idimUP,idimDW
  real(8),allocatable                   :: vvinit(:)
  real(8),allocatable                   :: alfa_(:),beta_(:)
  integer                               :: i,j,m
  integer                               :: iph,i_el
  real(8)                               :: norm2
  real(8),dimension(:),pointer          :: state_cvec
  real(8)                               :: state_e


contains

  !+------------------------------------------------------------------+
  !                            PHONON
  !PURPOSE  : Evaluate the Green's function of the phonons 
  ! D = -<x(\tau)x(0)> with x = (b + b^+)
  ! note: as x is hermitian, particle and holes contributions
  ! are identical so work out only one lanczos tridiag. work out the 
  ! reduction for both values of isign in the same call.
  !+------------------------------------------------------------------+
  subroutine build_gf_phonon()
    write(LOGfile,"(A)")"Get phonon Green function:"
    if(MPIMASTER)call start_timer()
    select case(ed_diag_type)
    case default
       call lanc_build_gf_phonon_main()
    case ("full")
       call full_build_gf_phonon_main()
    end select
    if(MPIMASTER)call stop_timer(unit=LOGfile)
    !
  end subroutine build_gf_phonon


  !################################################################
  !################################################################
  !################################################################
  !################################################################


  subroutine lanc_build_gf_phonon_main()
    integer,dimension(Ns_Ud)    :: iDimUps,iDimDws
    integer                     :: Nups(Ns_Ud)
    integer                     :: Ndws(Ns_Ud)
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
       call get_Nup(isector,Nups)
       call get_Ndw(isector,Ndws)
       if(MpiMaster.AND.ed_verbose>=3)write(LOGfile,"(A,I6,20I4)")'From sector:',isector,Nups,Ndws
       !
       idim = getdim(isector)
       call get_DimUp(isector,iDimUps)
       call get_DimDw(isector,iDimDws)
       iDimUp = product(iDimUps)
       iDimDw = product(iDimDws)
       !
       if(MpiMaster)then
          if(ed_verbose==3)write(LOGfile,"(A,I12)")'Apply x:',isector
          !
          allocate(vvinit(idim));vvinit=0d0
          !
          do i=1,iDim
             iph = (i-1)/(iDimUp*iDimDw) + 1
             i_el = mod(i-1,iDimUp*iDimDw) + 1
             !
             !apply destruction operator
             if(iph>1) then
                j = i_el + ((iph-1)-1)*iDimUp*iDimDw
                vvinit(j) = vvinit(j) + sqrt(dble(iph-1))*state_cvec(i)
             endif
             !
             !apply creation operator
             if(iph<DimPh) then
                j = i_el + ((iph+1)-1)*iDimUp*iDimDw
                vvinit(j) = vvinit(j) + sqrt(dble(iph))*state_cvec(i)
             endif
          enddo
       else
          allocate(vvinit(1));vvinit=0.d0
       endif
       !
       call tridiag_Hv_sector(isector,vvinit,alfa_,beta_,norm2)
       call add_to_lanczos_phonon(norm2,state_e,alfa_,beta_)
       deallocate(alfa_,beta_)
       if(allocated(vvinit))deallocate(vvinit)
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
  end subroutine lanc_build_gf_phonon_main


  !################################################################


  subroutine add_to_lanczos_phonon(vnorm2,Ei,alanc,blanc)
    real(8)                                    :: vnorm2,Ei,Ej,Egs,pesoF,pesoAB,pesoBZ,de,peso
    integer                                    :: nlanc
    real(8),dimension(:)                       :: alanc
    real(8),dimension(size(alanc))             :: blanc 
    real(8),dimension(size(alanc),size(alanc)) :: Z
    real(8),dimension(size(alanc))             :: diag,subdiag
    integer                                    :: i,j
    complex(8)                                 :: iw
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
       if(beta*dE > 1d-3)impDmats_ph(0)=impDmats_ph(0) - peso*2*(1d0-exp(-beta*dE))/dE 
       do i=1,Lmats
          impDmats_ph(i)=impDmats_ph(i) - peso*(1d0-exp(-beta*dE))*2d0*dE/(vm(i)**2+dE**2)
       enddo
       do i=1,Lreal
          impDreal_ph(i)=impDreal_ph(i) + peso*(1d0-exp(-beta*dE))*(1d0/(dcmplx(vr(i),eps) - dE) - 1d0/(dcmplx(vr(i),eps) + dE))
       enddo
    enddo
  end subroutine add_to_lanczos_phonon


  !################################################################
  !################################################################
  !################################################################
  !################################################################


  subroutine full_build_gf_phonon_main()
    integer,dimension(Ns_Ud)    :: iDimUps,iDimDws
    integer,dimension(2,Ns_Orb) :: Nud
    integer                     :: Iud(2)
    real(8)                     :: C1
    integer                     :: i,j,ll,ll2,isector
    integer                     :: idim
    real(8)                     :: Ei,Ej,peso
    real(8)                     :: expterm,de
    complex(8)                  :: iw 
    !
    do isector=1,Nsectors !loop over <i| total particle number
       call eta(isector,Nsectors,LOGfile)
       idim = getdim(isector)
       call get_DimUp(isector,iDimUps)
       call get_DimDw(isector,iDimDws)
       iDimUp = product(iDimUps)
       iDimDw = product(iDimDws)
       !
       do i=1,idim 
          do j=1,idim
             C1=0d0
             expterm=exp(-beta*espace(isector)%e(i))+exp(-beta*espace(isector)%e(j))
             if(expterm<cutoff)cycle
             do ll=1,idim
                iph = (ll-1)/(iDimUp*iDimDw) + 1
                i_el = mod(ll-1,iDimUp*iDimDw) + 1
                !
                !creation operator
                if(iph<DimPh) then
                   ll2 = i_el + ((iph+1)-1)*iDimUp*iDimDw
                   C1 = C1 + espace(isector)%M(ll2,i)*sqrt(dble(iph))*espace(isector)%M(ll,j)
                endif
                !
                !destruction operator
                if(iph>1) then
                   ll2 = i_el + ((iph-1)-1)*iDimUp*iDimDw
                   C1 = C1 + espace(isector)%M(ll2,i)*sqrt(dble(iph-1))*espace(isector)%M(ll,j)
                endif
             enddo
             Ei=espace(isector)%e(i)
             Ej=espace(isector)%e(j)
             de=Ei-Ej
             peso = C1**2.d0/zeta_function
             !
             !Matsubara (bosonic) frequency
             if(beta*dE > 1d-3)impDmats_ph(0)=impDmats_ph(0) - peso*2*exp(-beta*Ej)*(1d0-exp(-beta*dE))/dE
             do m=1,Lmats
                impDmats_ph(m)=impDmats_ph(m) - peso*exp(-beta*Ej)*2*dE/(vm(m)**2 + de**2)
             enddo
             !
             !Real-frequency: Retarded = Commutator = response function
             do m=1,Lreal
                iw=dcmplx(vr(m),eps)
                impDreal_ph(m)=impDreal_ph(m) + peso*(exp(-beta*Ei) - exp(-beta*Ej))/(iw+de)
             enddo
             !
          enddo
       enddo
    enddo
  end subroutine full_build_gf_phonon_main


END MODULE ED_GF_PHONON













