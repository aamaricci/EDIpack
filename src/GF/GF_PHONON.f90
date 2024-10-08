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
  USE ED_SECTOR
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
  real(8),dimension(:),allocatable      :: state_cvec
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
    if(MPIMASTER)call start_timer(unit=LOGfile)
    call lanc_build_gf_phonon_main()
    if(MPIMASTER)call stop_timer
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
          call es_return_cvector(MpiComm,state_list,istate,state_cvec)
       else
          call es_return_cvector(state_list,istate,state_cvec)
       endif
#else
       call es_return_cvector(state_list,istate,state_cvec)
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
       if(allocated(state_cvec))deallocate(state_cvec)
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
       if(beta*dE > 1d-3)impDmats(0)=impDmats(0) - peso*2*(1d0-exp(-beta*dE))/dE 
       do i=1,Lmats
          impDmats(i)=impDmats(i) - peso*(1d0-exp(-beta*dE))*2d0*dE/(vm(i)**2+dE**2)
       enddo
       do i=1,Lreal
          impDreal(i)=impDreal(i) + peso*(1d0-exp(-beta*dE))*(1d0/(dcmplx(vr(i),eps) - dE) - 1d0/(dcmplx(vr(i),eps) + dE))
       enddo
    enddo
  end subroutine add_to_lanczos_phonon



END MODULE ED_GF_PHONON













