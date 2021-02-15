MODULE ED_CHI_FUNCTIONS
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
  !
  USE ED_CHI_SPIN
  USE ED_CHI_DENS
  USE ED_CHI_PAIR
  USE ED_CHI_EXCT
  !
  implicit none
  private 

  public :: buildChi_impurity

contains


  !+------------------------------------------------------------------+
  ! SUSCEPTIBILITY CALCULATIONS
  !+------------------------------------------------------------------+
  subroutine buildChi_impurity()
    !
    call allocate_grids
    !
    !
    !BUILD SPIN SUSCEPTIBILITY
    spinChi_tau=zero
    spinChi_w=zero
    spinChi_iv=zero
    if(chispin_flag)call build_chi_spin()
    !
    !
    !BUILD CHARGE SUSCEPTIBILITY
    densChi_tau=zero
    densChi_w=zero
    densChi_iv=zero
    if(chidens_flag)call build_chi_dens()
    !
    !
    !BUILD PAIR SUSCEPTIBILITY
    pairChi_tau=zero
    pairChi_w=zero
    pairChi_iv=zero
    if(chipair_flag)call build_chi_pair()
    !
    !BUILD EXCITON SUSCEPTIBILITY
    exctChi_tau=zero
    exctChi_w=zero
    exctChi_iv=zero
    if(chiexct_flag)call build_chi_exct()
    !
    !
    !PRINTING:
    if(MPIMASTER.AND.(any([chispin_flag,chidens_flag,chipair_flag,chiexct_flag])))call ed_print_impChi()
    !
    call deallocate_grids
    !
  end subroutine buildChi_impurity




end MODULE ED_CHI_FUNCTIONS
