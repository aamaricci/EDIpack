MODULE ED_GREENS_FUNCTIONS
  USE SF_CONSTANTS, only:one,xi,zero,pi
  USE SF_TIMER  
  USE SF_IOTOOLS, only: str,free_unit,reg,free_units,txtfy
  USE SF_LINALG,  only: inv,eigh,eye
  USE ED_INPUT_VARS
  USE ED_VARS_GLOBAL
  USE ED_IO                     !< this contains the routine to print GF,Sigma and G0
  USE ED_EIGENSPACE
  USE ED_BATH
  USE ED_BATH_FUNCTIONS
  USE ED_SETUP
  USE ED_HAMILTONIAN
  USE ED_AUX_FUNX
  !
  USE ED_GF_ELECTRON
  USE ED_GF_PHONON
  !
  implicit none
  private 

  public :: buildGf_impurity

contains


  !+------------------------------------------------------------------+
  ! GF CALCULATIONS
  !+------------------------------------------------------------------+
  subroutine buildGF_impurity()
    !
    call allocate_grids
    !
    impGmats=zero
    impGreal=zero
    !
    impSmats = zero
    impSreal = zero
    !
    impG0mats=zero
    impG0real=zero
    !
    impDmats_ph=zero
    impDreal_ph=zero
    !
    write(LOGfile,"(A)")"Get impurity Greens functions:"
    call build_gf_normal()
    if(DimPh>1) call build_gf_phonon()
    call build_sigma_normal()
    !
    if(MPIMASTER)then
       if(ed_print_Sigma)call ed_print_impSigma()
       if(ed_print_G) then
          call ed_print_impG()
          if(DimPh>1)call ed_print_impD()
       endif
       if(ed_print_G0)call ed_print_impG0()
    endif
    !
    call deallocate_grids
    !
  end subroutine buildgf_impurity



end MODULE ED_GREENS_FUNCTIONS
