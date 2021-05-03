subroutine ed_get_dens_main(dens)
  real(8),dimension(Norb) :: dens
  dens = ed_dens
end subroutine ed_get_dens_main

subroutine ed_get_mag_main(mag) 
  real(8),dimension(Norb) :: mag
  mag = (ed_dens_up-ed_dens_dw)
end subroutine ed_get_mag_main

subroutine ed_get_docc_main(docc) 
  real(8),dimension(Norb) :: docc
  docc = ed_docc
end subroutine ed_get_docc_main

subroutine ed_get_doubles_main(docc)
  real(8),dimension(4) :: docc
  docc = [ed_Dust,ed_Dund,ed_Dse,ed_Dph]
end subroutine ed_get_doubles_main


subroutine ed_get_eimp_main(eimp)
  real(8),dimension(4) :: eimp
  eimp = [ed_Epot,ed_Eint,ed_Ehartree,ed_Eknot]
end subroutine ed_get_eimp_main







subroutine ed_get_dens_lattice(y_ineq,Nlat)
  integer                      :: Nlat
  real(8),dimension(Nlat,Norb) :: y_ineq
  y_ineq=0d0    
  if(allocated(dens_ineq))then
     if(Nlat>size(dens_ineq,1)) stop "ed_get_dens error: required N_sites > evaluated N_sites"
     y_ineq=dens_ineq
  end if
end subroutine ed_get_dens_lattice

subroutine ed_get_mag_lattice(y_ineq,Nlat)
  integer                      :: Nlat
  real(8),dimension(Nlat,Norb) :: y_ineq
  y_ineq=0d0
  if(allocated(mag_ineq))then
     if(Nlat>size(mag_ineq,1)) stop "ed_get_mag error: required N_sites > evaluated N_sites"
     y_ineq=mag_ineq
  endif
end subroutine ed_get_mag_lattice

subroutine ed_get_docc_lattice(y_ineq,Nlat) 
  integer                      :: Nlat
  real(8),dimension(Nlat,Norb) :: y_ineq
  y_ineq=0d0
  if(allocated(docc_ineq))then
     if(Nlat>size(docc_ineq,1)) stop "ed_get_docc error: required N_sites > evaluated N_sites"
     y_ineq=docc_ineq
  endif
end subroutine ed_get_docc_lattice

subroutine ed_get_doubles_lattice(y_ineq,Nlat)
  integer                      :: Nlat
  real(8),dimension(Nlat,4)    :: y_ineq
  y_ineq=0d0
  if(allocated(dd_ineq))then
     if(Nlat>size(dd_ineq,1)) stop "ed_get_doubles error: required N_sites > evaluated N_sites"
     y_ineq=dd_ineq(:,:)
  endif
end subroutine ed_get_doubles_lattice

subroutine ed_get_eimp_lattice(y_ineq,Nlat)
  integer                      :: Nlat
  real(8),dimension(Nlat,4)    :: y_ineq
  y_ineq=0d0
  if(allocated(e_ineq))then
     if(Nlat>size(e_ineq,1)) stop "ed_get_eimp error: required N_sites > evaluated N_sites"
     y_ineq=e_ineq
  endif
end subroutine ed_get_eimp_lattice


