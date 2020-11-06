subroutine ed_get_mag_1(mag) 
  real(8),dimension(Norb) :: mag
  mag = (ed_dens_up-ed_dens_dw)
end subroutine ed_get_mag_1

subroutine ed_get_mag_2(mag,iorb) 
  real(8)   :: mag
  integer   :: iorb
  if(iorb>Norb)stop "ed_get_dens_up error: orbital index > N_orbital"
  mag = ed_dens_up(iorb)-ed_dens_dw(iorb)
end subroutine ed_get_mag_2

subroutine ed_get_mag_lattice_1(yii,Nlat)
  integer                      :: Nlat
  real(8),dimension(Nlat,Norb) :: yii
  yii=0d0
  if(allocated(mii))then
     if(Nlat>size(mii,1)) stop "ed_get_mag error: required N_sites > evaluated N_sites"
     yii=mii
  endif
end subroutine ed_get_mag_lattice_1

subroutine ed_get_mag_lattice_2(yii,Nlat,iorb) 
  integer                 :: Nlat
  real(8),dimension(Nlat) :: yii
  integer                 :: iorb
  if(iorb>Norb)stop "ed_get_mag error: orbital index > N_orbital"
  yii=0d0
  if(allocated(mii))then
     if(Nlat>size(mii,1)) stop "ed_get_mag error: required N_sites > evaluated N_sites"
     yii=mii(:,iorb)
  endif
end subroutine ed_get_mag_lattice_2
