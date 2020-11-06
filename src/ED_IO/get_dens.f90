subroutine ed_get_dens_1(dens)
  real(8),dimension(Norb) :: dens
  dens = ed_dens
end subroutine ed_get_dens_1

subroutine ed_get_dens_2(dens,iorb)
  real(8)   :: dens
  integer   :: iorb
  if(iorb>Norb)stop "ed_get_dens error: orbital index > N_orbital"
  dens = ed_dens(iorb)
end subroutine ed_get_dens_2

subroutine ed_get_dens_lattice_1(yii,Nlat)
  integer                      :: Nlat
  real(8),dimension(Nlat,Norb) :: yii
  yii=0d0    
  if(allocated(nii))then
     if(Nlat>size(nii,1)) stop "ed_get_dens error: required N_sites > evaluated N_sites"
     yii=nii
  end if
end subroutine ed_get_dens_lattice_1

subroutine ed_get_dens_lattice_2(yii,Nlat,iorb)
  integer                 :: Nlat
  real(8),dimension(Nlat) :: yii
  integer                 :: iorb
  if(iorb>Norb)stop "ed_get_dens error: orbital index > N_orbital"
  yii=0d0
  if(allocated(nii))then
     if(Nlat>size(nii,1)) stop "ed_get_dens error: required N_sites > evaluated N_sites"
     yii=nii(:,iorb)
  endif
end subroutine ed_get_dens_lattice_2
