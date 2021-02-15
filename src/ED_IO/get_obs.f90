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

subroutine ed_get_dust_main(docc)
  real(8) :: docc
  docc = ed_Dust
end subroutine ed_get_dust_main

subroutine ed_get_dund_main(docc)
  real(8) :: docc
  docc = ed_Dund
end subroutine ed_get_dund_main

subroutine ed_get_dse_main(docc)
  real(8) :: docc
  docc = ed_Dse
end subroutine ed_get_dse_main

subroutine ed_get_dph_main(docc)
  real(8) :: docc
  docc = ed_Dph
end subroutine ed_get_dph_main


subroutine ed_get_eimp_main(eimp)
  real(8),dimension(4) :: eimp
  eimp = [ed_Epot,ed_Eint,ed_Ehartree,ed_Eknot]
end subroutine ed_get_eimp_main

subroutine ed_get_epot_main(eimp)
  real(8) :: eimp
  eimp = ed_Epot
end subroutine ed_get_epot_main

subroutine ed_get_eint_main(eimp)
  real(8) :: eimp
  eimp = ed_Eint
end subroutine ed_get_eint_main

subroutine ed_get_ehartree_main(eimp)
  real(8) :: eimp
  eimp = ed_Ehartree
end subroutine ed_get_ehartree_main

subroutine ed_get_eknot_main(eimp)
  real(8) :: eimp
  eimp = ed_Eknot
end subroutine ed_get_eknot_main






subroutine ed_get_dens_lattice(yii,Nlat)
  integer                      :: Nlat
  real(8),dimension(Nlat,Norb) :: yii
  yii=0d0    
  if(allocated(nii))then
     if(Nlat>size(nii,1)) stop "ed_get_dens error: required N_sites > evaluated N_sites"
     yii=nii
  end if
end subroutine ed_get_dens_lattice

subroutine ed_get_mag_lattice(yii,Nlat)
  integer                      :: Nlat
  real(8),dimension(Nlat,Norb) :: yii
  yii=0d0
  if(allocated(mii))then
     if(Nlat>size(mii,1)) stop "ed_get_mag error: required N_sites > evaluated N_sites"
     yii=mii
  endif
end subroutine ed_get_mag_lattice

subroutine ed_get_docc_lattice(yii,Nlat) 
  integer                      :: Nlat
  real(8),dimension(Nlat,Norb) :: yii
  yii=0d0
  if(allocated(dii))then
     if(Nlat>size(dii,1)) stop "ed_get_docc error: required N_sites > evaluated N_sites"
     yii=dii
  endif
end subroutine ed_get_docc_lattice

subroutine ed_get_doubles_lattice(yii,Nlat)
  integer                      :: Nlat
  real(8),dimension(Nlat,4)    :: yii
  yii=0d0
  if(allocated(ddii))then
     if(Nlat>size(ddii,1)) stop "ed_get_doubles error: required N_sites > evaluated N_sites"
     yii=ddii(:,:)
  endif
end subroutine ed_get_doubles_lattice

subroutine ed_get_dust_lattice(yii,Nlat)
  integer                 :: Nlat
  real(8),dimension(Nlat) :: yii
  yii=0d0
  if(allocated(ddii))then
     if(Nlat>size(ddii,1)) stop "ed_get_dust error: required N_sites > evaluated N_sites"
     yii=ddii(:,1)
  endif
end subroutine ed_get_dust_lattice

subroutine ed_get_dund_lattice(yii,Nlat)
  integer                 :: Nlat
  real(8),dimension(Nlat) :: yii
  yii=0d0
  if(allocated(ddii))then
     if(Nlat>size(ddii,1)) stop "ed_get_dund error: required N_sites > evaluated N_sites"
     yii=ddii(:,2)
  endif
end subroutine ed_get_dund_lattice

subroutine ed_get_dse_lattice(yii,Nlat)
  integer                 :: Nlat
  real(8),dimension(Nlat) :: yii
  yii=0d0
  if(allocated(ddii))then
     if(Nlat>size(ddii,1)) stop "ed_get_dse error: required N_sites > evaluated N_sites"
     yii=ddii(:,3)
  endif
end subroutine ed_get_dse_lattice

subroutine ed_get_dph_lattice(yii,Nlat)
  integer                 :: Nlat
  real(8),dimension(Nlat) :: yii
  yii=0d0
  if(allocated(ddii))then
     if(Nlat>size(ddii,1)) stop "ed_get_dph error: required N_sites > evaluated N_sites"
     yii=ddii(:,4)
  endif
end subroutine ed_get_dph_lattice

subroutine ed_get_eimp_lattice(yii,Nlat)
  integer                      :: Nlat
  real(8),dimension(Nlat,4)    :: yii
  yii=0d0
  if(allocated(eii))then
     if(Nlat>size(eii,1)) stop "ed_get_eimp error: required N_sites > evaluated N_sites"
     yii=eii
  endif
end subroutine ed_get_eimp_lattice

subroutine ed_get_epot_lattice(yii,Nlat)
  integer                 :: Nlat
  real(8),dimension(Nlat) :: yii
  yii=0d0
  if(allocated(eii))then
     if(Nlat>size(eii,1)) stop "ed_get_epot error: required N_sites > evaluated N_sites"
     yii=eii(:,1)
  endif
end subroutine ed_get_epot_lattice

subroutine ed_get_eint_lattice(yii,Nlat)
  integer                 :: Nlat
  real(8),dimension(Nlat) :: yii
  yii=0d0
  if(allocated(eii))then
     if(Nlat>size(eii,1)) stop "ed_get_eint error: required N_sites > evaluated N_sites"
     yii=eii(:,2)
  endif
end subroutine ed_get_eint_lattice

subroutine ed_get_ehartree_lattice(yii,Nlat)
  integer                 :: Nlat
  real(8),dimension(Nlat) :: yii
  yii=0d0
  if(allocated(eii))then
     if(Nlat>size(eii,1)) stop "ed_get_ehartree error: required N_sites > evaluated N_sites"
     yii=eii(:,3)
  endif
end subroutine ed_get_ehartree_lattice

subroutine ed_get_eknot_lattice(yii,Nlat)
  integer                 :: Nlat
  real(8),dimension(Nlat) :: yii
  yii=0d0
  if(allocated(eii))then
     if(Nlat>size(eii,1)) stop "ed_get_knot error: required N_sites > evaluated N_sites"
     yii=eii(:,4)
  endif
end subroutine ed_get_eknot_lattice


