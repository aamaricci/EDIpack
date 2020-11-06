subroutine ed_get_docc_1(docc) 
  real(8),dimension(Norb) :: docc
  docc = ed_docc
end subroutine ed_get_docc_1

subroutine ed_get_docc_2(docc,iorb) 
  real(8)   :: docc
  integer   :: iorb
  if(iorb>Norb)stop "ed_get_docc error: orbital index > N_orbital"
  docc = ed_docc(iorb)
end subroutine ed_get_docc_2

subroutine ed_get_docc_lattice_1(yii,Nlat) 
  integer                      :: Nlat
  real(8),dimension(Nlat,Norb) :: yii
  yii=0d0
  if(allocated(dii))then
     if(Nlat>size(dii,1)) stop "ed_get_docc error: required N_sites > evaluated N_sites"
     yii=dii
  endif
end subroutine ed_get_docc_lattice_1

subroutine ed_get_docc_lattice_2(yii,Nlat,iorb) 
  integer                 :: Nlat
  real(8),dimension(Nlat) :: yii
  integer                 :: iorb
  if(iorb>Norb)stop "ed_get_docc error: orbital index > N_orbital"
  yii=0d0
  if(allocated(dii))then
     if(Nlat>size(dii,1)) stop "ed_get_docc error: required N_sites > evaluated N_sites"
     yii=dii(:,iorb)
  endif
end subroutine ed_get_docc_lattice_2
