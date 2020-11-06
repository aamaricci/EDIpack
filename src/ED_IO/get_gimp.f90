!NORMAL, MATSUBARA GREEN'S FUNCTIONS
subroutine ed_get_gimp_matsubara(Gmats)
  complex(8),dimension(Nspin,Nspin,Norb,Norb,Lmats),intent(inout) :: Gmats
  Gmats(:,:,:,:,:) = impGmats(:,:,:,:,:)
end subroutine ed_get_gimp_matsubara

subroutine ed_get_gimp_matsubara_lattice(Gmats,Nsites)
  integer                                                                :: Nsites
  complex(8),dimension(Nsites,Nspin,Nspin,Norb,Norb,Lmats),intent(inout) :: Gmats
  Gmats(1:Nsites,:,:,:,:,:) = Gmatsii(1:Nsites,:,:,:,:,:)
end subroutine ed_get_gimp_matsubara_lattice

subroutine ed_get_gimp_matsubara_site(Gmats,ilat)
  integer                                                         :: ilat
  complex(8),dimension(Nspin,Nspin,Norb,Norb,Lmats),intent(inout) :: Gmats
  Gmats(:,:,:,:,:) = Gmatsii(ilat,:,:,:,:,:)
end subroutine ed_get_gimp_matsubara_site





subroutine ed_get_gimp_realaxis(Greal)
  complex(8),dimension(Nspin,Nspin,Norb,Norb,Lreal),intent(inout) :: Greal
  Greal(:,:,:,:,:) = impGreal(:,:,:,:,:)
end subroutine ed_get_gimp_realaxis

subroutine ed_get_gimp_realaxis_lattice(Greal,Nsites)
  integer                                                                :: Nsites
  complex(8),dimension(Nsites,Nspin,Nspin,Norb,Norb,Lreal),intent(inout) :: Greal
  Greal(1:Nsites,:,:,:,:,:) = Grealii(1:Nsites,:,:,:,:,:)
end subroutine ed_get_gimp_realaxis_lattice

subroutine ed_get_gimp_realaxis_site(Greal,ilat)
  integer                                                         :: ilat
  complex(8),dimension(Nspin,Nspin,Norb,Norb,Lreal),intent(inout) :: Greal
  Greal(:,:,:,:,:) = Grealii(ilat,:,:,:,:,:)
end subroutine ed_get_gimp_realaxis_site
