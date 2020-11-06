!NORMAL, MATSUBARA GREEN'S FUNCTIONS
subroutine ed_get_g0imp_matsubara(Gmats)
  complex(8),dimension(Nspin,Nspin,Norb,Norb,Lmats),intent(inout) :: Gmats
  Gmats(:,:,:,:,:) = impG0mats(:,:,:,:,:)
end subroutine ed_get_g0imp_matsubara

subroutine ed_get_g0imp_matsubara_lattice(Gmats,Nsites)
  integer                                                                :: Nsites
  complex(8),dimension(Nsites,Nspin,Nspin,Norb,Norb,Lmats),intent(inout) :: Gmats
  Gmats(1:Nsites,:,:,:,:,:) = G0matsii(1:Nsites,:,:,:,:,:)
end subroutine ed_get_g0imp_matsubara_lattice

subroutine ed_get_g0imp_matsubara_site(Gmats,ilat)
  integer                                                         :: ilat
  complex(8),dimension(Nspin,Nspin,Norb,Norb,Lmats),intent(inout) :: Gmats
  Gmats(:,:,:,:,:) = G0matsii(ilat,:,:,:,:,:)
end subroutine ed_get_g0imp_matsubara_site



!NORMAL, REAL GREEN'S FUNCTION
subroutine ed_get_g0imp_realaxis(Greal)
  complex(8),dimension(Nspin,Nspin,Norb,Norb,Lreal),intent(inout) :: Greal
  Greal(:,:,:,:,:) = impG0real(:,:,:,:,:)
end subroutine ed_get_g0imp_realaxis

subroutine ed_get_g0imp_realaxis_lattice(Greal,Nsites)
  integer                                                                :: Nsites
  complex(8),dimension(Nsites,Nspin,Nspin,Norb,Norb,Lreal),intent(inout) :: Greal
  Greal(1:Nsites,:,:,:,:,:) = G0realii(1:Nsites,:,:,:,:,:)
end subroutine ed_get_g0imp_realaxis_lattice

subroutine ed_get_g0imp_realaxis_site(Greal,ilat)
  integer                                                         :: ilat
  complex(8),dimension(Nspin,Nspin,Norb,Norb,Lreal),intent(inout) :: Greal
  Greal(:,:,:,:,:) = G0realii(ilat,:,:,:,:,:)
end subroutine ed_get_g0imp_realaxis_site
