!NORMAL, MATSUBARA GREEN'S FUNCTIONS
subroutine ed_get_gimp_matsubara_main(Gmats)
  complex(8),dimension(Nspin,Nspin,Norb,Norb,Lmats),intent(inout) :: Gmats
  Gmats(:,:,:,:,:) = impGmats(:,:,:,:,:)
end subroutine ed_get_gimp_matsubara_main

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


!NORMAL, REAL SELF-ENERGY
subroutine ed_get_gimp_realaxis_main(Greal)
  complex(8),dimension(Nspin,Nspin,Norb,Norb,Lreal),intent(inout) :: Greal
  Greal(:,:,:,:,:) = impGreal(:,:,:,:,:)
end subroutine ed_get_gimp_realaxis_main

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


!NORMAL, MATSUBARA SELF-ENEGRGY
subroutine ed_get_sigma_matsubara_main(Smats)
  complex(8),dimension(Nspin,Nspin,Norb,Norb,Lmats),intent(inout) :: Smats
  Smats(:,:,:,:,:) = impSmats(:,:,:,:,:)
end subroutine ed_get_sigma_matsubara_main

subroutine ed_get_sigma_matsubara_lattice(Smats,Nsites)
  integer                                                                :: Nsites
  complex(8),dimension(Nsites,Nspin,Nspin,Norb,Norb,Lmats),intent(inout) :: Smats
  Smats(1:Nsites,:,:,:,:,:) = Smatsii(1:Nsites,:,:,:,:,:)
end subroutine ed_get_sigma_matsubara_lattice

subroutine ed_get_sigma_matsubara_site(Smats,ilat)
  integer                                                         :: ilat
  complex(8),dimension(Nspin,Nspin,Norb,Norb,Lmats),intent(inout) :: Smats
  Smats(:,:,:,:,:) = Smatsii(ilat,:,:,:,:,:)
end subroutine ed_get_sigma_matsubara_site


!NORMAL, REAL SELF-ENERGY
subroutine ed_get_sigma_realaxis_main(Sreal)
  complex(8),dimension(Nspin,Nspin,Norb,Norb,Lreal),intent(inout) :: Sreal
  Sreal(:,:,:,:,:) = impSreal(:,:,:,:,:)
end subroutine ed_get_sigma_realaxis_main

subroutine ed_get_sigma_realaxis_lattice(Sreal,Nsites)
  integer                                                                :: Nsites
  complex(8),dimension(Nsites,Nspin,Nspin,Norb,Norb,Lreal),intent(inout) :: Sreal
  Sreal(1:Nsites,:,:,:,:,:) = Srealii(1:Nsites,:,:,:,:,:)
end subroutine ed_get_sigma_realaxis_lattice

subroutine ed_get_sigma_realaxis_site(Sreal,ilat)
  integer                                                         :: ilat
  complex(8),dimension(Nspin,Nspin,Norb,Norb,Lreal),intent(inout) :: Sreal
  Sreal(:,:,:,:,:) = Srealii(ilat,:,:,:,:,:)
end subroutine ed_get_sigma_realaxis_site
