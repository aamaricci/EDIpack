!NORMAL, MATSUBARA SELF-ENEGRGY
subroutine ed_get_sigma_matsubara(Smats)
  complex(8),dimension(Nspin,Nspin,Norb,Norb,Lmats),intent(inout) :: Smats
  Smats(:,:,:,:,:) = impSmats(:,:,:,:,:)
end subroutine ed_get_sigma_matsubara

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
subroutine ed_get_sigma_realaxis(Sreal)
  complex(8),dimension(Nspin,Nspin,Norb,Norb,Lreal),intent(inout) :: Sreal
  Sreal(:,:,:,:,:) = impSreal(:,:,:,:,:)
end subroutine ed_get_sigma_realaxis

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
