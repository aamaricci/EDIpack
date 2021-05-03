!NORMAL, MATSUBARA GREEN'S FUNCTIONS
subroutine ed_get_gimp_matsubara_main(Gmats)
  complex(8),dimension(Nspin,Nspin,Norb,Norb,Lmats),intent(inout) :: Gmats
  Gmats(:,:,:,:,:) = impGmats(:,:,:,:,:)
end subroutine ed_get_gimp_matsubara_main

subroutine ed_get_gimp_matsubara_lattice(Gmats,Nsites)
  integer                                                                :: Nsites
  complex(8),dimension(Nsites,Nspin,Nspin,Norb,Norb,Lmats),intent(inout) :: Gmats
  Gmats(1:Nsites,:,:,:,:,:) = Gmats_ineq(1:Nsites,:,:,:,:,:)
end subroutine ed_get_gimp_matsubara_lattice



!NORMAL, REAL SELF-ENERGY
subroutine ed_get_gimp_realaxis_main(Greal)
  complex(8),dimension(Nspin,Nspin,Norb,Norb,Lreal),intent(inout) :: Greal
  Greal(:,:,:,:,:) = impGreal(:,:,:,:,:)
end subroutine ed_get_gimp_realaxis_main

subroutine ed_get_gimp_realaxis_lattice(Greal,Nsites)
  integer                                                                :: Nsites
  complex(8),dimension(Nsites,Nspin,Nspin,Norb,Norb,Lreal),intent(inout) :: Greal
  Greal(1:Nsites,:,:,:,:,:) = Greal_ineq(1:Nsites,:,:,:,:,:)
end subroutine ed_get_gimp_realaxis_lattice


!NORMAL, MATSUBARA SELF-ENEGRGY
subroutine ed_get_sigma_matsubara_main(Smats)
  complex(8),dimension(Nspin,Nspin,Norb,Norb,Lmats),intent(inout) :: Smats
  Smats(:,:,:,:,:) = impSmats(:,:,:,:,:)
end subroutine ed_get_sigma_matsubara_main

subroutine ed_get_sigma_matsubara_lattice(Smats,Nsites)
  integer                                                                :: Nsites
  complex(8),dimension(Nsites,Nspin,Nspin,Norb,Norb,Lmats),intent(inout) :: Smats
  Smats(1:Nsites,:,:,:,:,:) = Smats_ineq(1:Nsites,:,:,:,:,:)
end subroutine ed_get_sigma_matsubara_lattice


!NORMAL, REAL SELF-ENERGY
subroutine ed_get_sigma_realaxis_main(Sreal)
  complex(8),dimension(Nspin,Nspin,Norb,Norb,Lreal),intent(inout) :: Sreal
  Sreal(:,:,:,:,:) = impSreal(:,:,:,:,:)
end subroutine ed_get_sigma_realaxis_main

subroutine ed_get_sigma_realaxis_lattice(Sreal,Nsites)
  integer                                                                :: Nsites
  complex(8),dimension(Nsites,Nspin,Nspin,Norb,Norb,Lreal),intent(inout) :: Sreal
  Sreal(1:Nsites,:,:,:,:,:) = Sreal_ineq(1:Nsites,:,:,:,:,:)
end subroutine ed_get_sigma_realaxis_lattice
