!ED_IO:
subroutine get_sigma_matsubara_site_c(sigma,d) bind(c, name='get_sigma_matsubara_site')
  integer(c_int64_t)                                                          :: d(5)
  complex(c_double_complex),dimension(d(1),d(2),d(3),d(4),d(5)),intent(inout) :: sigma
  call assert_shape(sigma,[Nspin,Nspin,Norb,Norb,Lmats],"get_sigma_matsubara","Sigma")
  call ed_get_sigma_matsubara(sigma)
end subroutine get_sigma_matsubara_site_c
!
subroutine get_sigma_matsubara_ineq_c(sigma,d) bind(c, name='get_sigma_matsubara_ineq')
  integer(c_int64_t)                                                               :: d(6)
  complex(c_double_complex),dimension(d(1),d(2),d(3),d(4),d(5),d(6)),intent(inout) :: sigma
  integer(c_int)                                                                   :: Nsites
  Nsites=size(sigma,1)
  call assert_shape(sigma,[Nsites,Nspin,Nspin,Norb,Norb,Lmats],"get_sigma_matsubara","Sigma")
  call ed_get_sigma_matsubara(sigma,Nsites)
end subroutine get_sigma_matsubara_ineq_c



subroutine get_sigma_realaxis_site_c(sigma,d)  bind(c, name='get_sigma_realaxis_site')
  integer(c_int64_t)                                                          :: d(5)
  complex(c_double_complex),dimension(d(1),d(2),d(3),d(4),d(5)),intent(inout) :: sigma
  call assert_shape(sigma,[Nspin,Nspin,Norb,Norb,Lreal],"get_sigma_realaxis","sigma")
  call ed_get_sigma_realaxis(sigma)
end subroutine get_sigma_realaxis_site_c
!
subroutine get_sigma_realaxis_ineq_c(sigma,d)  bind(c, name='get_sigma_realaxis_ineq')
  integer(c_int64_t)                                                               :: d(6)
  complex(c_double_complex),dimension(d(1),d(2),d(3),d(4),d(5),d(6)),intent(inout) :: sigma
  integer(c_int)                                                                   :: Nsites
  Nsites=size(sigma,1)
  call assert_shape(sigma,[Nsites,Nspin,Nspin,Norb,Norb,Lreal],"get_sigma_realaxis","sigma")
  call ed_get_sigma_realaxis(sigma,Nsites)
end subroutine get_sigma_realaxis_ineq_c



subroutine get_gimp_matsubara_site_c(gimp,d) bind(c, name='get_gimp_matsubara_site')
  integer(c_int64_t)                                                          :: d(5)
  complex(c_double_complex),dimension(d(1),d(2),d(3),d(4),d(5)),intent(inout) :: gimp
  call assert_shape(gimp,[Nspin,Nspin,Norb,Norb,Lmats],"get_gimp_matsubara","gimp")
  call ed_get_gimp_matsubara(gimp)
end subroutine get_gimp_matsubara_site_c
!
subroutine get_gimp_matsubara_ineq_c(gimp,d) bind(c, name='get_gimp_matsubara_ineq')
  integer(c_int64_t)                                                               :: d(6)
  complex(c_double_complex),dimension(d(1),d(2),d(3),d(4),d(5),d(6)),intent(inout) :: gimp
  integer(c_int)                                                                   :: Nsites
  Nsites=size(gimp,1)
  call assert_shape(gimp,[Nsites,Nspin,Nspin,Norb,Norb,Lmats],"get_gimp_matsubara","gimp")
  call ed_get_gimp_matsubara(gimp,Nsites)
end subroutine get_gimp_matsubara_ineq_c



subroutine get_gimp_realaxis_site_c(gimp,d)  bind(c, name='get_gimp_realaxis_site')
  integer(c_int64_t)                                                          :: d(5)
  complex(c_double_complex),dimension(d(1),d(2),d(3),d(4),d(5)),intent(inout) :: gimp
  call assert_shape(gimp,[Nspin,Nspin,Norb,Norb,Lreal],"get_gimp_realaxis","gimp")
  call ed_get_gimp_realaxis(gimp)
end subroutine get_gimp_realaxis_site_c
!
subroutine get_gimp_realaxis_ineq_c(gimp,d)  bind(c, name='get_gimp_realaxis_ineq')
  integer(c_int64_t)                                                               :: d(6)
  complex(c_double_complex),dimension(d(1),d(2),d(3),d(4),d(5),d(6)),intent(inout) :: gimp
  integer(c_int)                                                                   :: Nsites
  Nsites=size(gimp,1)
  call assert_shape(gimp,[Nsites,Nspin,Nspin,Norb,Norb,Lreal],"get_gimp_realaxis","gimp")
  call ed_get_gimp_realaxis(gimp,Nsites)
end subroutine get_gimp_realaxis_ineq_c



subroutine get_dens_site_c(arg,arg_dim1) bind(c, name='get_dens_site')
  integer(c_int),value                    :: arg_dim1
  real(c_double),dimension(arg_dim1)      :: arg
  call assert_shape(arg,[Norb],"get_dens_site","arg")
  call ed_get_dens(arg)
end subroutine get_dens_site_c
!
subroutine get_dens_ineq_c(arg,arg_dim1,arg_dim2,nlat) bind(c, name='get_dens_ineq')
  integer(c_int),value                          :: arg_dim1,arg_dim2,nlat
  real(c_double),dimension(arg_dim1,arg_dim2)   :: arg
  call assert_shape(arg,[Nlat,Norb],"get_dens_ineq","arg")  
  call ed_get_dens(arg,nlat)
end subroutine get_dens_ineq_c



subroutine get_mag_site_c(arg,arg_dim1) bind(c, name='get_mag_site')
  integer(c_int),value               :: arg_dim1
  real(c_double),dimension(arg_dim1) :: arg
  call assert_shape(arg,[Norb],"get_mag_site","arg")  
  call ed_get_mag(arg)
end subroutine get_mag_site_c
!
subroutine get_mag_ineq_c(arg,arg_dim1,arg_dim2,nlat) bind(c, name='get_mag_ineq')
  integer(c_int),value                        :: arg_dim1,arg_dim2,nlat
  real(c_double),dimension(arg_dim1,arg_dim2) :: arg
  call assert_shape(arg,[Nlat,Norb],"get_mag_ineq","arg")  
  call ed_get_mag(arg,nlat)
end subroutine get_mag_ineq_c


subroutine get_docc_site_c(arg,arg_dim1) bind(c, name='get_docc_site')
  integer(c_int),value                    :: arg_dim1
  real(c_double),dimension(arg_dim1)      :: arg
  call assert_shape(arg,[Norb],"get_doc_site","arg")  
  call ed_get_docc(arg)
end subroutine get_docc_site_c
!
subroutine get_docc_ineq_c(arg,arg_dim1,arg_dim2,nlat) bind(c, name='get_docc_ineq')
  integer(c_int),value                        :: arg_dim1,arg_dim2,nlat
  real(c_double),dimension(arg_dim1,arg_dim2) :: arg
  call assert_shape(arg,[Nlat,Norb],"get_docc_ineq","arg")  
  call ed_get_docc(arg,nlat)
end subroutine get_docc_ineq_c



subroutine get_eimp_site_c(arg) bind(c, name='get_eimp_site')
  real(c_double),dimension(4) :: arg
  call assert_shape(arg,[4],"get_eimp_site","arg")
  call ed_get_eimp(arg)
end subroutine get_eimp_site_c
!
subroutine get_eimp_ineq_c(arg,nlat) bind(c, name='get_eimp_ineq')
  integer(c_int),value                        :: nlat
  real(c_double),dimension(Nlat,4) :: arg
  call assert_shape(arg,[Nlat,4],"get_eimp_ineq","arg")  
  call ed_get_eimp(arg,nlat)
end subroutine get_eimp_ineq_c



subroutine get_doubles_site_c(arg) bind(c, name='get_doubles_site')
  real(c_double),dimension(4) :: arg
  call assert_shape(arg,[4],"get_doubles_site","arg")
  call ed_get_doubles(arg)
end subroutine get_doubles_site_c
!
subroutine get_doubles_ineq_c(arg,nlat) bind(c, name='get_doubles_ineq')
  integer(c_int),value             :: nlat
  real(c_double),dimension(Nlat,4) :: arg
  call assert_shape(arg,[Nlat,4],"get_doubles_ineq","arg")  
  call ed_get_doubles(arg,nlat)
end subroutine get_doubles_ineq_c







