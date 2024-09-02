!ED_MAIN:
subroutine init_solver_site_c(bath,dim_bath) bind(c, name='init_solver_site')
  integer(c_int),value                                :: dim_bath
  real(c_double),dimension(dim_bath),intent(inout)    :: bath
  call ed_init_solver(bath)
end subroutine init_solver_site_c
!
subroutine init_solver_ineq_c(bath,d1,d2) bind(c, name='init_solver_ineq')
  integer(c_int),value                           :: d1,d2
  real(c_double),dimension(d1,d2),intent(inout)  :: bath
  call ed_init_solver(bath)
end subroutine init_solver_ineq_c


subroutine solve_site_c(bath,dim_bath_1,hloc,dim_hloc_1,dim_hloc_2,dim_hloc_3,dim_hloc_4,sflag) bind(c, name='solve_site')
  integer(c_int),value                                                               :: dim_bath_1,dim_hloc_1,dim_hloc_2,dim_hloc_3,dim_hloc_4
  integer(c_int),value                                                               :: sflag
  real(c_double),dimension(dim_bath_1),intent(in)                                    :: bath
  real(c_double),dimension(dim_hloc_1,dim_hloc_2,dim_hloc_3,dim_hloc_4),intent(in)   :: hloc
  call assert_shape(hloc,[Nspin,Nspin,Norb,Norb],"solve","hloc")
  call ed_solve(bath,hloc,sflag=i2l(sflag))
end subroutine solve_site_c
!
subroutine solve_ineq_c(bath,dim_bath_1,dim_bath_2,hloc,dim_hloc_1,dim_hloc_2,dim_hloc_3,dim_hloc_4,dim_hloc_5,mpi_lanc) bind(c, name='solve_ineq')
  integer(c_int),value                                                                        :: dim_bath_1,dim_bath_2,dim_hloc_1,dim_hloc_2,dim_hloc_3,dim_hloc_4,dim_hloc_5
  integer(c_int),value                                                                        :: mpi_lanc
  real(c_double),dimension(dim_bath_1,dim_bath_2),intent(in)                                  :: bath
  real(c_double),dimension(dim_hloc_1,dim_hloc_2,dim_hloc_3,dim_hloc_4,dim_hloc_5),intent(in) :: hloc
  integer                                                                                     :: Nineq
  Nineq = size(bath,1)
  call assert_shape(Hloc,[Nineq,Nspin,Nspin,Norb,Norb],"solve","hloc")  
  call ed_solve(bath,hloc,i2l(mpi_lanc))
end subroutine solve_ineq_c
