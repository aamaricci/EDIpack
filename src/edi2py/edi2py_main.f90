!ED_MAIN:
subroutine init_solver_site_c(bath,dim_bath) bind(c, name='init_solver_site')
  integer(c_int64_t),dimension(1),intent(in)                 :: dim_bath
  real(c_double),dimension(dim_bath(1)),intent(inout)        :: bath
  call ed_init_solver(bath)
end subroutine init_solver_site_c
!
subroutine init_solver_ineq_c(bath,dim_bath) bind(c, name='init_solver_ineq')
  integer(c_int64_t),dimension(2),intent(in)                           :: dim_bath
  real(c_double),dimension(dim_bath(1),dim_bath(2)),intent(inout)      :: bath
  call ed_init_solver(bath)
end subroutine init_solver_ineq_c


subroutine solve_site_c(bath,dim_bath,hloc,dim_hloc,sflag) bind(c, name='solve_site')
  integer(c_int64_t),dimension(1),intent(in)                                             :: dim_bath
  integer(c_int64_t),dimension(4),intent(in)                                             :: dim_hloc
  integer(c_int),value                                                                   :: sflag
  real(c_double),dimension(dim_bath(1)),intent(in)                                       :: bath
  real(c_double),dimension(dim_hloc(1),dim_hloc(2),dim_hloc(3),dim_hloc(4)),intent(in)   :: hloc
  call assert_shape(hloc,[Nspin,Nspin,Norb,Norb],"solve","hloc")
  call ed_solve(bath,hloc,sflag=i2l(sflag))
end subroutine solve_site_c
!
subroutine solve_ineq_c(bath,dim_bath,hloc,dim_hloc,mpi_lanc) bind(c, name='solve_ineq')
  integer(c_int64_t),dimension(2),intent(in)                                                           :: dim_bath
  integer(c_int64_t),dimension(5),intent(in)                                                           :: dim_hloc
  integer(c_int),value                                                                                 :: mpi_lanc
  real(c_double),dimension(dim_bath(1),dim_bath(2)),intent(in)                                         :: bath
  real(c_double),dimension(dim_hloc(1),dim_hloc(2),dim_hloc(3),dim_hloc(4),dim_hloc(5)),intent(in)     :: hloc
  integer                                                                                              :: Nineq
  Nineq = size(bath,1)
  call assert_shape(Hloc,[Nineq,Nspin,Nspin,Norb,Norb],"solve","hloc")  
  call ed_solve(bath,hloc,i2l(mpi_lanc))
end subroutine solve_ineq_c
