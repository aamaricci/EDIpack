!USE ED_BATH_FIT:
subroutine chi2_fitgf_site_c(func,dim_func,bath,dim_bath,ispin,iorb) bind(c, name='chi2_fitgf_site')
  integer(c_int64_t)                                                                                             :: dim_func(5),dim_bath(1)
  real(c_double),dimension(dim_bath(1)),intent(inout)                                                            :: bath
  complex(c_double_complex),dimension(dim_func(1),dim_func(2),dim_func(3),dim_func(4),dim_func(5)),intent(in)    :: func
  integer(c_int),value                                                                                           :: ispin
  integer(c_int),value                                                                                           :: iorb     
  call assert_shape(func,[Nspin,Nspin,Norb,Norb,Lmats],"chi2_fitgf","func")
  if(iorb>0)then
     call ed_chi2_fitgf(func,bath,ispin,iorb)
  else
     call ed_chi2_fitgf(func,bath,ispin)
  endif
end subroutine chi2_fitgf_site_c

subroutine chi2_fitgf_ineq_c(func,dim_func,bath,dim_bath,hloc,dim_hloc,ispin) bind(c, name='chi2_fitgf_ineq')
  integer(c_int64_t)                                                                                                          :: dim_func(6),dim_bath(2),dim_hloc(5)
  real(c_double),dimension(dim_bath(1),dim_bath(2)),intent(inout)                                                             :: bath
  complex(c_double_complex),dimension(dim_func(1),dim_func(2),dim_func(3),dim_func(4),dim_func(5),dim_func(6)),intent(inout)  :: func
  real(c_double),dimension(dim_hloc(1),dim_hloc(2),dim_hloc(3),dim_hloc(4),dim_hloc(5)),intent(in)                            :: hloc
  integer(c_int),intent(in)                                                                                                   :: ispin
  integer(c_int)                                                                                                              :: Nineq
  Nineq = size(bath,1)
  call assert_shape(func,[Nineq,Nspin,Nspin,Norb,Norb,Lmats],"chi2_fitgf","func")
  call assert_shape(Hloc,[Nineq,Nspin,Nspin,Norb,Norb],"chi2_fitgf","hloc")
  call ed_chi2_fitgf(func,bath,hloc,ispin)
end subroutine chi2_fitgf_ineq_c
