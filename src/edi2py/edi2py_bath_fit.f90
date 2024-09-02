!USE ED_BATH_FIT:
subroutine chi2_fitgf_site_c(func,dim_func_1,dim_func_2,dim_func_3,dim_func_4,dim_func_5,bath,dim_bath,ispin,iorb) bind(c, name='chi2_fitgf_site')
  integer(c_int),value                                                                                       :: dim_func_1,dim_func_2,dim_func_3,dim_func_4,dim_func_5,dim_bath
  real(c_double),dimension(dim_bath),intent(inout)                                                           :: bath
  complex(c_double_complex),dimension(dim_func_1,dim_func_2,dim_func_3,dim_func_4,dim_func_5),intent(in)     :: func
  integer(c_int),value                                                                                       :: ispin
  integer(c_int),value                                                                                       :: iorb     
  call assert_shape(func,[Nspin,Nspin,Norb,Norb,Lmats],"chi2_fitgf","func")
  if(iorb>0)then
     call ed_chi2_fitgf(func,bath,ispin,iorb)
  else
     call ed_chi2_fitgf(func,bath,ispin)
  endif
end subroutine chi2_fitgf_site_c

subroutine chi2_fitgf_ineq_c(func,dim_func_1,dim_func_2,dim_func_3,dim_func_4,dim_func_5,dim_func_6,&
                             bath,dim_bath_1,dim_bath_2,&
                             hloc,dim_hloc_1,dim_hloc_2,dim_hloc_3,dim_hloc_4,dim_hloc_5,&
                             ispin) bind(c, name='chi2_fitgf_ineq')
  
  integer(c_int),value                                                          ::    dim_func_1,dim_func_2,dim_func_3,dim_func_4,dim_func_5,dim_func_6
  integer(c_int),value                                                          ::    dim_hloc_1,dim_hloc_2,dim_hloc_3,dim_hloc_4,dim_hloc_5
  integer(c_int),value                                                          ::    dim_bath_1,dim_bath_2
  real(c_double),dimension(dim_bath_1,dim_bath_2),intent(inout)                 :: bath
  complex(c_double_complex),dimension(dim_func_1,dim_func_2,dim_func_3,dim_func_4,dim_func_5,dim_func_6),intent(inout) :: func
  real(c_double),dimension(dim_hloc_1,dim_hloc_2,dim_hloc_3,dim_hloc_4,dim_hloc_5),intent(in)         :: hloc
  integer(c_int),intent(in)                       :: ispin
  integer(c_int)                                  :: Nineq
  Nineq = size(bath,1)
  call assert_shape(func,[Nineq,Nspin,Nspin,Norb,Norb,Lmats],"chi2_fitgf","func")
  call assert_shape(Hloc,[Nineq,Nspin,Nspin,Norb,Norb],"chi2_fitgf","hloc")
  call ed_chi2_fitgf(func,bath,hloc,ispin)
end subroutine chi2_fitgf_ineq_c
