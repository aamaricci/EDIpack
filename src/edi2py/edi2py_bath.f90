!ED_BATH:
integer(c_int) function get_bath_dimension_c() result(Nb) bind(c, name='get_bath_dimension')
  Nb=ed_get_bath_dimension()
end function get_bath_dimension_c

!H_REPLICA SETUP
subroutine init_Hreplica_direct_nn_c(Hloc) bind(c, name='init_Hreplica_direct_nn')
  real(c_double),dimension(Nspin,Nspin,Norb,Norb) :: Hloc
  call ed_set_Hreplica(hloc)
end subroutine init_Hreplica_direct_nn_c

subroutine init_Hreplica_direct_so_c(Hloc) bind(c, name='init_Hreplica_direct_so')
  real(c_double),dimension(Nspin*Norb,Nspin*Norb) :: Hloc 
  call ed_set_Hreplica(hloc)
end subroutine init_Hreplica_direct_so_c

subroutine init_Hreplica_symmetries_site_c(Hvec,lambdavec,Nsym) bind(c, name='init_Hreplica_symmetries_site')
  integer(c_int),value                                     :: Nsym
  real(c_double),dimension(Nspin,Nspin,Norb,Norb,Nsym)     :: Hvec
  real(c_double),dimension(Nsym)                           :: lambdavec
  call ed_set_Hreplica(Hvec,lambdavec)
end subroutine init_Hreplica_symmetries_site_C

subroutine init_Hreplica_symmetries_ineq_c(Hvec,lambdavec,Nlat,Nsym) bind(c, name='init_Hreplica_symmetries_ineq')
  integer(c_int),value                                     :: Nlat, Nsym
  real(c_double),dimension(Nspin,Nspin,Norb,Norb,Nsym)     :: Hvec
  real(c_double),dimension(Nlat,Nsym)                      :: lambdavec
  call ed_set_Hreplica(Hvec,lambdavec)
end subroutine init_Hreplica_symmetries_ineq_c




!BREAK_SYMMETRY_BATH
subroutine break_symmetry_bath_site_c(bath,dim_bath,field,sgn,sav) bind(c, name='break_symmetry_bath_site')
  integer(c_int64_t)                    :: dim_bath(1)
  real(c_double),dimension(dim_bath(1)) :: bath
  real(c_double),value                  :: field
  real(c_double),value                  :: sgn
  integer(c_int),value                  :: sav
  call ed_break_symmetry_bath(bath,field,sgn,i2l(sav))
end subroutine break_symmetry_bath_site_c
!
subroutine break_symmetry_bath_ineq_c(bath,dim_bath,field,sgn,sav) bind(c, name='break_symmetry_bath_ineq')
  integer(c_int64_t)                                :: dim_bath(2)
  real(c_double),dimension(dim_bath(1),dim_bath(2)) :: bath
  real(c_double),value                              :: field
  real(c_double),value                              :: sgn
  integer(c_int),value                              :: sav
  call ed_break_symmetry_bath(bath,field,sgn,i2l(sav))
end subroutine break_symmetry_bath_ineq_c



!SPIN_SYMMETRIZE_BATH
subroutine spin_symmetrize_bath_site_c(bath,dim_bath,sav) bind(c, name='spin_symmetrize_bath_site')
  integer(c_int64_t)                    :: dim_bath(1)
  integer(c_int),value                  :: sav
  real(c_double),dimension(dim_bath(1)) :: bath
  call ed_spin_symmetrize_bath(bath,i2l(sav))
end subroutine spin_symmetrize_bath_site_c
!
subroutine spin_symmetrize_bath_ineq_c(bath,dim_bath,sav) bind(c, name='spin_symmetrize_bath_ineq')
  integer(c_int64_t)                                :: dim_bath(2)
  integer(c_int),value                              :: sav
  real(c_double),dimension(dim_bath(1),dim_bath(2)) :: bath
  call ed_spin_symmetrize_bath(bath,i2l(sav))
end subroutine spin_symmetrize_bath_ineq_c



!ORB_SYMMETRIZE_BATH
subroutine orb_symmetrize_bath_site_c(bath,dim_bath,sav) bind(c, name='orb_symmetrize_bath_site')
  integer(c_int64_t)                    :: dim_bath(1)
  integer(c_int),value                  :: sav
  real(c_double),dimension(dim_bath(1)) :: bath
  call ed_orb_symmetrize_bath(bath,i2l(sav))
end subroutine orb_symmetrize_bath_site_c
!
subroutine orb_symmetrize_bath_ineq_c(bath,dim_bath,sav) bind(c, name='orb_symmetrize_bath_ineq')
  integer(c_int64_t)                                :: dim_bath(2)
  integer(c_int),value                              :: sav
  real(c_double),dimension(dim_bath(1),dim_bath(2)) :: bath
  call ed_orb_symmetrize_bath(bath,i2l(sav))
end subroutine orb_symmetrize_bath_ineq_c


!ORB_EQUALITY_BATH
subroutine orb_equality_bath_site_c(bath,dim_bath,indx,sav) bind(c, name='orb_equality_bath_site')
  integer(c_int64_t)                       :: dim_bath(1)
  real(c_double),dimension(dim_bath(1))    :: bath
  integer(c_int),value                     :: indx,sav
  call ed_orb_equality_bath(bath,indx,i2l(sav))
end subroutine orb_equality_bath_site_c
!
subroutine orb_equality_bath_ineq_c(bath,dim_bath,indx,sav) bind(c, name='orb_equality_bath_ineq')
  integer(c_int64_t)                                   :: dim_bath(2)
  real(c_double),dimension(dim_bath(1),dim_bath(2))    :: bath
  integer(c_int),value                                 :: indx,sav
  call ed_orb_equality_bath(bath,indx,i2l(sav))
end subroutine orb_equality_bath_ineq_c





!PH_SYMMETRIZE_BATH
subroutine ph_symmetrize_bath_site_c(bath,dim_bath,sav) bind(c, name='ph_symmetrize_bath_site')
  integer(c_int64_t)                    :: dim_bath(1)
  real(c_double),dimension(dim_bath(1)) :: bath
  integer(c_int),value                  :: sav
  call ed_ph_symmetrize_bath(bath,i2l(sav))
end subroutine ph_symmetrize_bath_site_c
!
subroutine ph_symmetrize_bath_ineq_c(bath,dim_bath,sav) bind(c, name='ph_symmetrize_bath_ineq')
  integer(c_int64_t)                                :: dim_bath(2)
  real(c_double),dimension(dim_bath(1),dim_bath(2)) :: bath
  integer(c_int),value                              :: sav
  call ed_ph_symmetrize_bath(bath,i2l(sav))
end subroutine ph_symmetrize_bath_ineq_c



!PH_TRANS_BATH
subroutine ph_trans_bath_site_c(bath,dim_bath,sav) bind(c, name='ph_trans_bath_site')
  integer(c_int64_t)                    :: dim_bath(1)
  real(c_double),dimension(dim_bath(1)) :: bath
  integer(c_int),value                  :: sav
  call ed_ph_trans_bath(bath,i2l(sav))
end subroutine ph_trans_bath_site_c
!
subroutine ph_trans_bath_ineq_c(bath,dim_bath,sav) bind(c, name='ph_trans_bath_ineq')
  integer(c_int64_t)                                :: dim_bath(2)
  real(c_double),dimension(dim_bath(1),dim_bath(1)) :: bath
  integer(c_int),value                              :: sav
  call ed_ph_trans_bath(bath,i2l(sav))
end subroutine ph_trans_bath_ineq_c



!BATH COMPONENT ROUTINES
subroutine get_bath_component_dimension_c(instr,Ndim) bind(c, name='get_bath_component_dimension')
  character(kind=c_char), dimension(1)   :: instr
  character(len=1)                       :: typ
  integer(c_int64_t)                     :: Ndim(3)
  typ(1:1)=instr(1)
  Ndim = ed_get_bath_component_dimension(typ)
end subroutine get_bath_component_dimension_c

subroutine get_bath_component(array,dim_array,bath,dim_bath,instr) bind(c, name='get_bath_component')
  integer(c_int64_t)                                               :: dim_bath(1),dim_array(3)
  real(c_double),dimension(dim_array(1),dim_array(2),dim_array(3)) :: array
  real(c_double),dimension(dim_bath(1))                            :: bath
  character(kind=c_char), dimension(1)                             :: instr
  character(len=1)                                                 :: typ
  typ(1:1)=instr(1)
  call ed_get_bath_component(array,bath,typ)
end subroutine get_bath_component

subroutine set_bath_component(array,dim_array,bath,dim_bath,instr) bind(c, name='set_bath_component')
  integer(c_int64_t)                                               :: dim_bath(1),dim_array(3)
  real(c_double),dimension(dim_array(1),dim_array(2),dim_array(3)) :: array
  real(c_double),dimension(dim_bath(1))                            :: bath
  character(kind=c_char), dimension(1)                             :: instr
  character(len=1)                                                 :: typ
  typ(1:1)=instr(1)
  call ed_set_bath_component(array,bath,typ)
end subroutine set_bath_component

subroutine copy_bath_component(bathIN,dim_bathin,bathOUT,dim_bathout,instr) bind(c, name='copy_bath_component')
  real(c_double),dimension(dim_bathin)      :: bathIN
  real(c_double),dimension(dim_bathout)     :: bathOUT
  character(kind=c_char), dimension(1)      :: instr
  character(len=1)                          :: typ
  integer(c_int),value                      :: dim_bathin, dim_bathout
  typ(1:1)=instr(1)
  call ed_copy_bath_component(bathIN,bathOUT,typ)
end subroutine copy_bath_component
