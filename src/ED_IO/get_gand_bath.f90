function delta_bath_user(x,bath_) result(Delta)
  complex(8),dimension(:),intent(in)                  :: x
  type(effective_bath)                                :: dmft_bath_
  complex(8),dimension(Nspin,Nspin,Norb,Norb,size(x)) :: Delta
  real(8),dimension(:)                                :: bath_
  logical                                             :: check
  check= check_bath_dimension(bath_)
  if(.not.check)stop "delta_bath_mats_main_ error: wrong bath dimensions"
  call allocate_dmft_bath(dmft_bath_)
  call set_dmft_bath(bath_,dmft_bath_)
  Delta = delta_bath_function(x,dmft_bath_)
  call deallocate_dmft_bath(dmft_bath_)
end function delta_bath_user


function g0and_bath_user(x,bath_) result(G0and)
  complex(8),dimension(:),intent(in)                  :: x
  type(effective_bath)                                :: dmft_bath_
  complex(8),dimension(Nspin,Nspin,Norb,Norb,size(x)) :: G0and
  real(8),dimension(:)                                :: bath_
  logical                                             :: check
  check= check_bath_dimension(bath_)
  if(.not.check)stop "g0and_bath_mats_main_ error: wrong bath dimensions"
  call allocate_dmft_bath(dmft_bath_)
  call set_dmft_bath(bath_,dmft_bath_)
  G0and = g0and_bath_function(x,dmft_bath_)
  call deallocate_dmft_bath(dmft_bath_)
end function g0and_bath_user



function invg0_bath_user(x,bath_) result(G0and)
  complex(8),dimension(:),intent(in)                  :: x
  type(effective_bath)                                :: dmft_bath_
  complex(8),dimension(Nspin,Nspin,Norb,Norb,size(x)) :: G0and
  real(8),dimension(:)                                :: bath_
  logical                                             :: check
  check= check_bath_dimension(bath_)
  if(.not.check)stop "invg0_bath_mats_main_ error: wrong bath dimensions"
  call allocate_dmft_bath(dmft_bath_)
  call set_dmft_bath(bath_,dmft_bath_)
  G0and = invg0_bath_function(x,dmft_bath_)
  call deallocate_dmft_bath(dmft_bath_)
end function invg0_bath_user

