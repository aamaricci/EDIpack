!Phononic hamiltonian H_ph = w0 b^+ b
  do iph=1,DimPh
     htmp = w0_ph*(iph-1)
     call sp_insert_element(spH0_ph,htmp,iph,iph)
  end do

