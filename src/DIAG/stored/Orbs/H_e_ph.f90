  !We build the electronic part of the electron-phonon interaction: Sum_iorb g_iorb n_iorb
  do i=MpiIstart,MpiIend
     call state2indices(i,[DimUps,DimDws],Indices)
     do iud=1,Ns_Ud
        mup = Hsector%H(iud)%map(Indices(iud))
        mdw = Hsector%H(iud+Ns_Ud)%map(Indices(iud+Ns_ud))
        Nups(iud,:) = Bdecomp(mup,Ns_Orb) ![Norb,1+Nbath]
        Ndws(iud,:) = Bdecomp(mdw,Ns_Orb)
     enddo
     Nup = Breorder(Nups)
     Ndw = Breorder(Ndws)
     !
     htmp = zero
     do iorb=1,Norb
        htmp = htmp + g_ph(iorb)*(nup(iorb)+ndw(iorb) - 1.d0)
     enddo
     !
     select case(MpiStatus)
     case (.true.)
        call sp_insert_element(MpiComm,spH0e_eph,htmp,i,i)
     case (.false.)
        call sp_insert_element(spH0e_eph,htmp,i,i)
     end select
     !
  enddo

  !Here we build the phononc part of the electron-phonon interaction: (b^+ + b)
  htmp = zero
  do iph=1,DimPh
     i = iph + 1
     if(i <= DimPh) then
        htmp = sqrt(dble(iph))
        call sp_insert_element(spH0ph_eph,htmp,iph,i)
     end if
     i = iph - 1
     if(i>0) then
        htmp = sqrt(dble(iph - 1))
        call sp_insert_element(spH0ph_eph,htmp,iph,i)
     end if
  end do













