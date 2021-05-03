  do i=1,Nloc
     i_el = mod(i-1,DimUp*DimDw) + 1
     iph = (i-1)/(DimUp*DimDw) + 1
     !
     call state2indices(i_el,[DimUps,DimDws],Indices)
     do iud=1,Ns_Ud
        mup = Hsector%H(iud)%map(Indices(iud))
        mdw = Hsector%H(iud+Ns_Ud)%map(Indices(iud+Ns_ud))
        Nups(iud,:) = Bdecomp(mup,Ns_Orb) ![Norb,1+Nbath]
        Ndws(iud,:) = Bdecomp(mdw,Ns_Orb)
     enddo
     Nup = Breorder(Nups)
     Ndw = Breorder(Ndws)
     !
     htmp=zero
     do iorb=1,Norb
        htmp = htmp + g_ph(iorb)*(nup(iorb)+ndw(iorb) - 1.d0) !electronin part
     enddo
     !
     do jj = 1,DimPh
        if(jj .eq. iph+1) then
           j = i_el + (jj-1)*DimUp*DimDw
           Hv(i) = Hv(i) + htmp*sqrt(dble(iph))*vin(j)
        endif
        !
        if(jj .eq. iph-1) then
           j = i_el + (jj-1)*DimUp*DimDw
           Hv(i) = Hv(i) + htmp*sqrt(dble(iph-1))*vin(j)
        endif
     enddo
  enddo
