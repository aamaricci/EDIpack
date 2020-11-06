  if(ph_type==1) then
     do i=1,Nloc
        i_el = mod(i-1,DimUp*DimDw) + 1
        iph = (i-1)/(DimUp*DimDw) + 1
        !
        iup = iup_index(i_el,DimUp)
        idw = idw_index(i_el,DimUp)
        !
        mup = Hsector%H(1)%map(iup)
        mdw = Hsector%H(2)%map(idw)
        !
        nup = bdecomp(mup,Ns)
        ndw = bdecomp(mdw,Ns)
        !
        htmp=zero
        do iorb=1,Norb
           htmp = htmp + g_ph(iorb)*(nup(iorb)+ndw(iorb) - 1.d0) !electronin part
        enddo
        !
        do jj = 1,DimPh
           if(jj .eq. iph+1) then !destruction of a phonon (read from right to left)
              j = i_el + (jj-1)*DimUp*DimDw
              Hv(i) = Hv(i) + htmp*sqrt(dble(iph))*vin(j)
           endif
           !
           if(jj .eq. iph-1) then  !creation of a phonon
              j = i_el + (jj-1)*DimUp*DimDw
              Hv(i) = Hv(i) + htmp*sqrt(dble(iph-1))*vin(j)
           endif
        enddo
     enddo
     !
  elseif(ph_type==2) then
     i = 1
     g_matrix = 0.d0     !matrix of electron-phonon coupling constants
     if(Norb>1) then
        do iorb=1,Norb
           do jorb=iorb+1,Norb
              g_matrix(iorb,jorb) = g_ph(i)
              g_matrix(jorb,iorb) = g_matrix(iorb,jorb)
              i = i + 1
           enddo
        enddo
     endif
     !
     do iph=1,DimPh
        do jdw=1,DimDw
           do jup=1,DimUp
              mup  = Hsector%H(1)%map(jup)
              nup  = bdecomp(mup,Ns)
              !
              mdw  = Hsector%H(2)%map(jdw)
              ndw  = bdecomp(mdw,Ns)
              !
              j    = jup + (jdw-1)*dimUp + (iph-1)*DimUp*DimDw
              !
              do iorb=1,Norb
                 do jorb=1,Norb
                    !Up spin part
                    Jcondition = &
                         (g_matrix(iorb,jorb)/=zero) .AND. &
                         (nup(jorb)==1) .AND. (nup(iorb)==0)
                    if (Jcondition) then
                       call c(jorb,mup,k1,sg1)
                       call cdg(iorb,k1,k2,sg2)
                       iup  = binary_search(Hsector%H(1)%map,k2)
                       idw  = jdw 
                       htmp = g_matrix(iorb,jorb)*sg1*sg2
                       !
                       do jj = 1,DimPh
                          if(jj .eq. iph+1) then !creation
                             i = iup + (idw-1)*DimUp + (jj-1)*DimUp*DimDw
                             Hv(i) = Hv(i) + htmp*sqrt(dble(iph))*vin(j)
                          endif
                          !
                          if(jj .eq. iph-1) then !annihilation
                             i = iup + (idw-1)*DimUp + (jj-1)*DimUp*DimDw
                             Hv(i) = Hv(i) + htmp*sqrt(dble(iph-1))*vin(j)
                          endif
                       enddo
                       !
                    endif
                    !Down spin part
                    Jcondition = &
                         (g_matrix(iorb,jorb)/=zero) .AND. &
                         (ndw(jorb)==1) .AND. (ndw(iorb)==0)
                    if (Jcondition) then
                       call c(jorb,mdw,k1,sg1)
                       call cdg(iorb,k1,k2,sg2)
                       idw = binary_search(Hsector%H(2)%map,k2)
                       iup = jup
                       htmp = g_matrix(iorb,jorb)*sg1*sg2
                       !
                       do jj = 1,DimPh
                          if(jj .eq. iph+1) then !creation
                             i = iup + (idw-1)*DimUp + (jj-1)*DimUp*DimDw
                             Hv(i) = Hv(i) + htmp*sqrt(dble(iph))*vin(j)
                          endif
                          !
                          if(jj .eq. iph-1) then !annihilation
                             i = iup + (idw-1)*DimUp + (jj-1)*DimUp*DimDw
                             Hv(i) = Hv(i) + htmp*sqrt(dble(iph-1))*vin(j)
                          endif
                       enddo
                       !
                    endif
                 enddo
              enddo
              !
           enddo
        enddo
     enddo
  endif
  
  
