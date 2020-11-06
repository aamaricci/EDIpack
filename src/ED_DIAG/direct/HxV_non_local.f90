  do j=1,Nloc
     j_el = mod(j-1,DimUp*DimDw) + 1 !electronic index
     iph = (j-1)/(DimUp*DimDw) + 1   !phononic index
     !
     jup = iup_index(j_el,DimUp)
     jdw = idw_index(j_el,DimUp)
     !
     mup = Hsector%H(1)%map(jup)
     mdw = Hsector%H(2)%map(jdw)
     !
     nup = bdecomp(mup,Ns)
     ndw = bdecomp(mdw,Ns)
     !
     ! SPIN-EXCHANGE (S-E) and PAIR-HOPPING TERMS
     !    S-E: J c^+_iorb_up c^+_jorb_dw c_iorb_dw c_jorb_up  (i.ne.j) 
     !    S-E: J c^+_{iorb} c^+_{jorb+Ns} c_{iorb+Ns} c_{jorb}
     if(Jhflag.AND.Jx/=0d0)then
        do iorb=1,Norb
           do jorb=1,Norb
              Jcondition=(&
                   (iorb/=jorb).AND.&
                   (nup(jorb)==1).AND.&
                   (ndw(iorb)==1).AND.&
                   (ndw(jorb)==0).AND.&
                   (nup(iorb)==0))
              if(Jcondition)then
                 call c(iorb,mdw,k1,sg1)  !DW
                 call cdg(jorb,k1,k2,sg2) !DW
                 idw=binary_search(Hsector%H(2)%map,k2)
                 call c(jorb,mup,k3,sg3)  !UP
                 call cdg(iorb,k3,k4,sg4) !UP
                 iup=binary_search(Hsector%H(1)%map,k4)
                 htmp = Jx*sg1*sg2*sg3*sg4
                 i = iup + (idw-1)*dimup + (iph-1)*DimUp*DimDw
                 !
                 Hv(i) = Hv(i) + htmp*vin(j)
                 !
              endif
           enddo
        enddo
     endif
     !
     ! PAIR-HOPPING (P-H) TERMS
     !    P-H: J c^+_iorb_up c^+_iorb_dw   c_jorb_dw   c_jorb_up  (i.ne.j) 
     !    P-H: J c^+_{iorb}  c^+_{iorb+Ns} c_{jorb+Ns} c_{jorb}
     if(Jhflag.AND.Jp/=0d0)then
        do iorb=1,Norb
           do jorb=1,Norb
              Jcondition=(&
                   (nup(jorb)==1).AND.&
                   (ndw(jorb)==1).AND.&
                   (ndw(iorb)==0).AND.&
                   (nup(iorb)==0))
              if(Jcondition)then
                 call c(jorb,mdw,k1,sg1)       !c_jorb_dw
                 call cdg(iorb,k1,k2,sg2)      !c^+_iorb_dw
                 idw = binary_search(Hsector%H(2)%map,k2)
                 call c(jorb,mup,k3,sg3)       !c_jorb_up
                 call cdg(iorb,k3,k4,sg4)      !c^+_iorb_up
                 iup = binary_search(Hsector%H(1)%map,k4)
                 htmp = Jp*sg1*sg2*sg3*sg4
                 i = iup + (idw-1)*dimup + (iph-1)*DimUp*DimDw
                 !
                 Hv(i) = Hv(i) + htmp*vin(j)
                 !
              endif
           enddo
        enddo
     endif
     !
  enddo
