  !>H_hyb: hopping terms for a given spin (imp <--> bath)
  do iorb=1,Ns_Ud               !==Norb
     do jdw=1,DimDws(iorb)
        mdw          = Hsector%H(iorb+Ns_Ud)%map(jdw)
        Ndws(iorb,:) = bdecomp(mdw,Ns_Orb) ![Norb,1+Nbath]
        !
        do kp=1,Nbath
           ialfa=1+kp
           if( (diag_hybr(Nspin,iorb,kp)/=0d0) &
                .AND. (Ndws(iorb,1)==1) .AND. (Ndws(iorb,ialfa)==0) )then
              call c(1,mdw,k1,sg1)
              call cdg(ialfa,k1,k2,sg2)
              idw  = binary_search(Hsector%H(iorb+Ns_Ud)%map,k2)
              htmp = diag_hybr(Nspin,iorb,kp)*sg1*sg2
              !
              call sp_insert_element(spH0dws(iorb),htmp,idw,jdw)
              !
           endif
           !
           if( (diag_hybr(Nspin,iorb,kp)/=0d0) &
                .AND. (Ndws(iorb,1)==0) .AND. (Ndws(iorb,ialfa)==1) )then
              call c(ialfa,mdw,k1,sg1)
              call cdg(1,k1,k2,sg2)
              idw  = binary_search(Hsector%H(iorb+Ns_Ud)%map,k2)
              htmp = diag_hybr(Nspin,iorb,kp)*sg1*sg2
              !
              call sp_insert_element(spH0dws(iorb),htmp,idw,jdw)
              !
           endif
        enddo
        !
     enddo
  enddo
