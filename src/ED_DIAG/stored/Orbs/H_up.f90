  !>H_hyb: hopping terms for a given spin (imp <--> bath)
  do iorb=1,Ns_Ud               !==Norb
     do jup=1,DimUps(iorb)
        mup          = Hsector%H(iorb)%map(jup)
        Nups(iorb,:) = bdecomp(mup,Ns_Orb) ![Norb,1+Nbath]
        !
        do kp=1,Nbath
           ialfa=1+kp
           if( (diag_hybr(1,iorb,kp)/=0d0) &
                .AND. (Nups(iorb,1)==1) .AND. (Nups(iorb,ialfa)==0) )then              
              call c(1,mup,k1,sg1)
              call cdg(ialfa,k1,k2,sg2)
              iup  = binary_search(Hsector%H(iorb)%map,k2)
              htmp = diag_hybr(1,iorb,kp)*sg1*sg2
              !
              call sp_insert_element(spH0ups(iorb),htmp,iup,jup)
              !
           endif
           !
           if( (diag_hybr(1,iorb,kp)/=0d0) &
                .AND. (Nups(iorb,1)==0) .AND. (Nups(iorb,ialfa)==1) )then
              call c(ialfa,mup,k1,sg1)
              call cdg(1,k1,k2,sg2)
              iup  = binary_search(Hsector%H(iorb)%map,k2)
              htmp = diag_hybr(1,iorb,kp)*sg1*sg2
              !
              call sp_insert_element(spH0ups(iorb),htmp,iup,jup)
              !
           endif
        enddo
        !
     enddo
  enddo

