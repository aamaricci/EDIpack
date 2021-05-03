  do jjdw=1,MpiQup
     do jjup=1,DimDw
        j = jjup + (jjdw-1)*DimDw
        call state2indices(j,[DimDws,DimUps],Jndices)
        !
        !>H_hyb: hopping terms for a given spin (imp <--> bath)
        do iorb=1,Ns_Ud           
           mdw  = Hsector%H(iorb+Ns_Ud)%map( Jndices(iorb) )
           Ndws(iorb,:) = bdecomp(mdw,Ns_Orb)
           !
           do kp=1,Nbath
              ialfa=1+kp
              !
              if( (diag_hybr(Nspin,iorb,kp)/=0d0) &
                   .AND. (Ndws(iorb,1)==1) .AND. (Ndws(iorb,ialfa)==0) )then
                 call c(1,mdw,k1,sg1)
                 call cdg(ialfa,k1,k2,sg2)                 
                 Indices       = Jndices
                 Indices(iorb) = binary_search(Hsector%H(iorb+Ns_Ud)%map,k2)
                 call indices2state(Indices,[DimDws,DimUps],i)
                 htmp = diag_hybr(Nspin,iorb,kp)*sg1*sg2
                 !
                 Hvt(i) = Hvt(i) + htmp*Vt(j)
                 !
              endif
              !
              if( (diag_hybr(Nspin,iorb,kp)/=0d0) &
                   .AND. (Ndws(iorb,1)==0) .AND. (Ndws(iorb,ialfa)==1) )then
                 call c(ialfa,mdw,k1,sg1)
                 call cdg(1,k1,k2,sg2)
                 Indices       = Jndices
                 Indices(iorb) = binary_search(Hsector%H(iorb+Ns_Ud)%map,k2)
                 call indices2state(Indices,[DimDws,DimUps],i)
                 htmp = diag_hybr(Nspin,iorb,kp)*sg1*sg2
                 !
                 Hvt(i) = Hvt(i) + htmp*Vt(j)
                 !
              endif
              !
           enddo
        enddo
        !
     enddo
  enddo






