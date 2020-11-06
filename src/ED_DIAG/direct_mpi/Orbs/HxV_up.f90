  do iph=1,DimPh
     do jjdw=1,MpiQdw
        do jjup=1,DimUp
           j = jjup + (jjdw-1)*DimUp
           call state2indices(j,[DimUps,DimDws],Jndices)
           j = j + (iph-1)*DimUp*MpiQdw
           !
           !>H_hyb: hopping terms for a given spin (imp <--> bath)
           do iorb=1,Ns_Ud
              mup = Hsector%H(iorb)%map(Jndices(iorb))
              Nups(iorb,:) = Bdecomp(mup,Ns_Orb)
              !
              do kp=1,Nbath
                 ialfa=1+kp
                 !
                 if( (diag_hybr(1,iorb,kp)/=0d0) &
                      .AND. (Nups(iorb,1)==1) .AND. (Nups(iorb,ialfa)==0) )then
                    call c(1,mup,k1,sg1)
                    call cdg(ialfa,k1,k2,sg2)
                    Indices       = Jndices
                    Indices(iorb) = binary_search(Hsector%H(iorb)%map,k2)
                    call indices2state(Indices,[DimUps,DimDws],i)
                    htmp = diag_hybr(1,iorb,kp)*sg1*sg2
                    !
                    i = i + (iph-1)*DimUp*MpiQdw
                    Hv(i) = Hv(i) + htmp*Vin(j)
                    !
                 endif
                 !
                 if( (diag_hybr(1,iorb,kp)/=0d0) &
                      .AND. (Nups(iorb,1)==0) .AND. (Nups(iorb,ialfa)==1) )then
                    call c(ialfa,mup,k1,sg1)
                    call cdg(1,k1,k2,sg2)
                    Indices       = Jndices
                    Indices(iorb) = binary_search(Hsector%H(iorb)%map,k2)
                    call indices2state(Indices,[DimUps,DimDws],i)
                    htmp = diag_hybr(1,iorb,kp)*sg1*sg2
                    !
                    i = i + (iph-1)*DimUp*MpiQdw
                    Hv(i) = Hv(i) + htmp*Vin(j)
                    !
                 endif
              enddo
           enddo
           !
        enddo
     enddo
  enddo





