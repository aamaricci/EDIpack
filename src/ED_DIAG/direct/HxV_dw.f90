  do iph=1,DimPh
     do jup=1,DimUp
        do jdw=1,DimDw
           mdw  = Hsector%H(2)%map(jdw)
           ndw  = bdecomp(mdw,Ns)
           j    = jup + (jdw-1)*dimUp + (iph-1)*DimUp*DimDw
           !
           !
           !> H_imp: Off-diagonal elements, i.e. non-local part. 
           !remark: iorb=jorb cant have simultaneously n=0 and n=1 (Jcondition)
           do iorb=1,Norb
              do jorb=1,Norb
                 Jcondition = &
                      (impHloc(Nspin,Nspin,iorb,jorb)/=zero) .AND. &
                      (ndw(jorb)==1) .AND. (ndw(iorb)==0)
                 if (Jcondition) then
                    call c(jorb,mdw,k1,sg1)
                    call cdg(iorb,k1,k2,sg2)
                    idw = binary_search(Hsector%H(2)%map,k2)
                    iup = jup
                    i   = iup + (idw-1)*DimUp + (iph-1)*DimUp*DimDw
                    htmp = impHloc(Nspin,Nspin,iorb,jorb)*sg1*sg2
                    !
                    Hv(i) = Hv(i) + htmp*vin(j)
                    !
                 endif
              enddo
           enddo
           !
           !
           !> H_Bath: inter-orbital bath hopping contribution.
           if(bath_type=="replica") then
              do kp=1,Nbath
                 do iorb=1,Norb
                    do jorb=1,Norb
                       !
                       ialfa = getBathStride(iorb,kp)
                       ibeta = getBathStride(jorb,kp)
                       Jcondition = &
                            (hbath_tmp(Nspin,Nspin,iorb,jorb,kp)/=zero) .AND. &
                            (ndw(ibeta)==1) .AND. (ndw(ialfa)==0)
                       !
                       if (Jcondition)then
                          call c(ibeta,mdw,k1,sg1)
                          call cdg(ialfa,k1,k2,sg2)
                          idw = binary_search(Hsector%H(2)%map,k2)
                          iup = jup
                          i   = iup + (idw-1)*DimUp + (iph-1)*DimUp*DimDw
                          htmp = hbath_tmp(Nspin,Nspin,iorb,jorb,kp)*sg1*sg2
                          !
                          Hv(i) = Hv(i) + htmp*vin(j)
                          !
                       endif
                    enddo
                 enddo
              enddo
           end if
           !
           !
           !>H_hyb: hopping terms for a given spin (imp <--> bath)
           do iorb=1,Norb
              do kp=1,Nbath
                 ialfa=getBathStride(iorb,kp)
                 !
                 if( (diag_hybr(Nspin,iorb,kp)/=0d0) .AND. &
                      (ndw(iorb)==1) .AND. (ndw(ialfa)==0) )then
                    call c(iorb,mdw,k1,sg1)
                    call cdg(ialfa,k1,k2,sg2)
                    idw = binary_search(Hsector%H(2)%map,k2)
                    iup = jup
                    i   = iup + (idw-1)*DimUp + (iph-1)*DimUp*DimDw
                    htmp=diag_hybr(Nspin,iorb,kp)*sg1*sg2
                    !
                    Hv(i) = Hv(i) + htmp*vin(j)
                    !
                 endif
                 if( (diag_hybr(Nspin,iorb,kp)/=0d0) .AND. &
                      (ndw(iorb)==0) .AND. (ndw(ialfa)==1) )then
                    call c(ialfa,mdw,k1,sg1)
                    call cdg(iorb,k1,k2,sg2)
                    idw = binary_search(Hsector%H(2)%map,k2)
                    iup = jup
                    i   = iup + (idw-1)*DimUp + (iph-1)*DimUp*DimDw
                    htmp=diag_hybr(Nspin,iorb,kp)*sg1*sg2
                    !
                    Hv(i) = Hv(i) + htmp*vin(j)
                    !
                 endif
              enddo
           enddo
           !
           !
        end do
     enddo
  enddo
