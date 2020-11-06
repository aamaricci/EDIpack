  if(ph_type==1) then
     !phonon coupling to electron density operators
     do i=1,Nloc
        i_el = mod(i-1,DimUp*MpiQdw) + 1
        iph = (i-1)/(DimUp*MpiQdw) + 1
        !
        iup = iup_index(i_el+mpiIshift,DimUp)
        idw = idw_index(i_el+mpiIshift,DimUp)
        !
        mup = Hsector%H(1)%map(iup)
        mdw = Hsector%H(2)%map(idw)
        !
        nup = bdecomp(mup,Ns)
        ndw = bdecomp(mdw,Ns)
        !
        htmp=zero
        do iorb=1,Norb
           htmp = htmp + g_ph(iorb)*(nup(iorb)+ndw(iorb) - 1.d0)
        enddo
        !
        do jj = 1,DimPh
           if(jj .eq. iph+1) then
              j = i_el + (jj-1)*DimUp*MpiQdw 
              Hv(i) = Hv(i) + htmp*sqrt(dble(iph))*vin(j)
           endif
           !
           if(jj .eq. iph-1) then
              j = i_el + (jj-1)*DimUp*MpiQdw
              Hv(i) = Hv(i) + htmp*sqrt(dble(iph-1))*vin(j)
           endif
        enddo
     enddo
  elseif(ph_type==2) then
     !phonon coupling to orbital hybridization
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
     !Up spin component
     do iph=1,DimPh
        do jdw=1,MpiQdw
           do jup=1,DimUp
              mup  = Hsector%H(1)%map(jup)
              nup  = bdecomp(mup,Ns)
              !
              j    = jup + (jdw-1)*dimUp + (iph-1)*DimUp*MpiQdw
              !
              do iorb=1,Norb
                 do jorb=1,Norb
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
                          if(jj .eq. iph+1) then
                             i    = iup + (idw-1)*dimUp + (jj-1)*DimUp*MpiQdw
                             Hv(i) = Hv(i) + htmp*sqrt(dble(iph))*vin(j)
                          endif
                          !
                          if(jj .eq. iph-1) then
                             i    = iup + (idw-1)*dimUp + (jj-1)*DimUp*MpiQdw
                             Hv(i) = Hv(i) + htmp*sqrt(dble(iph-1))*vin(j)
                          endif
                       enddo
                       !
                    endif
                 enddo
              enddo
           enddo
        enddo
     enddo
     !Down spin component
     allocate(vt(mpiQup*DimDw*DimPh))
     allocate(Hvt(mpiQup*DimDw*DimPh))
     vt=0d0
     Hvt=0d0
     do iph=1,DimPh
        i_start = 1 + (iph-1)*DimUp*MpiQdw
        i_end = iph*DimUp*MpiQdw
        i_start2 = 1 + (iph-1)*MpiQup*DimDw
        i_end2 = iph*MpiQup*DimDw
        call vector_transpose_MPI(DimUp,MpiQdw,Vin(i_start:i_end),DimDw,MpiQup,vt(i_start2:i_end2)) !Vin^T --> Vt
     enddo
     !
     do iph=1,DimPh
        do jdw=1,MpiQup
           do jup=1,DimDw
              mdw  = Hsector%H(2)%map(jup)
              ndw  = bdecomp(mdw,Ns)
              !
              j    = jup + (jdw-1)*DimDw + (iph-1)*DimDw*MpiQup
              !
              do iorb=1,Norb
                 do jorb=1,Norb
                    Jcondition = &
                    (g_matrix(iorb,jorb)/=zero) .AND. &
                    (ndw(jorb)==1) .AND. (ndw(iorb)==0)
                    if (Jcondition) then
                       call c(jorb,mdw,k1,sg1)
                       call cdg(iorb,k1,k2,sg2)
                       iup = binary_search(Hsector%H(2)%map,k2)
                       idw = jdw
                       htmp = g_matrix(iorb,jorb)*sg1*sg2
                       !
                       do jj = 1,DimPh
                          if(jj .eq. iph+1) then
                             i = iup + (idw-1)*DimDw + (jj-1)*DimDw*MpiQup 
                             Hvt(i) = Hvt(i) + htmp*sqrt(dble(iph))*vt(j)
                          endif
                          !
                          if(jj .eq. iph-1) then
                             i = iup + (idw-1)*DimDw + (jj-1)*DimDw*MpiQup
                             Hvt(i) = Hvt(i) + htmp*sqrt(dble(iph-1))*vt(j)
                          endif
                       enddo
                       !
                    endif
                 enddo
              enddo
           enddo
        enddo
     enddo
     !
     deallocate(vt) ; allocate(vt(DimUp*mpiQdw*DimPh)) ;vt=0d0         !reallocate Vt
     do iph=1,DimPh
        i_start = 1 + (iph-1)*DimUp*MpiQdw
        i_end = iph*DimUp*MpiQdw
        i_start2 = 1 + (iph-1)*MpiQup*DimDw
        i_end2 = iph*MpiQup*DimDw
        call vector_transpose_MPI(DimDw,mpiQup,Hvt(i_start2:i_end2),DimUp,mpiQdw,vt(i_start:i_end)) !Hvt^T --> Vt
     enddo
     Hv = Hv + Vt
     !
     deallocate(vt,Hvt)
  endif
  
  
  
