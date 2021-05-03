!We build the electronic part of the electron-phonon interaction: Sum_iorb g_iorb n_iorb
  !Coupling to combinations of electron density operators
  if(ph_type==1) then
     do i=MpiIstart,MpiIend
        iup = iup_index(i,DimUp)
        idw = idw_index(i,DimUp)
        !
        mup = Hsector%H(1)%map(iup)
        mdw = Hsector%H(2)%map(idw)
        !
        nup = bdecomp(mup,Ns)
        ndw = bdecomp(mdw,Ns)
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
     !
  !Coupling to orbital hybridization
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
     ! Up spin component
     do jup=1,DimUp
        mup  = Hsector%H(1)%map(jup)
        Nup  = Bdecomp(mup,Ns)
        !
        do iorb=1,Norb
           do jorb=1,Norb
              Jcondition = &
                   (g_matrix(iorb,jorb)/=zero) .AND. &
                   (Nup(jorb)==1) .AND. (Nup(iorb)==0)
              if (Jcondition) then
                 call c(jorb,mup,k1,sg1)
                 call cdg(iorb,k1,k2,sg2)
                 iup = binary_search(Hsector%H(1)%map,k2)
                 htmp = g_matrix(iorb,jorb)*sg1*sg2
                 !
                 call sp_insert_element(spH0e_eph,htmp,iup,jup)
                 !
              endif
           enddo
        enddo 
        !
     enddo
     !
     ! Down spin component
     do jdw=1,DimDw
        mdw  = Hsector%H(2)%map(jdw)
        Ndw  = bdecomp(mdw,Ns)
        !
        do iorb=1,Norb
           do jorb=1,Norb
              Jcondition = &
                   (g_matrix(iorb,jorb)/=zero) .AND. &
                   (Ndw(jorb)==1) .AND. (Ndw(iorb)==0)
              if (Jcondition) then
                 call c(jorb,mdw,k1,sg1)
                 call cdg(iorb,k1,k2,sg2)
                 idw = binary_search(Hsector%H(2)%map,k2)
                 htmp = g_matrix(iorb,jorb)*sg1*sg2
                 !
                 call sp_insert_element(spH0edw_eph,htmp,idw,jdw)
                 !
              endif
           enddo
        enddo
        !
     enddo
  endif


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













