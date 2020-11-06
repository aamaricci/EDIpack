  do i=1,Nloc
     i_el = mod(i-1,DimUp*MpiQdw) + 1
     !
     call state2indices(i_el+mpiIshift,[DimUps,DimDws],Indices)
     do iud=1,Ns_Ud
        mup = Hsector%H(iud)%map(Indices(iud))
        mdw = Hsector%H(iud+Ns_Ud)%map(Indices(iud+Ns_ud))
        Nups(iud,:) = Bdecomp(mup,Ns_Orb) ![1+Nbath]*Norb
        Ndws(iud,:) = Bdecomp(mdw,Ns_Orb) ![1+Nbath]*Norb
     enddo
     Nup = Breorder(Nups)
     Ndw = Breorder(Ndws)
     !
     !
     ! HxV_imp: Diagonal Elements, i.e. local part
     htmp = zero
     do iorb=1,Norb
        htmp = htmp + impHloc(1,1,iorb,iorb)*Nup(iorb)
        htmp = htmp + impHloc(Nspin,Nspin,iorb,iorb)*Ndw(iorb)
        htmp = htmp - xmu*(Nup(iorb)+Ndw(iorb))
     enddo
     !
     !
     !> H_Int: Kanamori interaction part. non-local S-E and P-H terms commented below.
     !
     ! density-density interaction: same orbital, opposite spins:
     !  = \sum_\a U_\a*(n_{\a,up}*n_{\a,dw})
     do iorb=1,Norb
        htmp = htmp + Uloc(iorb)*Nup(iorb)*Ndw(iorb)
     enddo
     if(Norb>1)then
        !density-density interaction: different orbitals, opposite spins:
        ! =   U'   *     sum_{i/=j} [ n_{i,up}*n_{j,dw} + n_{j,up}*n_{i,dw} ]
        ! =  (Uloc-2*Jh)*sum_{i/=j} [ n_{i,up}*n_{j,dw} + n_{j,up}*n_{i,dw} ]
        do iorb=1,Norb
           do jorb=iorb+1,Norb
              htmp = htmp + Ust*(Nup(iorb)*Ndw(jorb) + Nup(jorb)*Ndw(iorb))
           enddo
        enddo
        !density-density interaction: different orbitals, parallel spins
        ! = \sum_{i<j}    U''     *[ n_{i,up}*n_{j,up} + n_{i,dw}*n_{j,dw} ]
        ! = \sum_{i<j} (Uloc-3*Jh)*[ n_{i,up}*n_{j,up} + n_{i,dw}*n_{j,dw} ]
        do iorb=1,Norb
           do jorb=iorb+1,Norb
              htmp = htmp + (Ust-Jh)*(Nup(iorb)*Nup(jorb) + Ndw(iorb)*Ndw(jorb))
           enddo
        enddo
     endif
     !if using the Hartree-shifted chemical potential: mu=0 for half-filling
     !sum up the contributions of hartree terms:
     if(hfmode)then
        do iorb=1,Norb
           htmp = htmp - 0.5d0*Uloc(iorb)*(Nup(iorb)+Ndw(iorb)) + 0.25d0*uloc(iorb)
        enddo
        if(Norb>1)then
           do iorb=1,Norb
              do jorb=iorb+1,Norb
                 htmp=htmp-0.5d0*Ust*(Nup(iorb)+Ndw(iorb)+Nup(jorb)+Ndw(jorb))+0.25d0*Ust
                 htmp=htmp-0.5d0*(Ust-Jh)*(Nup(iorb)+Ndw(iorb)+Nup(jorb)+Ndw(jorb))+0.25d0*(Ust-Jh)
              enddo
           enddo
        endif
     endif
     !
     !
     !> H_Bath: local bath energy contribution.
     !diagonal bath hamiltonian: +energy of the bath=\sum_a=1,Norb\sum_{l=1,Nbath}\e^a_l n^a_l
     do iorb=1,size(bath_diag,2)
        do kp=1,Nbath
           ialfa = getBathStride(iorb,kp)
           htmp =htmp + bath_diag(1    ,iorb,kp)*Nup(ialfa) !UP
           htmp =htmp + bath_diag(Nspin,iorb,kp)*Ndw(ialfa) !DW
        enddo
     enddo
     !
     hv(i) = hv(i) + htmp*vin(i)
     !
  enddo



