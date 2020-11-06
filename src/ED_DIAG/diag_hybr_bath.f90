  
  if(allocated(diag_hybr))deallocate(diag_hybr)
  if(allocated(bath_diag))deallocate(bath_diag)
  select case (bath_type)
  case default
     Nfoo = size(dmft_bath%e,2)
     allocate(diag_hybr(Nspin,Norb,Nbath));diag_hybr=0d0
     allocate(bath_diag(Nspin,Nfoo,Nbath));bath_diag=0d0
     do ibath=1,Nbath
        do ispin=1,Nspin             
           do iorb=1,Norb
              diag_hybr(ispin,iorb,ibath)=dmft_bath%v(ispin,iorb,ibath)
           enddo
           do iorb=1,Nfoo
              bath_diag(ispin,iorb,ibath)=dmft_bath%e(ispin,iorb,ibath)
           enddo
        enddo
     enddo
  case ("replica")
     allocate(diag_hybr(Nspin,Norb,Nbath));diag_hybr=0d0
     allocate(bath_diag(Nspin,Norb,Nbath));bath_diag=0d0
     do ibath=1,Nbath
        Hbath_reconstructed = bath_from_sym(dmft_bath%item(ibath)%lambda)
        do ispin=1,Nspin
           do iorb=1,Norb
              diag_hybr(ispin,iorb,ibath)=dmft_bath%item(ibath)%v(ispin)
              bath_diag(ispin,iorb,ibath)=Hbath_reconstructed(ispin,ispin,iorb,iorb)
           enddo
        enddo
     enddo
  end select
