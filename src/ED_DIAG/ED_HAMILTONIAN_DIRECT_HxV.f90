! > SPARSE MAT-VEC DIRECT ON-THE-FLY PRODUCT 
MODULE ED_HAMILTONIAN_DIRECT_HxV
  USE ED_HAMILTONIAN_COMMON
  implicit none
  private


  !>Sparse Mat-Vec direct on-the-fly product 
  public  :: directMatVec_main
  public  :: directMatVec_orbs
#ifdef _MPI
  public  :: directMatVec_MPI_main
  public  :: directMatVec_MPI_orbs
#endif



contains


  subroutine directMatVec_main(Nloc,vin,Hv)
    integer                                        :: Nloc
    real(8),dimension(Nloc)                        :: vin
    real(8),dimension(Nloc)                        :: Hv
    real(8),dimension(:),allocatable               :: vt,Hvt
    integer,dimension(2*Ns_Ud)                     :: Indices,Jndices ![2-2*Norb]
    integer,dimension(Ns_Ud,Ns_Orb)                :: Nups,Ndws       ![1,Ns]-[Norb,1+Nbath]
    integer,dimension(Ns)                          :: Nup,Ndw
    real(8),dimension(Nspin,Nspin,Norb,Norb,Nbath) :: Hbath_tmp
    real(8),dimension(Norb,Norb)                   :: g_matrix

    if(.not.Hsector%status)stop "directMatVec_cc ERROR: Hsector NOT allocated"
    isector=Hsector%index
    !
    if(Nloc/=getdim(isector))stop "directMatVec_cc ERROR: Nloc != dim(isector)"
    !
    !Get diagonal hybridization, bath energy
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
          Hbath_tmp(:,:,:,:,ibath) = bath_from_sym(dmft_bath%item(ibath)%lambda)
          do ispin=1,Nspin
             do iorb=1,Norb
                diag_hybr(ispin,iorb,ibath)=dmft_bath%item(ibath)%v
                bath_diag(ispin,iorb,ibath)=Hbath_tmp(ispin,ispin,iorb,iorb,ibath)
             enddo
          enddo
       enddo
    end select
    !
    Hv=0d0
    !
    !-----------------------------------------------!
    !LOCAL HAMILTONIAN PART: H_loc*vin = vout
    include "direct/HxV_local.f90"
    !
    !UP HAMILTONIAN TERMS
    include "direct/HxV_up.f90"
    !    
    !DW HAMILTONIAN TERMS
    include "direct/HxV_dw.f90"
    !
    if(DimPh>1) then
       !PHONON TERMS
       include "direct/HxV_ph.f90"
       !ELECTRON-PHONON INTERACTION
       include "direct/HxV_eph.f90"
    end if
    !
    !NON-LOCAL HAMILTONIAN PART: H_non_loc*vin = vout
    if(Jhflag)then
       include "direct/HxV_non_local.f90"
    endif
    !-----------------------------------------------!
    !
    deallocate(diag_hybr,bath_diag)
    return
  end subroutine directMatVec_main




  subroutine directMatVec_orbs(Nloc,vin,Hv)
    integer                                        :: Nloc
    real(8),dimension(Nloc)                        :: vin
    real(8),dimension(Nloc)                        :: Hv
    real(8),dimension(:),allocatable               :: vt,Hvt
    integer                                        :: isector
    integer,dimension(2*Ns_Ud)                     :: Indices,Jndices ![2-2*Norb]
    integer,dimension(Ns_Ud,Ns_Orb)                :: Nups,Ndws       ![1,Ns]-[Norb,1+Nbath]
    integer,dimension(Ns)                          :: Nup,Ndw
    real(8),dimension(Nspin,Nspin,Norb,Norb,Nbath) :: Hbath_tmp
    !
    if(.not.Hsector%status)stop "directMatVec_cc ERROR: Hsector NOT allocated"
    isector=Hsector%index
    !
    if(Nloc/=getdim(isector))stop "directMatVec_cc ERROR: Nloc != dim(isector)"
    !
    !Get diagonal hybridization, bath energy
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
          Hbath_tmp(:,:,:,:,ibath) = bath_from_sym(dmft_bath%item(ibath)%lambda)
          do ispin=1,Nspin
             do iorb=1,Norb
                diag_hybr(ispin,iorb,ibath)=dmft_bath%item(ibath)%v!(ispin)
                bath_diag(ispin,iorb,ibath)=Hbath_tmp(ispin,ispin,iorb,iorb,ibath)
             enddo
          enddo
       enddo
    end select

    !
    Hv=0d0
    !
    !-----------------------------------------------!
    !LOCAL HAMILTONIAN PART: H_loc*vin = vout
    include "direct/Orbs/HxV_local.f90" 
    !
    !UP HAMILTONIAN TERMS
    include "direct/Orbs/HxV_up.f90"
    !
    !DW HAMILTONIAN TERMS
    include "direct/Orbs/HxV_dw.f90"
    !
    if(DimPh>1)then
       !PHONON
       include "direct/Orbs/HxV_ph.f90"
       !ELECTRON-PHONON
       include "direct/Orbs/HxV_eph.f90"
    endif
    !-----------------------------------------------!
    !
    deallocate(diag_hybr,bath_diag)
    return
  end subroutine directMatVec_orbs







#ifdef _MPI
  subroutine directMatVec_MPI_main(Nloc,vin,Hv)
    integer                                        :: Nloc,N
    real(8),dimension(Nloc)                        :: Vin
    real(8),dimension(Nloc)                        :: Hv
    real(8),dimension(:),allocatable               :: vt,Hvt
    !
    integer,dimension(2*Ns_Ud)                     :: Indices,Jndices ![2-2*Norb]
    integer,dimension(Ns_Ud,Ns_Orb)                :: Nups,Ndws       ![1,Ns]-[Norb,1+Nbath]
    integer,dimension(Ns)                          :: Nup,Ndw
    real(8),dimension(Nspin,Nspin,Norb,Norb,Nbath) :: Hbath_tmp
    real(8),dimension(Norb,Norb)                   :: g_matrix
    !
    integer                                        :: i_start,i_end,i_start2,i_end2
    !
    if(.not.Hsector%status)stop "directMatVec_cc ERROR: Hsector NOT allocated"
    isector=Hsector%index
    !
    !Get diagonal hybridization, bath energy
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
          Hbath_tmp(:,:,:,:,ibath) = bath_from_sym(dmft_bath%item(ibath)%lambda)
          do ispin=1,Nspin
             do iorb=1,Norb
                diag_hybr(ispin,iorb,ibath)=dmft_bath%item(ibath)%v!(ispin)
                bath_diag(ispin,iorb,ibath)=Hbath_tmp(ispin,ispin,iorb,iorb,ibath)
             enddo
          enddo
       enddo
    end select
    !
    if(MpiComm==MPI_UNDEFINED.OR.MpiComm==Mpi_Comm_Null)&
         stop "directMatVec_MPI_cc ERRROR: MpiComm = MPI_UNDEFINED"
    ! if(.not.MpiStatus)stop "directMatVec_MPI_cc ERROR: MpiStatus = F"
    !
    Hv=0d0
    !
    !-----------------------------------------------!
    !LOCAL HAMILTONIAN PART: H_loc*vin = vout
    include "direct_mpi/HxV_local.f90"
    !
    !UP HAMILTONIAN TERMS: MEMORY CONTIGUOUS
    include "direct_mpi/HxV_up.f90"
    !
    !DW HAMILTONIAN TERMS: MEMORY NON-CONTIGUOUS
    mpiQup=DimUp/MpiSize
    if(MpiRank<mod(DimUp,MpiSize))MpiQup=MpiQup+1
    do iph=1,DimPh
       allocate(vt(mpiQup*DimDw))
       allocate(Hvt(mpiQup*DimDw))
       vt=0d0
       Hvt=0d0
       i_start = 1 + (iph-1)*DimUp*MpiQdw
       i_end = iph*DimUp*MpiQdw
       !
       call vector_transpose_MPI(DimUp,MpiQdw,Vin(i_start:i_end),DimDw,MpiQup,vt) !Vin^T --> Vt
       include "direct_mpi/HxV_dw.f90"
       deallocate(vt) ; allocate(vt(DimUp*mpiQdw)) ;vt=0d0         !reallocate Vt
       call vector_transpose_MPI(DimDw,mpiQup,Hvt,DimUp,mpiQdw,vt) !Hvt^T --> Vt
       Hv(i_start:i_end) = Hv(i_start:i_end) + Vt
       !
       deallocate(vt,Hvt)
    end do
    !
    if(DimPh>1) then
       !PHONON TERMS
       include "direct_mpi/HxV_ph.f90"
       !ELECTRON-PHONON INTERACTION
       include "direct_mpi/HxV_eph.f90"
    end if
    !
    !NON-LOCAL HAMILTONIAN PART: H_non_loc*vin = vout
    if(Jhflag)then
       N = 0
       call AllReduce_MPI(MpiComm,Nloc,N)
       !
       allocate(vt(N)) ; vt = 0d0
       call allgather_vector_MPI(MpiComm,vin,vt)
       !
       include "direct_mpi/HxV_non_local.f90"
       !
       deallocate(Vt)
    endif
    !-----------------------------------------------!
    !
    deallocate(diag_hybr,bath_diag)
    return
  end subroutine directMatVec_MPI_main



  subroutine directMatVec_MPI_Orbs(Nloc,vin,Hv)
    integer                                        :: Nloc
    real(8),dimension(Nloc)                        :: Vin
    real(8),dimension(Nloc)                        :: Hv
    real(8),dimension(:),allocatable               :: vt,Hvt
    !
    integer,dimension(2*Ns_Ud)                     :: Indices,Jndices ![2-2*Norb]
    integer,dimension(Ns_Ud,Ns_Orb)                :: Nups,Ndws       ![1,Ns]-[Norb,1+Nbath]
    integer,dimension(Ns)                          :: Nup,Ndw
    real(8),dimension(Nspin,Nspin,Norb,Norb,Nbath) :: Hbath_tmp
    integer                                        :: i_start,i_end
    !
    if(.not.Hsector%status)stop "directMatVec_cc ERROR: Hsector NOT allocated"
    isector=Hsector%index
    !
    !Get diagonal hybridization, bath energy
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
          Hbath_tmp(:,:,:,:,ibath) = bath_from_sym(dmft_bath%item(ibath)%lambda)
          do ispin=1,Nspin
             do iorb=1,Norb
                diag_hybr(ispin,iorb,ibath)=dmft_bath%item(ibath)%v!(ispin)
                bath_diag(ispin,iorb,ibath)=Hbath_tmp(ispin,ispin,iorb,iorb,ibath)
             enddo
          enddo
       enddo
    end select
    !
    !    
    if(MpiComm==MPI_UNDEFINED.OR.MpiComm==Mpi_Comm_Null)&
         stop "directMatVec_MPI_cc ERRROR: MpiComm = MPI_UNDEFINED"
    !
    !
    Hv=0d0
    !
    !-----------------------------------------------!
    !LOCAL HAMILTONIAN PART: H_loc*vin = vout
    include "direct_mpi/Orbs/HxV_local.f90"
    !
    !NON-LOCAL TERMS:
    ! 
    !UP HAMILTONIAN TERMS: MEMORY CONTIGUOUS
    include "direct_mpi/Orbs/HxV_up.f90"
    !
    !DW HAMILTONIAN TERMS: MEMORY NON-CONTIGUOUS
    mpiQup=DimUp/MpiSize
    if(MpiRank<mod(DimUp,MpiSize))MpiQup=MpiQup+1
    !
    do iph=1,DimPh
       allocate(vt(mpiQup*DimDw)) ;vt=0d0
       allocate(Hvt(mpiQup*DimDw));Hvt=0d0
       i_start = 1 + (iph-1)*DimUp*MpiQdw
       i_end = iph*DimUp*MpiQdw
       !
       call vector_transpose_MPI(DimUp,MpiQdw,Vin(i_start:i_end),DimDw,MpiQup,vt) !Vin^T --> Vt
       include "direct_mpi/Orbs/HxV_dw.f90"
       deallocate(vt) ; allocate(vt(DimUp*mpiQdw)) ;vt=0d0         !reallocate Vt
       call vector_transpose_MPI(DimDw,mpiQup,Hvt,DimUp,mpiQdw,vt) !Hvt^T --> Vt
       Hv(i_start:i_end) = Hv(i_start:i_end) + Vt
       !
       deallocate(Hvt,vt)
    enddo
    !
    if(DimPh>1)then
       !PHONON TERMS
       include "direct_mpi/Orbs/HxV_ph.f90"
       !ELECTRON-PHONON INTERACTION
       include "direct_mpi/Orbs/HxV_eph.f90"
    endif
    !-----------------------------------------------!
    !
    deallocate(diag_hybr,bath_diag)
    return
  end subroutine directMatVec_MPI_Orbs
#endif


END MODULE ED_HAMILTONIAN_DIRECT_HXV
