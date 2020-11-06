MODULE ED_BATH_FUNCTIONS
  USE SF_CONSTANTS, only: zero
  USE SF_IOTOOLS, only:free_unit,reg,file_length,txtfy
  USE SF_LINALG, only: eye,inv,zeye,inv_her
  USE ED_INPUT_VARS
  USE ED_VARS_GLOBAL
  USE ED_BATH
  USE ED_AUX_FUNX
  implicit none

  private


  !
  !\DELTA, THE HYBRIDIZATION FUNCTION
  interface delta_bath_function
     module procedure delta_bath_array
  end interface delta_bath_function
  !
  !NON-INTERACTING GREEN'S FUNCTION 
  interface g0and_bath_function
     module procedure g0and_bath_array
  end interface g0and_bath_function
  !
  !INVERSE NON-INTERACTING GREEN'S FUNCTION 
  interface invg0_bath_function
     module procedure invg0_bath_array
  end interface invg0_bath_function

  public :: delta_bath_function
  public :: g0and_bath_function
  public :: invg0_bath_function






contains



  function delta_bath_array(x,dmft_bath_) result(Delta)
    complex(8),dimension(:),intent(in)                  :: x
    type(effective_bath)                                :: dmft_bath_
    complex(8),dimension(Nspin,Nspin,Norb,Norb,size(x)) :: Delta
    integer                                             :: i,ih,L
    integer                                             :: iorb,jorb,ispin,jspin,ibath
    integer                                             :: io,jo
    real(8),dimension(Nbath)                            :: eps,dps,vps
    real(8),dimension(Norb,Nbath)                       :: vops
    !
    real(8),dimension(Nspin,Nbath)                      :: ehel
    real(8),dimension(Nspin,Nspin,Nbath)                :: whel
    real(8),dimension(Nspin,Nspin,Norb,Nbath)           :: wohel
    !
    complex(8),dimension(Nspin*Norb,Nspin*Norb)         :: invH_k
    complex(8),dimension(Nspin,Nspin,Norb,Norb)         :: invH_knn
    !
    !
    Delta=zero
    !
    L = size(x)
    !
    select case(bath_type)
    case default                !normal: only _{aa} are allowed (no inter-orbital local mixing)
       !
       !\Delta_{aa} = \sum_k [ V_{a}(k) * V_{a}(k)/(iw_n - E_{a}(k)) ]
       do ispin=1,Nspin
          do iorb=1,Norb
             eps = dmft_bath_%e(ispin,iorb,1:Nbath)
             vps = dmft_bath_%v(ispin,iorb,1:Nbath)
             do i=1,L
                Delta(ispin,ispin,iorb,iorb,i) = sum( vps(:)*vps(:)/(x(i) - eps(:)) )
             enddo
          enddo
       enddo
       !
       !
    case ("hybrid")             !hybrid: all _{ab} components allowed (inter-orbital local mixing present)
       !
       !\Delta_{ab} = \sum_k [ V_{a}(k) * V_{b}(k)/(iw_n - E(k)) ]
       do ispin=1,Nspin
          eps  = dmft_bath_%e(ispin,1     ,1:Nbath)
          vops = dmft_bath_%v(ispin,1:Norb,1:Nbath)
          do iorb=1,Norb
             do jorb=1,Norb
                do i=1,L
                   Delta(ispin,ispin,iorb,jorb,i) = sum( vops(iorb,:)*vops(jorb,:)/(x(i) - eps(:)) )
                enddo
             enddo
          enddo
       enddo
       !
    case ("replica")
       !
       invH_k=zero
       do i=1,L
          do ibath=1,Nbath
             invH_knn = bath_from_sym(dmft_bath_%item(ibath)%lambda)
             invH_k   = nn2so_reshape(invH_knn,Nspin,Norb)
             invH_k   = zeye(Nspin*Norb)*x(i) - invH_k
             call inv(invH_k)
             invH_knn = so2nn_reshape(invH_k,Nspin,Norb)
             do ispin=1,Nspin   !this must be diagonal in the spin channel
                Delta(ispin,ispin,:,:,i)=Delta(ispin,ispin,:,:,i) + &
                     (dmft_bath_%item(ibath)%v(ispin)**2)*invH_knn(ispin,ispin,:,:)
             enddo
          enddo
          !
       enddo
       !
    end select
  end function delta_bath_array






  function g0and_bath_array(x,dmft_bath_) result(G0and)
    complex(8),dimension(:),intent(in)                  :: x
    type(effective_bath)                                :: dmft_bath_
    complex(8),dimension(Nspin,Nspin,Norb,Norb,size(x)) :: G0and,Delta
    integer                                             :: iorb,jorb,ispin,jspin,io,jo,Nso,i,L
    real(8),dimension(size(x))                          :: det
    complex(8),dimension(size(x))                       :: fg,ff
    complex(8),dimension(Norb,Norb)                     :: fgorb
    !
    G0and = zero
    !
    L=size(x)
    !
    Delta = delta_bath_array(x,dmft_bath_)
    !
    select case(bath_type)
    case default                !normal: only _{aa} are allowed (no inter-orbital local mixing)
       do ispin=1,Nspin
          do iorb=1,Norb
             fg(:) = x(:) + xmu - impHloc(ispin,ispin,iorb,iorb) - Delta(ispin,ispin,iorb,iorb,:)
             G0and(ispin,ispin,iorb,iorb,:) = one/fg(:)
          enddo
       enddo
       !
       !
    case ("hybrid","replica")             !hybrid: all _{ab} components allowed (inter-orbital local mixing present)
       do ispin=1,Nspin         !Spin diagonal
          do i=1,L
             fgorb = (x(i)+xmu)*zeye(Norb) - impHloc(ispin,ispin,:,:) - Delta(ispin,ispin,:,:,i)
             call inv(fgorb)
             G0and(ispin,ispin,:,:,i)=fgorb
          enddo
       enddo
       !
       !
    end select
  end function g0and_bath_array





  function invg0_bath_array(x,dmft_bath_) result(G0and)
    complex(8),dimension(:),intent(in)                  :: x
    type(effective_bath)                                :: dmft_bath_
    complex(8),dimension(Nspin,Nspin,Norb,Norb,size(x)) :: G0and,Delta
    integer                                             :: i,iorb,jorb,ispin,jspin,io,jo,Nso,L
    !
    G0and = zero
    !
    L=size(x)
    !
    select case(bath_type)
    case default                !normal: only _{aa} are allowed (no inter-orbital local mixing)
       !
       Delta = delta_bath_array(x,dmft_bath_)
       do ispin=1,Nspin
          do iorb=1,Norb
             G0and(ispin,ispin,iorb,iorb,:) = x(:) + xmu - impHloc(ispin,ispin,iorb,iorb) - Delta(ispin,ispin,iorb,iorb,:)
          enddo
       enddo
       !
       !
    case ("hybrid","replica")             !hybrid: all _{ab} components allowed (inter-orbital local mixing present)
       !
       Delta = delta_bath_array(x,dmft_bath_)
       do ispin=1,Nspin
          do i=1,L
             G0and(ispin,ispin,:,:,i) = (x(i)+xmu)*zeye(Norb)-impHloc(ispin,ispin,:,:)-Delta(ispin,ispin,:,:,i)
          enddo
       enddo
       !
    end select
    !
  end function invg0_bath_array







END MODULE ED_BATH_FUNCTIONS
