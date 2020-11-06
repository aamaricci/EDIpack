!+-------------------------------------------------------------------+
!PURPOSE  : Deallocate the ED bath
!+-------------------------------------------------------------------+
subroutine deallocate_dmft_bath(dmft_bath_)
  type(effective_bath) :: dmft_bath_
  integer              :: ibath,isym
  if(.not.dmft_bath_%status)return
  if(allocated(dmft_bath_%e))   deallocate(dmft_bath_%e)
  if(allocated(dmft_bath_%v))   deallocate(dmft_bath_%v)
  if(bath_type=="replica")then
     do ibath=1,Nbath
        dmft_bath_%item(ibath)%N_dec= 0
        deallocate(dmft_bath_%item(ibath)%v)
        deallocate(dmft_bath_%item(ibath)%lambda)
     enddo
     deallocate(dmft_bath_%item)
  endif
  dmft_bath_%status=.false.
end subroutine deallocate_dmft_bath



!+-------------------------------------------------------------------+
!PURPOSE  : Allocate the ED bath
!+-------------------------------------------------------------------+
subroutine allocate_dmft_bath(dmft_bath_)
  type(effective_bath) :: dmft_bath_
  integer              :: Nsym,ibath
  if(dmft_bath_%status)call deallocate_dmft_bath(dmft_bath_)
  !
  select case(bath_type)
  case default
     !
     allocate(dmft_bath_%e(Nspin,Norb,Nbath))  !local energies of the bath
     allocate(dmft_bath_%v(Nspin,Norb,Nbath))  !same-spin hybridization 
     !
  case('hybrid')
     !
     allocate(dmft_bath_%e(Nspin,1,Nbath))     !local energies of the bath
     allocate(dmft_bath_%v(Nspin,Norb,Nbath))  !same-spin hybridization 
     !
  case('replica')
     !
     if(.not.allocated(lambda_impHloc))stop "lambda_impHloc not allocated in allocate_dmft_bath" !FIXME
     call deallocate_dmft_bath(dmft_bath_)     !
     allocate(dmft_bath_%item(Nbath))
     !
     Nsym=size(lambda_impHloc)
     !
     !ALLOCATE coefficients vectors
     do ibath=1,Nbath
        dmft_Bath_%item(ibath)%N_dec=Nsym
        allocate(dmft_bath_%item(ibath)%v(Nspin))
        allocate(dmft_bath_%item(ibath)%lambda(Nsym))
     enddo
     !
  end select
  !
  dmft_bath_%status=.true.
  !
end subroutine allocate_dmft_bath




!+-------------------------------------------------------------------+
!PURPOSE  : Reconstruct bath matrix from lambda vector
!+-------------------------------------------------------------------+
function bath_from_sym(lambdavec) result (Hbath)
  integer                                  :: Nsym,isym
  real(8),dimension(:)                     :: lambdavec
  real(8),dimension(Nspin,Nspin,Norb,Norb) :: Hbath
  !
  Nsym=size(lambdavec)
  !
  Hbath=zero
  !
  do isym=1,Nsym
     Hbath=Hbath+lambdavec(isym)*H_Basis(isym)%O
  enddo
  !
end function bath_from_sym




!+------------------------------------------------------------------+
!PURPOSE  : Initialize the DMFT loop, builindg H parameters and/or 
!reading previous (converged) solution
!+------------------------------------------------------------------+
subroutine init_dmft_bath(dmft_bath_)
  type(effective_bath) :: dmft_bath_
  integer,dimension(Nbath) :: Nlambdas
  integer              :: i,unit,flen,Nh,isym,Nsym
  integer              :: io,jo,iorb,ispin,jorb,jspin,ibath
  logical              :: IOfile
  real(8)              :: de
  real(8)              :: rescale(Nbath),offset_b(Nbath)
  !
  if(.not.dmft_bath_%status)stop "init_dmft_bath error: bath not allocated"
  !
  if(Nbath>1)then
     rescale=linspace(HWBAND/Nbath,HWBAND,Nbath)
  else
     rescale(1)=0.d0
  endif
  !
  select case(bath_type)
  case default
     !Get energies:
     dmft_bath_%e(:,:,1)    =-hwband
     dmft_bath_%e(:,:,Nbath)= hwband
     Nh=Nbath/2
     if(mod(Nbath,2)==0.and.Nbath>=4)then
        de=hwband/max(Nh-1,1)
        dmft_bath_%e(:,:,Nh)  = -1.d-1
        dmft_bath_%e(:,:,Nh+1)=  1.d-1
        do i=2,Nh-1
           dmft_bath_%e(:,:,i)   =-hwband + (i-1)*de
           dmft_bath_%e(:,:,Nbath-i+1)= hwband - (i-1)*de
        enddo
     elseif(mod(Nbath,2)/=0.and.Nbath>=3)then
        de=hwband/Nh
        dmft_bath_%e(:,:,Nh+1)= 0d0
        do i=2,Nh
           dmft_bath_%e(:,:,i)        =-hwband + (i-1)*de
           dmft_bath_%e(:,:,Nbath-i+1)= hwband - (i-1)*de
        enddo
     endif
     !Get spin-keep yhbridizations
     do i=1,Nbath
        dmft_bath_%v(:,:,i)=max(0.1d0,1d0/sqrt(dble(Nbath)))
     enddo
     !
  case('replica')
     !BATH V INITIALIZATION
     do ibath=1,Nbath
        do ispin=1,Nspin
           dmft_bath%item(ibath)%v(ispin)=max(0.1d0,1d0/sqrt(dble(Nbath)))
        enddo
     enddo
     !
     !BATH LAMBDAS INITIALIZATION
     do ibath=1,Nbath
        Nsym = dmft_bath%item(ibath)%N_dec
        do isym=1,Nsym
           if(is_diagonal(H_basis(isym)%O))then
              dmft_bath%item(ibath)%lambda(isym)=rescale(ibath)*lambda_impHloc(isym)
           else
              dmft_bath%item(ibath)%lambda(isym) =  lambda_impHloc(isym)
           endif
        enddo
     enddo
     !
  end select
  !
  !
  !
  !Read from file if exist:
  !
  inquire(file=trim(Hfile)//trim(ed_file_suffix)//".restart",exist=IOfile)
  if(IOfile)then
     write(LOGfile,"(A)")'Reading bath from file'//trim(Hfile)//trim(ed_file_suffix)//".restart"
     unit = free_unit()
     flen = file_length(trim(Hfile)//trim(ed_file_suffix)//".restart")
     !
     open(unit,file=trim(Hfile)//trim(ed_file_suffix)//".restart")
     !
     select case(bath_type)
     case default
        !
        read(unit,*)
        do i=1,min(flen,Nbath)
           read(unit,*)((&
                dmft_bath_%e(ispin,iorb,i),&
                dmft_bath_%v(ispin,iorb,i),&
                iorb=1,Norb),ispin=1,Nspin)
        enddo
        !
     case ('hybrid')
        read(unit,*)
        !
        do i=1,min(flen,Nbath)
           read(unit,*)(&
                dmft_bath_%e(ispin,1,i),&
                (&
                dmft_bath_%v(ispin,iorb,i),&
                iorb=1,Norb),&
                ispin=1,Nspin)
        enddo
        !
     case ('replica')
        !
        !read number of lambdas
        do ibath=1,Nbath
           read(unit,"(I3)")Nlambdas(ibath)
        enddo
        do ibath=1,Nbath
           !read V
           do ispin=1,Nspin
              read(unit,*)dmft_bath%item(ibath)%v(ispin)
           enddo
           !read lambdas
           read(unit,*)(dmft_bath%item(ibath)%lambda(jo),jo=1,Nlambdas(ibath))
        enddo
        !
        !
     end select
     close(unit)
  endif
end subroutine init_dmft_bath



!+-------------------------------------------------------------------+
!PURPOSE  : write out the bath to a given unit with 
! the following column formatting: 
! [(Ek_iorb,Vk_iorb)_iorb=1,Norb]_ispin=1,Nspin
!+-------------------------------------------------------------------+
subroutine write_dmft_bath(dmft_bath_,unit)
  type(effective_bath) :: dmft_bath_
  integer,optional     :: unit
  integer              :: unit_
  integer              :: i
  integer              :: io,jo,iorb,ispin,isym
  real(8)              :: hybr_aux
  real(8)              :: hrep_aux(Nspin*Norb,Nspin*Norb)
  !
  character(len=64)    :: string_fmt,string_fmt_first
  unit_=LOGfile;if(present(unit))unit_=unit
  if(.not.dmft_bath_%status)stop "write_dmft_bath error: bath not allocated"
  select case(bath_type)
  case default
     !
     write(unit_,"(90(A21,1X))")&
          ((&
          "#Ek_l"//reg(txtfy(iorb))//"_s"//reg(txtfy(ispin)),&
          "Vk_l"//reg(txtfy(iorb))//"_s"//reg(txtfy(ispin)),&
          iorb=1,Norb),ispin=1,Nspin)
     do i=1,Nbath
        write(unit_,"(90(ES21.12,1X))")((&
             dmft_bath_%e(ispin,iorb,i),&
             dmft_bath_%v(ispin,iorb,i),&
             iorb=1,Norb),ispin=1,Nspin)
     enddo
     !
  case('hybrid')
     !
     write(unit_,"(90(A21,1X))")(&
          "#Ek_s"//reg(txtfy(ispin)),&
          ("Vk_l"//reg(txtfy(iorb))//"_s"//reg(txtfy(ispin)),iorb=1,Norb),&
          ispin=1,Nspin)
     do i=1,Nbath
        write(unit_,"(90(ES21.12,1X))")(&
             dmft_bath_%e(ispin,1,i),&
             (dmft_bath_%v(ispin,iorb,i),iorb=1,Norb),&
             ispin=1,Nspin)
     enddo
     !
  case ('replica')
     !
     string_fmt      ="(A8,A5,"//str(Nspin*Norb)//"(ES8.4,1X))"
     !
     if(unit_==LOGfile)then
        if(Nspin*Norb.le.8)then
           write(unit_,"(A1)")" "
           write(unit_,"(A8,A3,a5,90(A15,1X))")"V","||"," ","Re(H) | Im(H)"
           do ibath=1,Nbath
              write(unit_,"(A1)")" "
              Hrep_aux   = nn2so_reshape( bath_from_sym(dmft_bath_%item(ibath)%lambda) ,Nspin,Norb)
              do ispin=1,Nspin
                 write(unit_,"(ES8.4)")dmft_bath_%item(ibath)%v(ispin)
              enddo
              do io=1,Nspin*Norb
                 write(unit_,string_fmt) "  "  ,"||  ",(hrep_aux(io,jo),jo=1,Nspin*Norb)
              enddo
           enddo
           write(unit_,"(A1)")" "
        else
           write(LOGfile,"(A)")"Bath matrix too large to print: printing the parameters (including eventual offset)."
           write(unit_,"(A8,A5,90(A8,1X))")"V"," ","lambdas"        
           do ibath=1,Nbath
              do ispin=1,Nspin
                 write(unit_,"(ES8.4)")dmft_bath_%item(ibath)%v(ispin)
              enddo
              write(unit_,"(A8,A5,90(ES8.4,1X))")"","|   ",&
                   (dmft_bath_%item(ibath)%lambda(io),io=1,dmft_bath_%item(ibath)%N_dec)
           enddo
        endif
     else
        do ibath=1,Nbath
           write(unit,"(I3)")dmft_bath_%item(ibath)%N_dec
        enddo
        do ibath=1,Nbath
           do ispin=1,Nspin
              write(unit,*)dmft_bath_%item(ibath)%v(ispin)
           enddo
           write(unit,*)(dmft_bath_%item(ibath)%lambda(jo),jo=1,dmft_bath_%item(ibath)%N_dec)
        enddo
     endif
     !
  end select
end subroutine write_dmft_bath






!+-------------------------------------------------------------------+
!PURPOSE  : save the bath to a given file using the write bath
! procedure and formatting: 
!+-------------------------------------------------------------------+
subroutine save_dmft_bath(dmft_bath_,file,used)
  type(effective_bath)      :: dmft_bath_
  character(len=*),optional :: file
  character(len=256)        :: file_
  logical,optional          :: used
  logical                   :: used_
  character(len=16)         :: extension
  integer                   :: unit_
  if(.not.dmft_bath_%status)stop "save_dmft_bath error: bath is not allocated"
  used_=.false.;if(present(used))used_=used
  extension=".restart";if(used_)extension=".used"
  file_=str(str(Hfile)//str(ed_file_suffix)//str(extension))
  if(present(file))file_=str(file)
  unit_=free_unit()
  open(unit_,file=str(file_))
  call write_dmft_bath(dmft_bath_,unit_)
  close(unit_)
end subroutine save_dmft_bath




!+-------------------------------------------------------------------+
!PURPOSE  : set the bath components from a given user provided 
! bath-array 
!+-------------------------------------------------------------------+
subroutine set_dmft_bath(bath_,dmft_bath_)
  real(8),dimension(:)   :: bath_
  type(effective_bath)   :: dmft_bath_
  integer                :: stride,io,jo,i
  integer                :: iorb,ispin,jorb,jspin,ibath
  logical                :: check
  !
  if(.not.dmft_bath_%status)stop "set_dmft_bath error: bath not allocated"
  check = check_bath_dimension(bath_)
  if(.not.check)stop "set_dmft_bath error: wrong bath dimensions"
  !
  select case(bath_type)
  case default
     !
     stride = 0
     do ispin=1,Nspin
        do iorb=1,Norb
           do i=1,Nbath
              io = stride + i + (iorb-1)*Nbath + (ispin-1)*Nbath*Norb
              dmft_bath_%e(ispin,iorb,i) = bath_(io)
           enddo
        enddo
     enddo
     stride = Nspin*Norb*Nbath
     do ispin=1,Nspin
        do iorb=1,Norb
           do i=1,Nbath
              io = stride + i + (iorb-1)*Nbath + (ispin-1)*Nbath*Norb
              dmft_bath_%v(ispin,iorb,i) = bath_(io)
           enddo
        enddo
     enddo
     !
     !
  case ('hybrid')
     !
     stride = 0
     do ispin=1,Nspin
        do i=1,Nbath
           io = stride + i + (ispin-1)*Nbath
           dmft_bath_%e(ispin,1,i) = bath_(io)
        enddo
     enddo
     stride = Nspin*Nbath
     do ispin=1,Nspin
        do iorb=1,Norb
           do i=1,Nbath
              io = stride + i + (iorb-1)*Nbath + (ispin-1)*Norb*Nbath
              dmft_bath_%v(ispin,iorb,i) = bath_(io)
           enddo
        enddo
     enddo
     !
     !
  case ('replica')
     !
     ! do ibath=1,Nbath
     !    dmft_bath_%item(ibath)%N_dec = 0
     !    dmft_bath_%item(ibath)%v     = 0d0
     !    dmft_bath_%item(ibath)%lambda= 0d0
     ! enddo
     ! !
     stride = 0
     !Get N_dec
     do ibath=1,Nbath
        stride = stride + 1
        dmft_bath_%item(ibath)%N_dec=NINT(bath_(stride))
     enddo
     !Get N_dec
     !get Lambdas
     do ibath=1,Nbath
        do ispin=1,Nspin
           stride = stride + 1
           dmft_bath_%item(ibath)%v(ispin) = bath_(stride)
        enddo
        dmft_bath_%item(ibath)%lambda=bath_(stride+1 :stride+dmft_bath_%item(ibath)%N_dec)
        stride=stride+dmft_bath_%item(ibath)%N_dec
     enddo
  end select
end subroutine set_dmft_bath



!+-------------------------------------------------------------------+
!PURPOSE  : copy the bath components back to a 1-dim array 
!+-------------------------------------------------------------------+
subroutine get_dmft_bath(dmft_bath_,bath_)
  type(effective_bath)   :: dmft_bath_
  real(8),dimension(:)   :: bath_
  real(8)                :: hrep_aux(Nspin*Norb,Nspin*Norb)
  integer                :: stride,io,jo,i
  integer                :: iorb,ispin,jorb,jspin,ibath,maxspin
  logical                :: check
  if(.not.dmft_bath_%status)stop "get_dmft_bath error: bath not allocated"
  check=check_bath_dimension(bath_)
  if(.not.check)stop "get_dmft_bath error: wrong bath dimensions"
  !
  select case(bath_type)
  case default
     !
     stride = 0
     do ispin=1,Nspin
        do iorb=1,Norb
           do i=1,Nbath
              io = stride + i + (iorb-1)*Nbath + (ispin-1)*Nbath*Norb
              bath_(io) = dmft_bath_%e(ispin,iorb,i) 
           enddo
        enddo
     enddo
     stride = Nspin*Norb*Nbath
     do ispin=1,Nspin
        do iorb=1,Norb
           do i=1,Nbath
              io = stride + i + (iorb-1)*Nbath + (ispin-1)*Nbath*Norb
              bath_(io) = dmft_bath_%v(ispin,iorb,i)
           enddo
        enddo
     enddo
     !
     !
  case ('hybrid')
     !
     stride = 0
     do ispin=1,Nspin
        do i=1,Nbath
           io = stride + i + (ispin-1)*Nbath
           bath_(io) =  dmft_bath_%e(ispin,1,i)
        enddo
     enddo
     stride = Nspin*Nbath
     do ispin=1,Nspin
        do iorb=1,Norb
           do i=1,Nbath
              io = stride + i + (iorb-1)*Nbath + (ispin-1)*Norb*Nbath
              bath_(io) =  dmft_bath_%v(ispin,iorb,i)
           enddo
        enddo
     enddo
     !
     !
  case ('replica')
     !
     stride = 0
     do ibath=1,Nbath
        stride = stride + 1
        bath_(stride)=dmft_bath_%item(ibath)%N_dec
     enddo
     do ibath=1,Nbath
        do ispin=1,Nspin
           stride = stride + 1
           bath_(stride)=dmft_bath_%item(ibath)%v(ispin)
        enddo
        bath_(stride+1 : stride+dmft_bath_%item(ibath)%N_dec)=dmft_bath_%item(ibath)%lambda
        stride=stride+dmft_bath_%item(ibath)%N_dec
     enddo
  end select
end subroutine get_dmft_bath



