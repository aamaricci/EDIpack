!ED_AUX_FUNX:
subroutine search_variable(var,ntmp,converged) bind(c, name='search_variable')
  real(c_double),dimension(1)         :: var(1)
  real(c_double),dimension(1)         :: ntmp(1)
  integer(c_int),dimension(1)         :: converged(1)
  logical                             :: bool
  converged(1)=0
  call ed_search_variable(var(1),ntmp(1),bool)
  if (bool) converged(1)=1
end subroutine search_variable




!AUX:
subroutine check_convergence(Xnew,dim_xnew,eps,N1,N2,oerr,convergence) bind(c, name='check_convergence')
  complex(c_double_complex)          :: Xnew(dim_xnew)
  integer(c_int),value               :: dim_xnew
  real(c_double),value               :: eps
  integer(c_int),value               :: N1,N2
  real(c_double),dimension(1)        :: oerr(1)
  integer(c_int),dimension(1)        :: convergence(1)  
  integer                            :: i,Msum
  real(c_double)                     :: err
  real(c_double)                     :: M,S
  complex(8),save,allocatable        :: Xold(:)
  integer,save                       :: success=0,check=1
  character(len=100)                 :: file_
  file_='error.err'
  Msum=size(Xnew)
  if(.not.allocated(Xold))then
     allocate(Xold(Msum))
     Xold=0.d0
  endif
  S=0.d0 ; M=0.d0
  do i=1,Msum
     M=M + abs(Xnew(i)-Xold(i))
     S=S + abs(Xnew(i))
  enddo
  err= M/S
  Xold=Xnew
  open(10,file=reg(file_),position="append")
  write(10,"(I5,ES15.7)")check,err
  close(10)
  if(err < eps)then
     success=success+1
  else
     success=0
  endif
  convergence(1)=0
  if(success > N1)convergence(1)=1
  if(check>=N2)then
     open(10,file="ERROR.README")
     write(10,*)""
     close(10)
     write(*,"(A,I4,A)")"Not converged after",N2," iterations."
  endif
  if(convergence(1)==1)then
     write(*,"(A,ES15.7)")bold_green("error="),err
  else
     if(err < eps)then
        write(*,"(A,ES15.7)")bold_yellow("error="),err
     else
        write(*,"(A,ES15.7)")bold_red("error="),err
     endif
  endif
  oerr(1)=err
  check=check+1
end subroutine check_convergence





