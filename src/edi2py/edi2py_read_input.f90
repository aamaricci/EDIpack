!ED_MAIN:
subroutine ed_read_input_c(instr) bind(c, name='read_input')
    character(kind=c_char), dimension(*), intent(IN) :: instr
    character(len=20), allocatable :: INPUTunit
    integer :: length
    integer :: i
    length=0
    INPUTunit=" "
    do
       if (instr(length+1) == C_NULL_CHAR) exit
       length = length + 1
    end do
    do i = 1, length
      INPUTunit(i:i) = instr(i)
    enddo
    INPUTunit=trim(INPUTunit)
    call ed_read_input(INPUTunit) 
end subroutine ed_read_input_c
