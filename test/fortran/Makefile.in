FC=mpif90

INC=$(shell pkg-config --cflags scifor edipack)
LIB=$(shell pkg-config --libs   scifor edipack)

FFLAG = -O2 -ffree-line-length-none -cpp -D_MPI

# Extends the implicit support of the Makefile to .f90 files
.SUFFIXES: .f90

all:
	@echo ""	
	$(FC) $(FFLAG) $(INC) $(LIB)  hm_bethe.f90 -o hm_bethe.exe

clean: 
	@echo "Cleaning:"
	@rm -f *.exe *.mod *.o *~
