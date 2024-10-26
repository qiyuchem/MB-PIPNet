#compiles pes_fit program
FC = ifort

FFLAGS = -r8 -O3 -mkl=sequential

OBJ = constants.o bemsa1b.o bemsa2b.o MolNN_shell.o pes_shell.o getpot.o

%.o : %.f90
	$(FC) -c $(FFLAGS) $<


getpot.x : $(OBJ)
	$(FC) $^ $(FFLAGS) $(LIBS) -o $@ 

clean: 
	rm *.o *.mod
