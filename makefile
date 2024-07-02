#Makefile for SX-9 or x86_64 @ ifort
FC=ifx
FFLAGS= -O2 -mavx2 -qopenmp -inline-level=2 -qmkl=sequential
OBJ= main.o chempot.o mkhm.o

.SUFFIXES:
main: $(OBJ)
	$(FC) $(OBJ) $(FFLAGS) -o main $(FFLIBS)
%.o: %.f90
	$(FC) $(FFLAGS) -c $<
.PHONY : clean bury new
comclean = *.o *.mod main *~
comnew = *.eps *.gnu *.dat fort.* *.optrpt *.L out*
clean :
	rm -rf $(comclean)
new :
	rm -f $(comnew)
bury :
	rm -rf $(comclean) $(comnew)
