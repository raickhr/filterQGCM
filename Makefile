FC=gfortran
FFLAGS=-L`nf-config --fflags`


OBJS = kinds.o gridMod.o configMod.o

filterQGCM: $(OBJS)
	$(FC) $(OBJS) $@.F90 -o $@ 

kinds.o: kinds.F90
	$(FC) -c kinds.F90

gridMod.o: gridMod.F90 kinds.o configMod.o
	$(FC) -c gridMod.F90

configMod.o: configMod.F90 kinds.o
	$(FC) -c configMod.F90

clean:
	rm -rf *.o *.mod filterQGCM

tidy:
	rm -rf *.mod *.o 