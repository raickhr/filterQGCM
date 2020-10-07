FC=mpif90
FFLAGS=`nf-config --fflags`
FLIBS=`nf-config --flibs`


OBJS = configMod.o fields.o gridMod.o kinds.o ncdf_wrapper.o netCDFio.o filter.o mpiMod.o workDiv.o operators.o constants.o

filterQGCM: $(OBJS)
	$(FC) $(FLIBS) $(OBJS) $@.F90 -o $@ 
	
kinds.o: kinds.F90
	$(FC) -c kinds.F90

constants.o: kinds.o
	$(FC) -c constants.F90

mpiMod.o: mpiMod.F90 kinds.o
	$(FC) -c mpiMod.F90

filter.o: filter.F90 fields.o kinds.o gridMod.o workDiv.o
	$(FC) -c filter.F90

netCDFio.o: kinds.o gridMod.o fields.o ncdf_wrapper.o
	$(FC) -c netCDFio.F90

ncdf_wrapper.o : ncdf_wrapper.F90
	$(FC) $(FFLAGS) -c ncdf_wrapper.F90

gridMod.o: gridMod.F90 kinds.o configMod.o mpiMod.o
	$(FC) -c gridMod.F90

configMod.o: configMod.F90 kinds.o constants.o
	$(FC) -c configMod.F90

fields.o: kinds.o gridMod.o configMod.o 
	$(FC) -c fields.F90

workDiv.o: gridMod.o mpiMod.o
	$(FC) -c workDiv.F90

operators.o: kinds.o constants.o gridMod.o
	$(FC) -c operators.F90

clean:
	rm -rf *.o *.mod filterQGCM

tidy:
	rm -rf *.mod *.o 