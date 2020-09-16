FC=gfortran
FFLAGS=-I`nf-config --fflags`
FLIBS=-L/usr/local/Cellar/netcdf/4.7.4_1/lib -lnetcdff -lnetcdf


OBJS = configMod.o fields.o gridMod.o kinds.o ncdf_wrapper.o netCDFio.o

filterQGCM: $(OBJS)
	$(FC) $(FLIBS) $(OBJS) $@.F90 -o $@ 

# mpiMod.o: 
# 	$(FC) -c mpiMod.F90

netCDFio.o: kinds.o gridMod.o fields.o ncdf_wrapper.o
	$(FC) -c netCDFio.F90

ncdf_wrapper.o : ncdf_wrapper.F90
	$(FC) $(FFLAGS) -c ncdf_wrapper.F90

kinds.o: kinds.F90
	$(FC) -c kinds.F90

gridMod.o: gridMod.F90 kinds.o configMod.o
	$(FC) -c gridMod.F90

configMod.o: configMod.F90 kinds.o
	$(FC) -c configMod.F90

fields.o: kinds.o gridMod.o configMod.o 
	$(FC) -c fields.F90

clean:
	rm -rf *.o *.mod filterQGCM

tidy:
	rm -rf *.mod *.o 