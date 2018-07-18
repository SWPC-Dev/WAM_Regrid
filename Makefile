

FC=gfortran
OPT= -ffree-line-length-none -O3

LIB=-L/home/joe/spack/opt/spack/linux-ubuntu16.04-x86_64/gcc-5.4.0/netcdf-fortran-4.4.4-llqdlrjpkkf7eatofynzx3djwbvl2ubt/lib -lnetcdff 
LIB+= -L/home/joe/spack/opt/spack/linux-ubuntu16.04-x86_64/gcc-5.4.0/netcdf-4.6.1-ohldsascx6453mdolilxykvhnufhl3q3/lib -lnetcdf

INC=-I/home/joe/spack/opt/spack/linux-ubuntu16.04-x86_64/gcc-5.4.0/netcdf-fortran-4.4.4-qbxkckbzlpyq7xadda3z3mbbw4jkfs22/include

.PHONY : clean \
         regrid

regrid: WAM_Reinterpolate.o
	${FC} ${OPT} WAM_Reinterpolate.o ${LIB} ${INC} -o $@

WAM_Reinterpolate.o: WAM_Reinterpolate.f90
	${FC} ${OPT} -c WAM_Reinterpolate.f90 ${LIB} ${INC} -o $@


clean: 
	rm *.o *.mod regrid
