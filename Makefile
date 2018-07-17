

FC=gfortran
OPT= -ffree-line-length-none -O0 -g -fbacktrace

LIB=-L/home/joe/spack/opt/spack/linux-ubuntu16.04-x86_64/gcc-5.4.0/netcdf-4.6.1-ohldsascx6453mdolilxykvhnufhl3q3/lib -lnetcdf -L/home/joe/spack/opt/spack/linux-ubuntu16.04-x86_64/gcc-5.4.0/netcdf-fortran-4.4.4-llqdlrjpkkf7eatofynzx3djwbvl2ubt/lib -lnetcdff

INC=-I/home/joe/spack/opt/spack/linux-ubuntu16.04-x86_64/gcc-5.4.0/netcdf-fortran-4.4.4-qbxkckbzlpyq7xadda3z3mbbw4jkfs22/include

regrid: WAM_Reinterpolate.f90
	${FC} ${OPT} WAM_Reinterpolate.f90 ${LIB} ${INC} -o $@
