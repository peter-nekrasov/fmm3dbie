gfortran -std=legacy -c ../../src/kernels/lap_kernels.f90 -o ../../src/kernels/lap_kernels.o
gfortran -std=legacy test_lap2d_sprime.f -o ./int1 ../../src/kernels/lap_kernels.o
./int1
