rm int1
gfortran -std=legacy -c ../../src/kernels_2d/lap_kernels_2d.f90 -o ../../src/kernels_2d/lap_kernels_2d.o
gfortran -std=legacy test_lap2d_gdn.f -o ./int1 ../../src/kernels_2d/lap_kernels_2d.o
./int1
