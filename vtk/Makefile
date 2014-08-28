F90 = ifort
OPT = -openmp -i8 -xHOST

smash : vtk-generator.F90
	$(F90) -o vtk-generator $(OPT) vtk-generator.F90
