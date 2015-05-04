PROG= vtk-generator cube-generator
F90 = ifort
OPT = -openmp -i8 -xHOST

all : $(PROG)

vtk-generator : vtk-generator.F90
	$(F90) -o vtk-generator $(OPT) vtk-generator.F90

cube-generator : cube-generator.F90
	$(F90) -o cube-generator $(OPT) cube-generator.F90

clean:
	rm -f vtk-generator cube-generator *.mod 
