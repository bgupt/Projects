MODULES = kinds.o vars.o rk4.o rhs.o
F90 = ifort
#F90 = gfortran

F90FLAGS = -O2 -g -openmp -openmp-report2

all:	evol.x

evol.x:	main.o $(MODULES)
	$(F90) $(F90FLAGS) main.o $(MODULES) -o evol.x

main.o:	main.f90 $(MODULES)
	$(F90) $(F90FLAGS) -c main.f90

rhs.o: rhs.f90 kinds.o vars.o
	$(F90) $(F90FLAGS) -c rhs.f90

rk4.o: rk4.f90 rhs.o kinds.o vars.o
	$(F90) $(F90FLAGS) -c rk4.f90

kinds.o:	kinds.f90
	$(F90) $(F90FLAGS) -c kinds.f90

vars.o:	vars.f90 kinds.o
	$(F90) $(F90FLAGS) -c vars.f90

clean:
	rm -rf *.o *.mod *.x
