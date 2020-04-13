include ../../make.inc

mpi := \
mpi.o



all: $(mpi) 

$(mpi):%.o : ./%.f90

./%.o :	./%.f90
	$(FC) $(openmp)  -o $@ -c $(<)  $(Optionl)
	
clean:
	rm -rf *.o *.mod 

