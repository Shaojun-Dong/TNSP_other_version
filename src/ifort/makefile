include ../../make.inc

ifort := ifort.o



all: $(ifort) 

$(ifort):%.o : ./%.f90

./%.o :	./%.f90
	$(FC) $(openmp)  -o $@ -c $(<) -I../Tensor-$(version) -I../Symmetry-$(Symmetry_version) $(Optionl)

	
clean:
	rm -rf *.o *.mod 

