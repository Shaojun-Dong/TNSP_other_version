include ../../make.inc




fTensor :=\
SymDimension.o\
SymTensor.o\
fTensor.o

PTensor :=\
SymDimension.o\
SymTensor.o\
PTensor.o



all: $(fTensor) 

$(fTensor):%.o : ./%.f90

./%.o :	./%.f90
	$(FC) $(openmp)  -o $@ -c $(<)  $(Optionl) -I../Tensor-$(version)  
	
clean:
	rm -rf *.o *.mod *.a

