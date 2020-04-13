include ../../make.inc


fTensor :=\
QuantumNumber.o\
SymDimension.o\
SymTensor.o\
U1.o\
Parity.o\
fTensor.o



all: $(fTensor) 

$(fTensor):%.o : ./%.f90

./%.o :	./%.f90
	$(FC) $(openmp) -cpp -o $@ -c $(<)  $(Optionl) -I../Tensor-$(version)  
	
clean:
	rm -rf *.o *.mod *.a

