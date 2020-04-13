include ../../make.inc




fTensor :=\
QuantumNumber.o\
SymDimension.o\
U1.o\
Parity.o\
SymTensor.o\
fTensor.o




all: $(fTensor) 

$(fTensor):%.o : ./%.f90

./%.o :	./%.f90
	$(FC) $(openmp) -cpp  -DIFORT=$(IFIFORT) -o $@ -c $(<)  $(Optionl) -I../Tensor-$(version)  
	
clean:
	rm -rf *.o *.mod *.a

