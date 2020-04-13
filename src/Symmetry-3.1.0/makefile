include ../../make.inc




fTensor :=\
QuantumNumber.o\
SymDimension.o\
SymTensor.o\
fTensor.o

PTensor :=\
QuantumNumber.o\
SymDimension.o\
SymTensor.o\
PTensor.o

U1fTensor :=\
QuantumNumber.o\
SymDimension.o\
SymTensor.o\
U1.o\
U1fTensor.o

PfTensor :=\
QuantumNumber.o\
SymDimension.o\
SymTensor.o\
Parity.o\
PfTensor.o


all: $(fTensor) $(U1fTensor) $(PfTensor)

$(fTensor):%.o : ./%.f90

./%.o :	./%.f90
	$(FC) $(openmp)  -o $@ -c $(<)  $(Optionl) -I../Tensor-$(version)  
	
clean:
	rm -rf *.o *.mod *.a

