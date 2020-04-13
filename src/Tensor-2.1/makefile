include ../../make.inc



Tensor := \
function.o     \
Dimension.o      \
Tensor_real.o\
Tensor.o   \
eigen.o

all: $(Tensor) 
	
$(Tensor):%.o : ./%.f90

./%.o :	./%.f90
	$(FC) $(openmp)  -o $@ -c $(<)  $(Optionl)
	
clean:
	rm -rf *.o *.mod *.a

