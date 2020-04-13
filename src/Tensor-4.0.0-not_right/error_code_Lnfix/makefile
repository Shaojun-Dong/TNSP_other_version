include $(HOME)/makefileConf

OBJ:= \
lnfix.o 

run:	$(OBJ)
	$(FC)  -fopenmp -o $@ $^  $(TensorDIRLink) $(TensorlibLink) $(PACKAGE)
./%.o :	./%.f90
	$(FC) -fopenmp  -o $@ -c $(<) $(TensorDIRLink)
clean:
	rm *.o *.mod 
cleanAll:
	rm *.o *.mod run
cleanEx:
	rm run
