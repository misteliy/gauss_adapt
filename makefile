test_gauss : test_gauss.o random.o determinante.o 
	gfortran -o test_gauss test_gauss.o random.o determinante.o  -lblas -llapack
	
test_gauss.o:test_gauss.f90 random.o determinante.o 
	gfortran -c test_gauss.f90
random.o:random.f90 
	gfortran -c random.f90
determinante.o:determinante.f90
	gfortran -c determinante.f90
