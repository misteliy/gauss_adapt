test_gauss : test_gauss.o random.o
	gfortran -o test_gauss test_gauss.o random.o -lblas -llapack
	
test_gauss.o:test_gauss.f90 random.o
	gfortran -c test_gauss.f90
random.o:random.f90 
	gfortran -c random.f90
