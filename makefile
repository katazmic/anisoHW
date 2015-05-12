
#GSL_L	=	/usr/local/lib/
#GSL_I	=	/usr/local/include/


GSL_L	=	/usr/lib/
GSL_I	=	/usr/include/



clean: rm *.o

shell_aniso_hw: work.o HW_aniso.o
	g++  -I$(GSL_I) -L$(GSL_L) -fopenmp -lgsl -lgslcblas -lm  work.o HW_aniso.o -o shell_aniso_hw

work.o: work.cpp HW_aniso.h
	g++ -c work.cpp

HW_aniso.o: HW_aniso.cpp HW_aniso.h
	g++ -c HW_aniso.cpp 
