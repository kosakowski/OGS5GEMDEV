# By Wenqing Wang
OPTFLAG = -O3
DEBFLAG = #-ggdb
C++    =  g++ $(DEBFLAG)  -Wall 

HFile = matrix_class.h mesh.h 
CFile = partition.cpp matrix_class.cpp mesh.cpp 
OFile = partition.o  matrix_class.o mesh.o

.SUFFIXES: .o .cpp .h 

.cpp.o:
	$(C++) -c -o $*.o $<

main: 	$(OFile) 
	$(C++) -o m2g $(OFile)

clean:  
	rm -f *.o
	rm m2g


m2g: partition.h  matrix_class.h mesh.h 

