CPPFLAGS = -O3 -fopenmp

all: CodeSim randomc cmaes

CMAES_main: all CMAES_main.cpp
	$(CXX) $(CPPFLAGS) CMAES_main.cpp mersenne.o userintf.o cmaes.o Bit.o -o CMAES_main.out


cmaes: cmaes.o cmaes.h cmaes_interface.h
	
randomc: randomc.h mersenne.o userintf.o

CodeSim: CodeSim.h Bit.o ConvoCode.o LT.h LTCode.h 

	

.PHONY: clean
clean: 
	@rm *.o