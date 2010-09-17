CPPFLAGS = -O3

all: CodeSim randomc cmaes

cmaes: cmaes.o cmaes.h cmaes_interface.h
	
randomc: randomc.h mersenne.o userintf.o

CodeSim: CodeSim.h Bit.o ConvoCode.o LT.h LTCode.h 

	

.PHONY: all clean
clean: 
	@rm *.o