CPPFLAGS = -O3 -fopenmp

all: CodeSim randomc cmaes

CMAES.out: all CMAES_main.cpp
	$(CXX) $(CPPFLAGS) CMAES_main.cpp mersenne.o userintf.o cmaes.o Bit.o -o $@

histogram.out: all LT_sim_histogram.cpp
	$(CXX) $(CPPFLAGS) LT_sim_histogram.cpp mersenne.o userintf.o cmaes.o Bit.o -o $@

exp_SVC.out : randomc CodeSim exp_SVC.cpp
	$(CXX) $(CPPFLAGS) exp_SVC.cpp mersenne.o userintf.o ConvoCode.o Bit.o -o $@

decoding_pattern_test.out: CodeSim decoding_pattern_test.cpp
	$(CXX) $(CPPFLAGS) decoding_pattern_test.cpp mersenne.o userintf.o Bit.o -o $@


cmaes: cmaes.o cmaes.h cmaes_interface.h
	
randomc: randomc.h mersenne.o userintf.o

CodeSim: CodeSim.h Bit.o ConvoCode.o LT.h LTCode.h 

	

.PHONY: clean
clean: 
	@rm *.o