CPPFLAGS = -O3 -fopenmp

all: CodeSim randomc cmaes

CMAES.out: all CMAES_main.cpp
	$(CXX) $(CPPFLAGS) CMAES_main.cpp mersenne.o userintf.o cmaes.o Bit.o -o $@

CMAES2.out: all CMAES2.cpp
	$(CXX) $(CPPFLAGS) CMAES2.cpp mersenne.o userintf.o cmaes.o Bit.o -o $@

CMAES3.out: all CMAES3.cpp
	$(CXX) $(CPPFLAGS) CMAES3.cpp mersenne.o userintf.o cmaes.o Bit.o -o $@
	
histogram.out: all LT_sim_histogram.cpp
	$(CXX) $(CPPFLAGS) LT_sim_histogram.cpp mersenne.o userintf.o cmaes.o Bit.o -o $@

exp_SVC.out : randomc CodeSim exp_SVC.cpp
	$(CXX) $(CPPFLAGS) exp_SVC.cpp mersenne.o userintf.o ConvoCode.o Bit.o -o $@

exp_LT_sim.out : randomc CodeSim exp_LT_sim.cpp
	$(CXX) $(CPPFLAGS) exp_LT_sim.cpp mersenne.o userintf.o ConvoCode.o Bit.o -o $@

decoding_pattern_test.out: CodeSim decoding_pattern_test.cpp
	$(CXX) $(CPPFLAGS) decoding_pattern_test.cpp mersenne.o userintf.o Bit.o -o $@

fitting_data_set_generator.out: fitting_data_set_generator.cpp
	$(CXX) $(CPPFLAGS) fitting_data_set_generator.cpp -o $@

cmaes: cmaes.o cmaes.h cmaes_interface.h
	
randomc: randomc.h mersenne.o userintf.o

CodeSim: CodeSim.h Bit.o ConvoCode.o LT.h LTCode.h 

commitall:
	git commit -a

push:
	git push github master
pull:
	git pull github master

.PHONY: clean commitall push pull
clean: 
	@rm *.o