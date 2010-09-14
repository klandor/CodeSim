CC = g++ -O3

program: CodeSim random

	
random: randomc.h mersenne.o userintf.o

CodeSim: CodeSim.h Bit.o ConvoCode.o LT.h LTCode.h 

Bit.o: Bit.cpp CodeSim.h
	$(CC) -c -o Bit.o Bit.cpp
	
ConvoCode.o: ConvoCode.cpp ConvoCode.h CodeSim.h
	$(CC) -c -o ConvoCode.o ConvoCode.cpp

mersenne.o: mersenne.cpp
	$(CC) -c -o mersenne.o mersenne.cpp
	
userintf.o: userintf.cpp
	$(CC) -c -o userintf.o userintf.cpp

clean: 
	@rm *.o