/*
 *  LT.h
 *  new_LT
 *
 *  Created by 刁培倫 on 2010/8/6.
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
 */

#include <vector>
#include "randomc.h"
#include "Code.h"
using namespace std;

namespace LT_class{


	class LT_sim : public CodeSim::Code{
	public:
		LT_sim();
		LT_sim(int k, int max_n, int tag_size, int *tags, double * omega, int seed);
		void init(int k, int max_n, int tag_size, int *tags, double * omega, int seed);
		void seqReceive(int t);
		void receive(int t);
		void decode();
		bool isDecoded(int t);
		double run();
		vector<char> encode(vector<char> a);
		vector<char> decode(vector<char> a);
		void reset();
	private:
		void codeGen(int t);
		
		long int Num_of_Input, Num_of_Output, Num_of_Degree, Num_of_Decoding, 
				ReceivedSize, generatedCode;
		CRandomMersenne Rnd;
		vector<long int>	d;//(Num_of_Output, 0);
		vector<vector<long int> >	edge;// = new vector<long int>[Num_of_Output];
		vector<long int>	R_M;//(Num_of_Output, 0);
		//vector<long int>	erasure;//(Num_of_Output, 0);
		vector<long int>	DE;//(Num_of_Input, 0);
		vector<int>			Degree_of_Edge;
		vector<double>		Omega;
		vector<int>			receivedMask;
	};
}