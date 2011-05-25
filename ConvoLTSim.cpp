/*
 *  ConvoLTSim.cpp
 *  new_LT
 *
 *  Created by 刁培倫 on 2011/5/25.
 *  Copyright 2011 __MyCompanyName__. All rights reserved.
 *
 */


#include "ConvoCode.h"
#include "LTCode.h"
#include "randomc.h"
#include <ctime>
#include <cmath>
//#define L 100000
#define MAX_BIT 1000000000L
#define BASE 1.05
#define STEPS 6
#define Delta 0.01
#include <omp.h>

using namespace CodeSim;
using namespace std;

int main(int argn, char **args){
	
	if (argn < 3) {
		cerr << "Usage: ConvoLTSim.out convo_code_file LT_file [punchering_table_string]" << endl;
		exit(-1);
	}
	
	int start_time = time(0);
	string puncher_table = "1";
	
	
	if (argn > 3) {
		puncher_table = args[3];
	}
	
	CRandomMersenne r(time(0));
	ConvoCode cc( args[1] );
	int Layer = cc.getK();
	unsigned long L = 80000/Layer;
	unsigned long Run = MAX_BIT / 80000;
	ifstream ifsLT(args[2]);
	if(ifsLT.fail()){
		cerr << "Error: file \""+string(args[2])+"\" does not exist." << endl;
		exit(-1);
	}
	int K, MaxN, Dsize, *Tags;
	double *Distribution;
	
	ifsLT >> K >> Dsize;
	
	Tags = new int[Dsize];
	Distribution = new double[Dsize];
	
	for (int i=0; i<Dsize; i++) {
		ifsLT >> Tags[i];
	}
	for (int i=0; i<Dsize; i++) {
		ifsLT >> Distribution[i];
	}
	
	
	Permutator<Bit> inter("interleaver-block-20x8.txt");
	
	
	cout << "ConvoCode file: \"" << args[1] <<"\" LT file: " << args[2] << endl;
	cout << "Epsilon";
	for (int i=0; i<Layer; i++) {
		cout << "\tLayer" << i+1;
	}
	cout << "\tTotalBits\n";
	
	
	for (int s=0; s<STEPS; s++) {
		vector<unsigned long> err(Layer, 0);
		unsigned long total=0;
		double Epsilon = (BASE+s*Delta) -1;
		//		for (int i=0; i<Layer; i++) {
		//			err[i]=0;
		//		}
		
		//while (1) 
		{
			#pragma omp parallel for num_threads(6)
			for (int i=0; i<Run; i++) {
				Codeword<Bit> a;
				a.reserve(Layer*L);
				for (int t=0; t<Layer*L; t++) {
					a.push_back(r.IRandomX(0, 1));
				}
				Codeword<Bit> b = cc.encode(a);
				b = inter.permutate(b);
				Codeword<Byte> b1 = BitToByteCoverter::convert(b);
				
				if(s == 0 && i==0)
					cout << "LT block:" << b1.size() << endl;
				LT_sim<Byte> lt(b1.size(), b1.size()*(1+Epsilon), Dsize, Tags, Distribution, r.BRandom());
				Codeword<Byte> c1 = lt.encode(b1);
				c1 = lt.decode(c1);
				Codeword<Bit> c2 = BitToByteCoverter::revert(c1);
				c2 = inter.depermutate(c2);
				
				Codeword<Bit> c = cc.decode(c2);
				
				for (int i=0; i<c.size(); i++) {
					if (!(a[i] == c[i])) {
						#pragma omp atomic
						err[i%Layer]++;
					}
				}
				#pragma omp atomic
				total += L;
			}
			
		
		}
		
		// enough # of errors
		
		cout << Epsilon;
		
		for (int i=0; i<Layer; i++) {
			cout << '\t' << err[i]/(double)total;
		}
		
		cout << '\t' << total << endl;
		
		//errRate *= pow(10, -1.0/STEPS_PER_ORDER);
		Epsilon += 0.01;
		
		
	}
	
	cout << "Time: " << time(0)-start_time << endl;
	
}