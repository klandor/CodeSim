/*
 *  SVC_UEP_exp.cpp
 *  new_LT
 *
 *  Created by 刁培倫 on 2011/4/23.
 *  Copyright 2011 __MyCompanyName__. All rights reserved.
 *
 */

#include "randomc.h"
#include "LTCode.h"
#include "ConvoCode.h"
#include <ctime>
#include <fstream>
#include <string>
#include <sstream>
#include <omp.h>

using namespace std;
using namespace CodeSim;

#define MAX_STEPS 25
#define MAX_LAYERSxGOPs 25
int STEPS = 25;
double Delta = 0.01;
int LAYERS = 3, GOPs = 5, PACKET_SIZE=500, Run=100;

int main(int argn, char **args) {
	int start_time = time(0);
	CRandomMersenne random(start_time);
	
	string precoder = "convo-4-6-10.txt", postcoder = "LT_10k_0.01x.txt";
	ConvoCode cc( precoder );
	
	ifstream ifs[MAX_LAYERSxGOPs], ifsStreamSize;
	ofstream ofs[MAX_LAYERSxGOPs];
	Codeword<Bit> a[MAX_LAYERSxGOPs];
	int streamSize[MAX_LAYERSxGOPs];
	
	
	// read LT parameters
	ifstream ifsLT(postcoder.c_str());
	if(ifsLT.fail()){
		cerr << "Error: file \""+postcoder+"\" does not exist." << endl;
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
	
	MaxN = K*1.5;
	
	
	
	// read in data
	short *p;
	
	ifsStreamSize.open("Stream_Size");
	if(ifsStreamSize.fail()){
		cerr << "Error: file \"Stream_Size\" does not exist." << endl;
		exit(-1);
	}
	for (int i=0; i<LAYERS*GOPs; i++) {
		ifsStreamSize >> streamSize[i];
		
		p = new short[streamSize[i]];
		
		ostringstream oss;
		oss << "Stream_GOP" << i/LAYERS << "_QL" << i%LAYERS;
		ifs[i].open( oss.str().c_str() );
		if(ifs[i].fail()){
			cerr << "Error: file \"" << oss.str() << "\" does not exist." << endl;
			exit(-1);
		}
		else {
			for (int j=0; j<streamSize[i]; j++) {
				ifs[i] >> p[j];
			}
			
			a[i].assign(p, p+streamSize[i]);
			if(a[i].size() != streamSize[i])
				cerr << "ERROR: stream " << i << ": size not matched" << endl;
			
		}
		ifs[i].close();
		delete [] p;
		
		// open ofs
		ofs[i].open( ("O"+oss.str()).c_str() );
		if(ofs[i].fail()){
			cerr << "Error: file \"" << "O"+oss.str() << "\" can not be opened." << endl;
			exit(-1);
		}
		
		
	}
	ifsStreamSize.close();
	
	
	// concacenate to LAYERS streams
	vector <Codeword<Bit> > b(LAYERS,Codeword<Bit>());
	
	for (int i=0; i<LAYERS*GOPs; i++) {
		b[2-(i%LAYERS)].insert(b[2-(i%LAYERS)].end(), a[i].begin(), a[i].end());
	}
	
	
	// mux into one
	int max_length=0;
	for (int i=0; i<LAYERS; i++) {
		if(b[i].size() > max_length){
			max_length = b[i].size();
		}
	}
	
	Codeword<Bit> c(LAYERS*max_length, 0);
	
	for (int i=0; i<c.size(); i++) {
		if(i/LAYERS < b[i%LAYERS].size())
			c[i] = b[i%LAYERS][i/LAYERS];
	}
	
	
	// precoder encoding
	Codeword<Bit> d = cc.encode(c);
	
	// interleaver
	Permutator<Bit> inter1("interleaver-block-20x8.txt", true);
	Codeword<Bit> id = inter1.permutate(d);
	
	
	// convert to Byte
	Codeword<Byte> e = BitToByteCoverter::convert(id);
	
	cout << e.size() << ' '<<MaxN << endl;
	// simulating start
	#pragma omp parallel for num_threads(4)
	for (int run=0; run<Run; run++) {
		vector<bool> error[MAX_STEPS][MAX_LAYERSxGOPs];
		//#pragma omp parallel for //num_threads(6)
		for (int d=0; d<STEPS; d++) {
			for (int i=0; i<LAYERS*GOPs; i++) {
				error[d][i].assign(streamSize[i], 0);
			}
			
			int seed;
			#pragma omp critical
			seed = random.BRandom();
			CRandomMersenne rnd(seed); // local random
			
			LTCode<Byte> lt(K, K*1.5, Dsize, Tags, Distribution, rnd.BRandom());
			Codeword<Byte> f = lt.encode(e);
			
			// packetize and channel
			for (int frame_base=0; frame_base<f.size(); frame_base+=PACKET_SIZE) {
				if( rnd.Random() < d*Delta ) {
					for (int i=0; i<PACKET_SIZE; i++) {
						if(frame_base + i == f.size())
							break;
						
						f[frame_base + i].setErased(true);
					}
				}
				
			}
			Codeword<Byte> df = lt.decode(f);
			Codeword<Bit> ff = BitToByteCoverter::revert(df);
			ff=inter1.depermutate(ff);
			Codeword<Bit> g = cc.decode(ff);
			
			// compare result
			for (int s=0; s<LAYERS; s++) { // s: stream
				int base=0;
				for (int i=0; i<GOPs; i++) {
					for (int b=0; b<streamSize[i*LAYERS+s]; b++) {
						if(g[(base+b)*LAYERS+(2-s)] == c[(base+b)*LAYERS+(2-s)])
							error[d][i*LAYERS+s][b]=false;
						else {
							error[d][i*LAYERS+s][b]=true;
						}
					}
					base += streamSize[i*LAYERS+s];
				}
				
			}
			
			
			
		}
		
		// output result
		#pragma omp critical
		{
			for (int i=0; i<LAYERS*GOPs; i++) {
				for (int d=0; d<STEPS; d++) {
					for (int b=0; b<streamSize[i]; b++) {
						ofs[i] << error[d][i][b] << ' ';
					}
					ofs[i] << '\n';
					
				}
				ofs[i] << '\n';
			}
		}
	}
	
	cout << "Time: " << time(0)-start_time << endl;
	system("PAUSE");
	
}
