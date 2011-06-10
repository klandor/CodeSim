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
#define GOPs 5
int STEPS = 6;
double Delta = 0.01, BASE = 1.05, errorRate=0.0;
int LAYERS = 3, PACKET_SIZE=1, Run=100;

int main(int argn, char **args) {
	int start_time = time(0);
	CRandomMersenne random(start_time);
	
	string precoder = "convo-2-3-6.txt", postcoder = "LT_ArithmWtAt0.2.txt"
		, interleaver = "interleaver-block-160x8.txt";
	
	ifstream ifs[MAX_LAYERSxGOPs], ifsStreamSize;
	ofstream ofs[MAX_LAYERSxGOPs];
	Codeword<Bit> a[MAX_LAYERSxGOPs];
	int streamSize[MAX_LAYERSxGOPs];
	
	if(argn >= 2){
		precoder = args[1];
	}
	
	if(argn >= 3){
		postcoder = args[2];
	}
	
	if(argn >= 4){
		interleaver = args[3];
	}
	
	ConvoCode cc( precoder );

	cout << "precoder: \"" << precoder << "\" postcoder: \""<< postcoder << '"';
	cout << " interleaver: \"" << interleaver << '"' << endl;
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
	unsigned long total_k = 0;
	
	
	// read in data
	short *p;
	
	ifsStreamSize.open("Stream_Size.txt");
	if(ifsStreamSize.fail()){
		cerr << "Error: file \"Stream_Size\" does not exist." << endl;
		exit(-1);
	}
	for (int i=0; i<LAYERS*GOPs; i++) {
		ifsStreamSize >> streamSize[i];
		total_k += streamSize[i];
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
	
	ofstream realRate("realRate.txt");
	
	
	// concacenate to LAYERS streams
	vector <Codeword<Bit> > b(LAYERS,Codeword<Bit>());
	
//	for (int i=0; i<LAYERS*GOPs; i++) {
//		b[2-(i%LAYERS)].insert(b[2-(i%LAYERS)].end(), a[i].begin(), a[i].end());
//	}
	
	
//	// mux into one
//	int max_length=0;
//	for (int i=0; i<LAYERS; i++) {
//		if(b[i].size() > max_length){
//			max_length = b[i].size();
//		}
//	}
	
	Codeword<Bit> c[GOPs];
	
	for (int i=0; i<GOPs; i++) {
		int max_length=0;
		for (int j=0; j<LAYERS; j++) {
			if(streamSize[i*LAYERS+j] > max_length)
			{
				max_length = streamSize[i*LAYERS+j];
			}
		}
		
		c[i].assign(LAYERS*max_length, 0);
		for (int j=0; j<c[i].size(); j++) {
			if(j/LAYERS < streamSize[i*LAYERS + (j%LAYERS)])
			{
				c[i][j] = a[i*LAYERS + (j%LAYERS)][j/LAYERS];
			}
		}
	}
	
	
	// precoder encoding
	Codeword<Bit> d[GOPs];
	for (int i=0; i<GOPs; i++) {
		d[i] = cc.encode(c[i]);
	}
	
	
	// interleaver
	Permutator<Bit> inter1(interleaver, true);
	Codeword<Bit> id[GOPs];
	for (int i=0; i<GOPs; i++) {
		id[i]  = inter1.permutate(d[i]);
	}
	
	// convert to Byte
	Codeword<Byte> e[GOPs];
	for (int i=0; i<GOPs; i++) {
		e[i]  = BitToByteCoverter::convert(id[i]);
	}
	
	for (int i=0; i<GOPs; i++) {
		cout << e[i].size() << endl;
	}
	
	// simulating start
	#pragma omp parallel for schedule(dynamic) num_threads(6)
	for (int run=0; run<Run; run++) {
		vector<bool> error[MAX_STEPS][MAX_LAYERSxGOPs];
		unsigned long total_n[MAX_STEPS];
		
		int seed;
		#pragma omp critical
		seed = random.BRandom();
		CRandomMersenne rnd(seed); // local random
		
		for (int d=0; d<STEPS; d++) {
			for (int i=0; i<LAYERS*GOPs; i++) {
				error[d][i].assign(streamSize[i], 0);
			}
			total_n[d] = 0;
			
			for (int gop =0; gop<GOPs; gop++) {
				vector<bool> mask;
				mask.reserve((size_t)(e[gop].size()/(1-errorRate)*2));
				long received =0;
				
				while (received < e[gop].size()*(BASE+d*Delta) ) {
					if(rnd.Random() > errorRate) { // received
						mask.insert(mask.end(), PACKET_SIZE, 0);
						received += PACKET_SIZE;
					}
					else {  // erased
						mask.insert(mask.end(), PACKET_SIZE, 1);
					}
				}
				
				total_n[d] += mask.size();
				LT_sim<Byte> lt(e[gop].size(), mask.size(),  Dsize, Tags, Distribution, rnd.BRandom());
				Codeword<Byte> f = lt.encode(e[gop]);
				
				if(f.size() != mask.size()) {
					cerr << "Error: Codeword size doesn't match." << endl;
					exit(-1);
				}
				
				for (int i=0; i<f.size(); i++) {
					f[i].setErased(mask[i]);
				}
				
				Codeword<Byte> df = lt.decode(f);
				Codeword<Bit> ff = BitToByteCoverter::revert(df);
				ff=inter1.depermutate(ff);
				Codeword<Bit> g = cc.decode(ff);
				
				// compare result
				for (int b=0; b<g.size(); b++) {
					if( b/LAYERS < streamSize[gop*LAYERS+b%LAYERS] && !(g[b] == c[gop][b]) ) {
						error[d][gop*LAYERS+b%LAYERS][b/LAYERS] = 1;
					}
				}
			}
			
			
//			LTCode<Byte> lt(K, K*1.5, Dsize, Tags, Distribution, rnd.BRandom());
//			Codeword<Byte> f = lt.encode(e);
//			
//			// packetize and channel
//			for (int frame_base=0; frame_base<f.size(); frame_base+=PACKET_SIZE) {
//				if( rnd.Random() < d*Delta ) {
//					for (int i=0; i<PACKET_SIZE; i++) {
//						if(frame_base + i == f.size())
//							break;
//						
//						f[frame_base + i].setErased(true);
//					}
//				}
//				
//			}
//			Codeword<Byte> df = lt.decode(f);
//			Codeword<Bit> ff = BitToByteCoverter::revert(df);
//			ff=inter1.depermutate(ff);
//			Codeword<Bit> g = cc.decode(ff);
//			
//			// compare result
//			for (int s=0; s<LAYERS; s++) { // s: stream
//				int base=0;
//				for (int i=0; i<GOPs; i++) {
//					for (int b=0; b<streamSize[i*LAYERS+s]; b++) {
//						if(g[(base+b)*LAYERS+(2-s)] == c[(base+b)*LAYERS+(2-s)])
//							error[d][i*LAYERS+s][b]=false;
//						else {
//							error[d][i*LAYERS+s][b]=true;
//						}
//					}
//					base += streamSize[i*LAYERS+s];
//				}
//				
//			}
			
			
			
		}
		
		// output result
		#pragma omp critical
		{
			for (int i=0; i<LAYERS*GOPs; i++) {
				for (int d=STEPS-1; d >=0; d--) {
					for (int b=0; b<streamSize[i]; b++) {
						ofs[i] << error[d][i][b] << ' ';
					}
					ofs[i] << '\n';
					
				}
				ofs[i] << endl;
				realRate << endl;
			}
			for (int d=0; d<STEPS; d++) {
				realRate << total_n[d]*8/(double)total_k << '\t';
			}
		}
	}
	
	cout << "Time: " << time(0)-start_time << endl;
	system("PAUSE");
	
}
