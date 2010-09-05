/*
 *  LTTest.cpp
 *  new_LT
 *
 *  Created by 刁培倫 on 2010/9/2.
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
 */

#include "LTCode.h"
using namespace CodeSim;

int main(){
	int Degree_of_Edge[10] = {1, 2, 3, 4, 5, 7, 9, 19, 59, 179};
	double Omega[10] = {7.9379E-02, 4.0129E-01, 1.0121E-01, 2.1679E-01, 5.0996E-02, 
		5.8338E-05, 3.9740E-02, 7.7470E-02, 2.1520E-02, 1.1547E-02};
	double sum[16];
	for (int i = 0; i< 16; i++)
		sum[i] = 0;
	
	CRandomMersenne r(12345);
	
	int K = 1000;
	for (int run=0; run<500; run++) {
		//cout << "Run: " << run << "\n";
		LTCode<Bit> lt(K, K*1.2, 10, Degree_of_Edge, Omega, r.BRandom());
		
		Codeword<Bit> in(5*K, 0);
		for (int i=1; i<5*K; i+=2) {
			in[i] = 1;
		}
		Codeword<Bit> c = lt.encode(in);
		
		for (int q=0; q<5; q++){
			for (int i=K*1.04; i<K*1.2; i++) {
				c[q*K*1.2 + i].setErased(true);
			}
		}
		
		for (int i = 0; i< 16; i++) {
			for (int q=0; q<5; q++){
				for (int j=K*(1.04+0.01*i); j<K*(1.05+0.01*i); j++) {
					c[q*K*1.2 +j].setErased(false);
				}
			}
			Codeword<Bit> out = lt.decode(c);
			
			
			//sum[i] += sim.run();
			//cout << sim.run();
			int s=0;
			for (int j = 0; j<5*K; j++) {
				if ( !out[j].isErased() && in[j] == out[j]) 
				{
					s++;
				}
				
			}
			sum[i] += 1- s/(double)( 5*K) ;
		}
	}
	for (int i = 0; i< 16; i++)
		cout << sum[i] / 500 << '\t';
	
	
    return 0;
	
}