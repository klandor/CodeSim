/*
 *  exp_SVC.cpp
 *  new_LT
 *
 *  Created by 刁培倫 on 2010/10/15.
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
 */



#include <iostream>
#include "LTCode.h"
#include "ConvoCode.h"
#include "randomc.h"
#include <ctime>
#include <omp.h>

#define K 60000
#define LT_K 1000
#define RUN 1
using namespace std;
using namespace CodeSim;

int main(){
	int Degree_of_Edge[10] = {1, 2, 3, 4, 5, 7, 9, 19, 59, 179};
	double Omega[10] = {7.9379E-02, 4.0129E-01, 1.0121E-01, 2.1679E-01, 5.0996E-02, 
		5.8338E-05, 3.9740E-02, 7.7470E-02, 2.1520E-02, 1.1547E-02};
	
	
	CRandomMersenne random(time(0));
	ConvoCode convo("convo-4-6-10.txt");
	Permutator<Bit> permut("Interleaver.txt", true);
	
	long ber[21][3];
	for (int i=0; i<21*3; i++) {
		ber[0][i] = 0;
	}
	
	#pragma omp parallel for num_threads(6)
	for (int run = 0; run < RUN; run++) {
		
		Codeword<Bit> m(K,0);
		for (Codeword<Bit>::iterator i = m.begin(); i!=m.end(); i++) {
			*i = random.IRandomX(0, 1);
		}
		LTCode<Bit> lt(LT_K, LT_K*1.2, 10, Degree_of_Edge, Omega, random.BRandom());
		
		
		Codeword<Bit> mid =  convo.encode(m) ;
		mid = permut.permutate(mid);
		mid = lt.encode(mid);
		
		for (int e=0; e<21; e++) 
		{
			
			for (int i=0; i<mid.size(); i++) {
				if ( (i%1200) < LT_K*(1+e*0.01)) {
					mid[i].setErased(false);
				}
				else {
					mid[i].setErased(true);
				}

			}
			
			Codeword<Bit> out = lt.decode(mid);
			out = permut.depermutate(out);
			out = convo.decode(out);
			
			for (int i=0; i<out.size(); i++) {
				if ( !(out[i]==m[i]) ) {
					#pragma omp atomic
					ber[e][i%3]++;
				}
			}
		}
		
	}
	
	cout << "BER:\n";
	for (int i=0; i<3; i++) {
		for (int j=0; j<21; j++) {
			cout << ber[j][i] / (K/3.0*RUN) << ' ';
		}
		cout << '\n';
	}
	//<< ber[0] / (double) (K*RUN/3) << ' '<< ber[1] / (double) (K*RUN/3) << ' '<< ber[2] / (double) (K*RUN/3) << '\n';
	
	return 0;

}