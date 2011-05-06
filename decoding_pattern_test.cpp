/*
 *  consecutive_test.cpp
 *  new_LT
 *
 *  Created by 刁培倫 on 2010/10/24.
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
 */


#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <ctime>    
#include <cstdlib>
#include <vector>
#include <math.h>
#include "randomc.h"
#include "LTCode.h"
#include <omp.h>
using namespace std;
using namespace CodeSim;

//#define K 1000
int K;

#define EPS 0.2
#define RUN 100
#define BASE 1.08
#define Delta 0.005
#define STEPS 6
//#define MAX_WINDOW_SIZE 100

int tags = 10, windowSize=50;

int		*Degree;//[10] = {1, 2, 3, 4, 5, 7, 9, 19, 59, 179};//{1,2,3,4,5,8,9,19,65,66};

double  *Omega;//[10] ;//= {7.9379E-02, 4.0129E-01, 1.0121E-01, 2.1679E-01, 5.0996E-02, 
//	5.8338E-05, 3.9740E-02, 7.7470E-02, 2.1520E-02, 1.1547E-02};
int *err_histo[STEPS]; //, l0a_max[STEPS][1000];

int main(){
	cin >> K >> tags;
	Degree = new int[tags];
	Omega = new double[tags];
	for (int i=0; i<tags; i++) {
		cin >> Degree[i];
	}
	for (int i=0; i<tags; i++) {
		cin >> Omega[i];
	}
	CRandomMersenne r(time(0));
//	Permutator<Bit> per("Interleaver.txt", true);
	
	for (int s=0; s<STEPS; s++) {
		
		err_histo[s] = new int[windowSize+1];
		for (int i=0; i<windowSize+1; i++) {
			err_histo[s][i]=0;
		}
	}
	
	#pragma omp parallel for num_threads(6)
	for (int run=0; run<RUN; run++) {
		LT_sim<Bit> lt(K, K*(BASE+STEPS*Delta), tags, Degree, Omega, r.BRandom());
		
		for (int s=0; s<STEPS; s++) {
//			int l0h[1000];
//			err_histo[s] = new int[windowSize+1];
//			for (int i=0; i<windowSize+1; i++) {
//				err_histo[s][i]=0;
//			}
			
			lt.seqReceive(K*(BASE+s*Delta)-1);
			Codeword<Bit> a = lt.getResult(); // decoding result
			
			// start compute histogram
			int errNO=0;
			
			// first 'windowSize' bits
			for (int p=0; p<windowSize; p++) {
				if (a[p].isErased()) {
					errNO ++;
				}
			}
			
			#pragma omp atomic
			err_histo[s][errNO]++;
			
			
			for (int p=windowSize; p< a.size(); p++) {
				if (a[p].isErased()) {
					errNO++;
				}
				if (a[p-windowSize].isErased()) {
					errNO --;
				}
				
				
				#pragma omp atomic
				err_histo[s][errNO]++;
			}
			
			// end compute histogram
	
		} // for STEPS
		
	}// for RUN
	
	// print out histogram
	for (int s=0; s<STEPS; s++) {
		cout << BASE + s*Delta << ":\t";
		for (int i=0; i<windowSize+1; i++) {
			cout << err_histo[s][i] / (double)(RUN*(K-windowSize+1)) << '\t' ;//<< l0a_max[i] << '\n';
		}
		cout << '\n';
	}

	
}