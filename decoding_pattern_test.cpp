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

#define K 1000
#define LEN 100
#define EPS 0.2
#define RUN 1000
#define WINDOW_SIZE 10

int tags = 10;

int		Degree[10] = //{1,2,3,4,5,8,9,19,65,66};
{1, 2, 3, 4, 5, 7, 9, 19, 59, 179};
double  Omega[10] = {7.9379E-02, 4.0129E-01, 1.0121E-01, 2.1679E-01, 5.0996E-02, 
	5.8338E-05, 3.9740E-02, 7.7470E-02, 2.1520E-02, 1.1547E-02};

int main(){
	CRandomMersenne r(time(0));
	int l0a[1000], l0a_max[1000], ber=0;
	for (int i=0; i<1000; i++) {
		l0a[i]=0;
		l0a_max[i]=0;
	}
	
	#pragma omp parallel for num_threads(6)
	for (int run=0; run<RUN; run++) {
		LTCode<Bit> lt(K, K*(1+EPS), tags, Degree, Omega, r.BRandom());
		int l0 =0, l0h[1000];
		for (int i=0; i<1000; i++) {
			l0h[i]=0;
			
		}
		Codeword<Bit> a(K*LEN,0);
		a = lt.encode(a);
		a = lt.decode(a);
		a.insert(a.end(), a.begin(), a.begin() + (WINDOW_SIZE-1));
		
		int errNO=0, errLen=0;
		for (int p=0; p<WINDOW_SIZE; p++) {
			if (a[p].isErased()) {
				errNO ++;
			}
		}
//		if(errNO/(double)winSize > errorDensityBound)
//			errLen=1;
		//cout << errNO << '\n';
		#pragma omp atomic
		l0a[errNO]++;
		#pragma omp atomic
		l0h[errNO]++;
		for (int p=WINDOW_SIZE; p< a.size(); p++) {
			if (a[p].isErased()) {
				errNO++;
			}
			if (a[p-WINDOW_SIZE].isErased()) {
				errNO --;
			}
			
			//if(errNO/(double)winSize > errorDensityBound)
//			{
//				errLen ++;
//			}
//			else {
//				if (errLen > 750) {
//					//fit +=failurePenalty[i];
//#pragma omp atomic
//					failureCount[i] += errLen / 750.0;
//				}
//				errLen = 0;
//			}
			#pragma omp atomic
			l0a[errNO]++;
			#pragma omp atomic
			l0h[errNO]++;
			//cout << errNO << '\n';
		}
// counting consecutive 0 or 1		
//		int l0 =0, l0h[1000];
//		for (int i=0; i<1000; i++) {
//			l0h[i]=0;
//			
//		}
//		for (int i=0; i<a.size(); i++) {
//			if (! a[i].isErased()) {
//				l0++;
//			}
//			else {
//				if (l0>10000) {
//					cout << "l0: " <<l0<< endl;
//				}
//				else {
//					ber++;
//					l0a[l0/10]+=1;
//					l0h[l0/10]+=1;
//					l0 = 0;
//				}
//			}
//
//		}
//		if (l0>0) {
//			l0a[l0/10]+=1;
//			l0h[l0/10]+=1;
//		}
//		
		for (int i =0; i<WINDOW_SIZE+1; i++) {
			#pragma omp critical
			if (l0h[i] > l0a_max[i]) {
				l0a_max[i] = l0h[i];
			}
		}
		
		
		
	}
	cout.precision(10);
	//cout << "BER: "<<ber/(double)(K*RUN)<< endl;
	
	for (int i=0; i<WINDOW_SIZE+1; i++) {
		cout << l0a[i] / (double)RUN << ' ' << l0a_max[i] << '\n';
	}
	
}