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
#define RUN 1000000
#define BASE 1.20
#define STEP_SIZE 0.01
#define STEPS 6
//#define MAX_WINDOW_SIZE 100

int tags = 10, windowSize=50;

int		*Degree;//[10] = {1, 2, 3, 4, 5, 7, 9, 19, 59, 179};//{1,2,3,4,5,8,9,19,65,66};

double  *Omega;//[10] ;//= {7.9379E-02, 4.0129E-01, 1.0121E-01, 2.1679E-01, 5.0996E-02, 
//	5.8338E-05, 3.9740E-02, 7.7470E-02, 2.1520E-02, 1.1547E-02};
int *l0a[RUN][STEPS]; //, l0a_max[STEPS][1000];

int main(){
	string s;
//	cerr << "Enter output filename: ";
//	getline(cin, s, '\n');
//	ofstream fout(s.c_str());
//	
//	cerr << "Enter comment: ";
//	getline(cin, s, '\n');
//	fout << "Comment: " << s << endl;
	
//	cerr << "Enter Window Size(>0): ";
//	cin >> windowSize;
	
//	cerr << "Enter degree distribution: ";
	cin >> tags;
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
	
//	for (int i=0; i<1000*STEPS; i++) {
//		l0a[0][i]=0;
//		l0a_max[0][i]=0;
//	}
	
	#pragma omp parallel for num_threads(6)
	for (int run=0; run<RUN; run++) {
		LT_sim<Bit> lt(K, K*(BASE+STEPS*STEP_SIZE), tags, Degree, Omega, r.BRandom());
		
		for (int s=0; s<STEPS; s++) {
//			int l0h[1000];
			l0a[run][s] = new int[windowSize+1];
			for (int i=0; i<windowSize+1; i++) {
				l0a[run][s][i]=0;
			}
			
			lt.seqReceive(K*(BASE+s*STEP_SIZE)-1);
			//a.insert(a.end(), a.begin(), a.begin() + (windowSize-1));
			Codeword<Bit> a = lt.getResult();
			int errNO=0;
			for (int p=0; p<windowSize; p++) {
				if (a[p].isErased()) {
					errNO ++;
				}
			}
			
			//#pragma omp atomic
			l0a[run][s][errNO]++;
//			l0h[errNO]++;
			for (int p=windowSize; p< a.size(); p++) {
				if (a[p].isErased()) {
					errNO++;
				}
				if (a[p-windowSize].isErased()) {
					errNO --;
				}
				
				
				//#pragma omp atomic
				l0a[run][s][errNO]++;
//				l0h[errNO]++;
				//cout << errNO << '\n';
			}
	
//			for (int i =0; i<windowSize+1; i++) {
//				#pragma omp critical
//				if (l0h[i] > l0a_max[s][i]) {
//					l0a_max[s][i] = l0h[i];
//				}
//			}
		}
		


		
		
		
	}
	//cout.precision(10);
	//cout << "BER: "<<ber/(double)(K*RUN)<< endl;
//	fout << "Average:\n";
	for (int s=0; s<STEPS; s++) {
		for (int run = 0; run<RUN; run++) {
			
			for (int i=0; i<windowSize+1; i++) {
				cout << l0a[run][s][i] / (double)((K-windowSize+1)) << '\t' ;//<< l0a_max[i] << '\n';
			}
			cout << '\n';
		}
		
		cout << "\n\n\n\n\n\n\n\n\n\n";
	}
//	fout << "\nMaximum:\n";
//	for (int s=0; s<STEPS; s++) {
//		for (int i=0; i<windowSize+1; i++) {
//			fout << l0a_max[s][i] / (double)(K-windowSize+1)<< ' ' ;//<< l0a_max[i] << '\n';
//		}
//		fout << '\n';
//	}
	
}