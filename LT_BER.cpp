#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <ctime>    
#include <cstdlib>
#include <vector>
#include <math.h>
#include "randomc.h"
#include "LT.h"
#include "statistics.h"
#include <omp.h>
#include <algorithm>

int K;

int Run;
//#define C 0.05
#define Delta 0.03
int STEPS;
#define MaxN (K*(1+Delta*(STEPS-1)))

using namespace std;
using namespace CodeSim;


//unsigned long ErrorCount[STEPS][16];
//double BER[STEPS];
vector< unsigned long > BER;
int Dsize;


int* Degree;
//double* Omega;
//int		Degree[10] = //{1,2,3,4,5,8,9,19,65,66};
//		{1, 2, 3, 4, 5, 7, 9, 19, 59, 179};
//double  Omega[10] = {7.9379E-02, 4.0129E-01, 1.0121E-01, 2.1679E-01, 5.0996E-02,
//				5.8338E-05, 3.9740E-02, 7.7470E-02, 2.1520E-02, 1.1547E-02
//};

//uniformRandom *Rnd;
int g_seed = (int)time(0);
CRandomMersenne Rnd(g_seed);


double* D;
double* SD;


double e = 1.05;					



int main(){
	//int i;
	//Rnd = new ran3();
	
	//D = Robust_Soliton_Distribution(K,C,Delta);
	
	//cout << "haha";
	cin >> K;
	cin >> Run;
	cin >> Dsize;
	Degree = new int[Dsize];
	D = new double[Dsize];
	SD = new double[Dsize];
	for (int i=0; i<Dsize; i++) {
		cin >> Degree[i];
	}
//	for (int i=0; i<Dsize; i++) {
//		cin >> D[i];
//	}	
	cin >> STEPS;
	
	
	vector<double> epsilons(STEPS, 0);
	for (int i=0; i<STEPS; i++) {
		cin >> epsilons[i];
	}	
	
//	vector<double> sum(STEPS,0), mean(STEPS,0), var(STEPS,0),
//	dev(STEPS,0), skew(STEPS,0), kurt(STEPS,0), FER(STEPS,0);
	while (cin >> D[0])
	{
//		cout << "distribution\t";
		for (int i =1; i<Dsize; i++) {
			cin >> D[i];
//			cout << D[i] << '\t';
			
		}
		
		BER.assign(STEPS, 0);
//		cout << "\nEpsilons\t";
//		for (int i = 0; i<STEPS; i++) {
//			cout << i*Delta << '\t';
//			//BER[i] = 0;
//			BER[i].assign(Run,0);
//			for(int j = 0; j< 16; j++)
//				ErrorCount[i][j] = 0;
//		}
//		cout << '\n';
		
//		int start_time = time(0);
		//int n=0;
		#pragma omp parallel for schedule(dynamic) num_threads(6)
		for(int run=0;run<Run;run++){
			int seed;
			#pragma omp critical
			seed = Rnd.BRandom();
			
			LT_sim<Bit> sim(K, K*(1+epsilons[STEPS-1]), Dsize, Degree, D, seed);
			//cout << "Run " << i <<'\n';
			//#pragma omp parallel for
			for (int i = 0; i< STEPS; i++) {
				sim.seqReceive( K*(1+epsilons[i]) -1);
				//sim.decode();
				double t = sim.failureRate();// = Encoder(K, K*(1.05+0.01*i), Dsize);
//				if (t<1) {
//					#pragma omp atomic
//					ErrorCount[i][(int)(t*16)]++;
//				}
//				else {
//					#pragma omp atomic
//					ErrorCount[i][15]++;
//				}
				
				BER[i] += t*K;
//				if (t>0) {
//					FER[i]++;
//				}
			}
			
			//cout << '\n';
		}
		//cout<<"Histogram";
//		for (int i = 0; i<16; i++) {
//			cout << (i+1)/16.0<<'\t';
//			for(int j = 0; j< STEPS; j++)
//				cout <<  ErrorCount[j][i] / (double)Run<< '\t';
//			cout << '\n';
//		}
		
//		#pragma omp parallel for num_threads(6)
//		for (int i = 0; i<STEPS; i++) {
//			computeStats(BER[i].begin( ), BER[i].end( ), sum[i], mean[i], var[i], dev[i], skew[i], kurt[i]);
//			sort(BER[i].begin( ), BER[i].end( ));
//		}
		
//		cout << "BER\t";
		for (int i = 0; i<STEPS; i++) {
			cout <<  BER[i] / (double)(K*Run)<< '\t';
		}
		cout << '\n';
//		cout << "BER 90%\t";
//		for (int i = 0; i<STEPS; i++) {
//			cout <<  BER[i][Run/10*9]-mean[i]<< '\t';
//		}
//		cout << '\n';
//		cout << "BER 10%\t";
//		for (int i = 0; i<STEPS; i++) {
//			cout <<  mean[i]-BER[i][Run/10]<< '\t';
//		}
//		cout << '\n';
//		cout << "FER\t";
//		for (int i = 0; i<STEPS; i++) {
//			cout <<  FER[i]/Run<< '\t';
//		}
//		cout << '\n';
//		cout << "Variance\t";
//		for (int i = 0; i<STEPS; i++) {
//			cout <<  var[i]<< '\t';
//		}
//		cout << '\n';
//		cout << "standard deviation\t";
//		for (int i = 0; i<STEPS; i++) {
//			cout <<  dev[i]<< '\t';
//		}
//		cout << '\n';
//		cout << "Skewness\t";
//		for (int i = 0; i<STEPS; i++) {
//			cout <<  skew[i]<< '\t';
//		}
//		cout << '\n';
//		cout << "kurtosis\t";
//		for (int i = 0; i<STEPS; i++) {
//			cout <<  kurt[i]<< '\t';
//		}
//		cout << '\n';
//		
//		cout << "Time\t" << time(0) - start_time<< "\tK\t"<< K << "\tRun\t"<< Run <<endl;
	}
	
	return 0;
}




