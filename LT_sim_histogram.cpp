#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <ctime>    
#include <cstdlib>
#include <vector>
#include <math.h>
//#include "yp_random.h"
#include "randomc.h"
#include "LT.h"
#include "statistics.h"
#include <omp.h>
#include <algorithm>

#define K 1000

#define Run 10000000
//#define C 0.05
#define Delta 0.03
#define STEPS 16
#define MaxN (K*(1+Delta*(STEPS-1)))

using namespace std;
using namespace CodeSim;

//int M_encode[MaxN][K];				// 記錄 完整的 matrix 
//int V_recover[K];					// 記錄 哪些 symbols 已經被 recover
//int S_recover = K;					// 記錄 有幾個 bit 已解碼  = sum(V_recover);
//int V_degree[MaxN];					// 記錄 M_encode 每個 code word 還有多少 degrees
//int V_ripp[K];						// 記錄 還沒被使用的 degree one symbol    
//int S_ripp = 0;						// 記錄 V_ripp 的個數 = nnz(V_ripp);    
//
//int M_code2sym[MaxN][K];			// 記錄每個 code 哪些 bit 有值
//int S_code2sym[MaxN];
//int M_sym2code[K][MaxN];            // 記錄反向連結 - 每個symbol連到哪幾個 codeword
//int S_sym2code[K];                  // 反向連結的 size 

unsigned long ErrorCount[STEPS][16];
//double BER[STEPS];
vector< vector<double> > BER;
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

double* normolize(double* d){
	int i;
	double z = 0;
	for(i=0;i<Dsize;i++){
		if(d[i]<0) d[i] = -d[i];
		z = z + d[i];
	}
	for(i=0;i<Dsize;i++) d[i] = d[i]/z;
	return d;
}


int main(){
	//int i;
	//Rnd = new ran3();
	
	//D = Robust_Soliton_Distribution(K,C,Delta);
	
	//cout << "haha";
	cin >> Dsize;
	Degree = new int[Dsize];
	D = new double[Dsize];
	SD = new double[Dsize];
	for (int i=0; i<Dsize; i++) {
		cin >> Degree[i];
	}
	for (int i=0; i<Dsize; i++) {
		cin >> D[i];
	}	
	
	BER.assign(STEPS, 0);
	vector<double> sum(STEPS,0), mean(STEPS,0), var(STEPS,0),
					dev(STEPS,0), skew(STEPS,0), kurt(STEPS,0), FER(STEPS,0);
	//while (cin >> D[0])
	{
		cout << "distribution\t";
		for (int i =0; i<Dsize; i++) {
			//cin >> D[i];
			cout << D[i] << '\t';
			
		}
		cout << "\nEpsilons\t";
		for (int i = 0; i<STEPS; i++) {
			cout << i*Delta << '\t';
			//BER[i] = 0;
			BER[i].assign(Run,0);
			for(int j = 0; j< 16; j++)
			ErrorCount[i][j] = 0;
		}
		cout << '\n';
		
		int start_time = time(0);
		//int n=0;
		#pragma omp parallel for num_threads(6)
		for(int run=0;run<Run;run++){
			LT_sim<Bit> sim(K, MaxN, Dsize, Degree, D, Rnd.BRandom());
			//cout << "Run " << i <<'\n';
			//#pragma omp parallel for
			for (int i = 0; i< STEPS; i++) {
				sim.seqReceive( K*(1+Delta*i) -1);
				//sim.decode();
				double t = sim.failureRate();// = Encoder(K, K*(1.05+0.01*i), Dsize);
				#pragma omp atomic
				ErrorCount[i][(int)(t*16)]++;
				
				BER[i][run]=t;
				if (t>0) {
					FER[i]++;
				}
			}
			
			//cout << '\n';
		}
		//cout<<"Histogram";
		for (int i = 0; i<16; i++) {
			cout << (i+1)/16.0<<'\t';
			for(int j = 0; j< STEPS; j++)
				cout <<  ErrorCount[j][i] / (double)Run<< '\t';
			cout << '\n';
		}
		
		#pragma omp parallel for num_threads(6)
		for (int i = 0; i<STEPS; i++) {
			computeStats(BER[i].begin( ), BER[i].end( ), sum[i], mean[i], var[i], dev[i], skew[i], kurt[i]);
			sort(BER[i].begin( ), BER[i].end( ));
		}
		
		cout << "BER\t";
		for (int i = 0; i<STEPS; i++) {
			cout <<  mean[i]<< '\t';
		}
		cout << '\n';
		cout << "BER 90%\t";
		for (int i = 0; i<STEPS; i++) {
			cout <<  BER[i][Run/10*9]-mean[i]<< '\t';
		}
		cout << '\n';
		cout << "BER 10%\t";
		for (int i = 0; i<STEPS; i++) {
			cout <<  mean[i]-BER[i][Run/10]<< '\t';
		}
		cout << '\n';
		cout << "FER\t";
		for (int i = 0; i<STEPS; i++) {
			cout <<  FER[i]/Run<< '\t';
		}
		cout << '\n';
		cout << "Variance\t";
		for (int i = 0; i<STEPS; i++) {
			cout <<  var[i]<< '\t';
		}
		cout << '\n';
		cout << "standard deviation\t";
		for (int i = 0; i<STEPS; i++) {
			cout <<  dev[i]<< '\t';
		}
		cout << '\n';
		cout << "Skewness\t";
		for (int i = 0; i<STEPS; i++) {
			cout <<  skew[i]<< '\t';
		}
		cout << '\n';
		cout << "kurtosis\t";
		for (int i = 0; i<STEPS; i++) {
			cout <<  kurt[i]<< '\t';
		}
		cout << '\n';
		
		cout << "Time\t" << time(0) - start_time<< "\tK\t"<< K << "\tRun\t"<< Run <<endl;
	}
	
	return 0;
}




