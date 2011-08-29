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

//#define K 1000
int K;
//#define Run (10000000000LL/K)
long Run;
//#define C 0.05
//#define Delta 0.03
double Delta;
//#define STEPS 16
int STEPS;
#define MaxN (K*(1+Delta*(STEPS-1)))
//int windowSize = 50;
long histo_bins[16] = {0, 1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096, 8192, 16384};
const int N_histo_bins = 16;

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

//unsigned long ErrorCount[STEPS][16];
vector< vector<double> > histoErrorCount;

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
	
	//cin >> K >> Run;
	cin >> K;
	Run = 1000000;//(1000000000LL/K);
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
	//cin >> STEPS >> Delta;
	STEPS = 101;
	Delta = 0.005;
	histoErrorCount.assign(STEPS, vector<double>() );
	BER.assign(STEPS, vector<double>());
	vector<double> sum(STEPS,0), mean(STEPS,0), var(STEPS,0),
					dev(STEPS,0), skew(STEPS,0), kurt(STEPS,0);
	
	double rhos[5] = {0, 0.0001, 0.001, 0.01, 0.1}, percents[5] = {.9, .99, .999, .9999, .1};
	const int N_rhos = 5, N_percents = 5;
	vector< vector<long> > BFailureCount(N_rhos, vector<long>(STEPS, 0));
	vector<double> averageEpsilon(N_rhos, 0);
	vector< vector<double> > rhoEpsilons(N_rhos, vector<double>(Run, 0));
	//while (cin >> D[0])
	{
		cout << "tag\t";
		for (int i =0; i<Dsize; i++) {
			//cin >> D[i];
			cout << Degree[i] << '\t';
			
		}
		cout << "\ndistribution\t";
		for (int i =0; i<Dsize; i++) {
			//cin >> D[i];
			cout << D[i] << '\t';
			
		}
		cout << "\nEpsilons\t";
		for (int i = 0; i<STEPS; i++) {
			cout << i*Delta << '\t';
			//BER[i] = 0;
			histoErrorCount[i].assign(N_histo_bins,0);
			BER[i].assign(Run,0);	
		}
		cout << '\n';
		
		if(histo_bins[N_histo_bins-1]<K)
			histo_bins[N_histo_bins-1]=K;
		double histo_bins_ratio[N_histo_bins];
		for (int i=0; i<N_histo_bins; i++) {
			histo_bins_ratio[i] = histo_bins[i]/double(K);
			if (histo_bins_ratio[i]>1) {
				histo_bins_ratio[i]=1;
			}
		}
		
		
		int start_time = time(0);
		//int n=0;
		#pragma omp parallel for schedule(dynamic) num_threads(PARALLEL_THREADS)
		for(long run=0;run<Run;run++){
			int seed;
			#pragma omp critical
			seed = Rnd.BRandom();
			
			LT_sim<Bit> sim(K, 10*K, Dsize, Degree, D, seed);

			long recievedSymbol = 0;
			vector<bool> thresholdReached(N_rhos, 0);
			
			for (int i = 0; i< STEPS; i++) {
				while (recievedSymbol <= K*(1+Delta*i) -.9) {
					sim.receive(recievedSymbol);
					recievedSymbol++;
					
					double t = sim.failureRate();
					for (int r=0; r<N_rhos; r++) {
						if (thresholdReached[r] == false && t<=rhos[r]+(1.0/K/10)) {
							#pragma omp atomic
							averageEpsilon[r] += (recievedSymbol-K)/(double)K/Run;
							rhoEpsilons[r][run] = (recievedSymbol-K)/(double)K; 
							thresholdReached[r] = true;
						}
					}
				}
				//sim.seqReceive( K*(1+Delta*i) -1);
				//sim.decode();
				double t = sim.failureRate();// = Encoder(K, K*(1.05+0.01*i), Dsize);

				for (int h=0; h<N_histo_bins; h++) {
					if (t<=histo_bins_ratio[h]+(1.0/K/10)) {
						#pragma omp atomic
						histoErrorCount[i][h]++;
						break;
					}
				}

				BER[i][run]=t;
				for (int j=0; j<N_rhos; j++) {
					if (t>rhos[j]+(1.0/K/10)) {
						#pragma omp atomic
						BFailureCount[j][i]++;
					}
				}

			}
			
			bool allReached = false;
			while (allReached == false) {
				allReached = true;
				sim.receive(recievedSymbol);
				recievedSymbol++;
				double t = sim.failureRate();
				for (int r=0; r<N_rhos; r++) {
					if (thresholdReached[r] == false && t<=rhos[r]+(1.0/K/10)) {
						#pragma omp atomic
						averageEpsilon[r] += (recievedSymbol-K)/(double)K/Run;
						rhoEpsilons[r][run] = (recievedSymbol-K)/(double)K;
						thresholdReached[r] = true;
					}
					
					allReached &= thresholdReached[r];
				}
			}
		}
		//cout<<"Histogram";
		for (int i = 0; i<N_histo_bins; i++) {
			cout << histo_bins_ratio[i]<<'\t';
			for(int j = 0; j< STEPS; j++)
				cout <<  histoErrorCount[j][i] / (double)Run<< '\t';
			cout << '\n';
		}
		
		#pragma omp parallel for schedule(dynamic) num_threads(PARALLEL_THREADS)
		for (int i = 0; i<STEPS; i++) {
			computeStats(BER[i].begin( ), BER[i].end( ), sum[i], mean[i], var[i], dev[i], skew[i], kurt[i]);
			sort(BER[i].begin( ), BER[i].end( ));
		}
		
		cout << "AverageFailureRatio\t";
		for (int i = 0; i<STEPS; i++) {
			cout <<  mean[i]<< '\t';
		}
		cout << '\n';
		for (int j=0; j<N_percents; j++) {
			cout << "FailureRatio @"<< percents[j]*100 <<"%\t";
			for (int i = 0; i<STEPS; i++) {
				cout <<  BER[i][Run*percents[j]-0.9]<< '\t';
			}
			cout << '\n';
		}
		
		for (int j=0; j<N_rhos; j++) {
			cout << "BlockFailureRate rho="<< rhos[j]*100 <<"%\t";
			for (int i = 0; i<STEPS; i++) {
				cout <<  BFailureCount[j][i]/(double)Run<< '\t';
			}
			cout << '\n';
		}
		
		for (int j=0; j<N_rhos; j++) {
			cout << "pdf rho="<< rhos[j]*100 <<"%\t";
			cout <<  1-(BFailureCount[j][0]/(double)Run)<< '\t';
			for (int i = 1; i<STEPS; i++) {
				cout <<  (BFailureCount[j][i-1]-BFailureCount[j][i])/(double)Run<< '\t';
			}
			cout << '\n';
		}
		
		for (int j=0; j<N_rhos; j++) {
			cout << "rho="<< rhos[j]*100 <<"%\t";
			cout <<  averageEpsilon[j]<< '\t';
			sort(rhoEpsilons[j].begin(), rhoEpsilons[j].end());
			cout << rhoEpsilons[j][Run*0.9-0.9] - averageEpsilon[j]<< '\t'; // 90% error bar
			cout << averageEpsilon[j] - rhoEpsilons[j][Run*0.1-0.9]<< '\t'; // 10% error bar
			cout << '\n';
		}
		
		
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
		
//		for (int i=0; i<windowSize+1; i++) {
//			cout << i;
//			
//			for (int j=0; j<STEPS; j++) {
//				cout << '\t' << l0a[j][i]/(double)(Run*(K-windowSize+1));
//			}
//			cout << endl;
//		}
	}
	
	return 0;
}




