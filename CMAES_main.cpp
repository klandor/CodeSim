#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <ctime>    
#include <cstdlib>
#include <math.h>
#include <iomanip>
#include "cmaes_interface.h"
#include "randomc.h"
#include "LT.h"
#include <omp.h>

// =================================================================
#define K 1000			// K size
#define MaxN 1050		// set Max code word number (set 2*K ,but we only use 1.2*K)
#define Run 1000			// how many simulations per fitness 
#define MAXFEC 10000		// set Max function evaluations in CMAES  
#define Lambda 10		// set parameter lambda in CMAES
#define INFO 1			// 1 : show the info during evolution , 0 : don't display
// =================================================================
using namespace std;
using namespace CodeSim;

int M_encode[MaxN][K];				// 記錄 完整的 matrix 
int V_recover[K];					// 記錄 哪些 symbols 已經被 recover
int S_recover = K;					// 記錄 有幾個 bit 已解碼  = sum(V_recover);
int V_degree[MaxN];					// 記錄 M_encode 每個 code word 還有多少 degrees
int V_ripp[K];						// 記錄 還沒被使用的 degree one symbol    
int S_ripp = 0;						// 記錄 V_ripp 的個數 = nnz(V_ripp);    

int M_code2sym[MaxN][K];			// 記錄每個 code 哪些 bit 有值
int S_code2sym[MaxN];
int M_sym2code[K][MaxN];            // 記錄反向連結 - 每個symbol連到哪幾個 codeword
int S_sym2code[K];                  // 反向連結的 size 

// =================================================================
// 依照要跑的 degree 去設定 
int 	Dsize = 10;
int		Set_tags[10] = {1, 2, 3, 4, 5, 7, 9, 19, 59, 179};
// =================================================================

int*	Tags;
double D[] = {0.122890589, 0.290453339, 0.320847326, 0.031809185, 0.1677545, 0.001198199, 
	0.046811719, 0.01644817, 0.001234801, 0.000552171};
double* SD;
double* Std;


CRandomMersenne *Rnd;

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


//  初始化設定參數  從uniform distribution 開始  STD 設為 0.025 
void Parameter_init(){
	Tags = new int[Dsize];
	//D  = new double[Dsize];
	
	SD  = new double[Dsize];
	Std    = new double[Dsize];
	for(int i =0;i<Dsize;i++){
		Tags[i] = Set_tags[i];
		//D[i] = 1/(double)Dsize;
		Std[i] 	 = 0.025;
	}
}




double fitfun(double* Indiv , int dim, bool &needResample){
	//int i;		
	double n=0;
	
	normolize(Indiv);
//	for(i=0;i<Dsize;i++) D[i] = Indiv[i];
//	// SD = 累計機率	
//	SD[0] = D[0];	
//	for(i=1;i<Dsize;i++)	SD[i] = SD[i-1]+D[i];
	int seed[Run];
	for (int i=0; i<Run; i++) {
		seed[i] = Rnd->BRandom();
	}
	
	
	#pragma omp parallel for num_threads(4) reduction(+:n)
	for(int i=0;i<Run;i++){
		LT_sim<Bit> lt(K, MaxN, Dsize, Set_tags, Indiv, seed[i]);
		lt.reset();
		lt.seqReceive(MaxN);
		lt.decode();
		//cout << "Run "<< i+1 << endl;
		n += lt.falureRate();
	}
	return n/Run;
}

/* the optimization loop */
int main(int argn, char **args) {	
	
	int i; 
	fstream fs;	
	Rnd = new CRandomMersenne(time(0));
  	time_t rawtime;
  	struct tm * timeinfo;
	
	cmaes_t* evo; /* an CMA-ES type struct or "object" */
	double *arFunvals, *const*pop, *xbest;
	
	// open file
	fs.open("LT_CMAES_Result.txt",fstream::out);
	// recored start time 
	time(&rawtime);
	timeinfo=localtime( &rawtime );
	fs<<"Start time : "<<asctime(timeinfo)<<endl;	
	
	Parameter_init();
	// write Tags and init distribtuion into file
	fs<<"Tags\n";
	for(i=0;i<Dsize;i++) fs<<Tags[i]<<"\t";
	fs<<"\nInitial distribution \n";
	for(i=0;i<Dsize;i++) fs<<D[i]<<"\t";
	fs<<"\nGen\tFEvals\tFitness\tFbest\tXbest\n";
	
	evo = new cmaes_t();
	
	/* Initialize everything into the struct evo, 0 means default */
	arFunvals = cmaes_init(evo, Dsize, D, Std, 0, Lambda, "non"); 
	evo->sp.stopMaxFunEvals = MAXFEC;	
	cout<<cmaes_SayHello(evo)<<endl;
	
	/* Iterate until stop criterion holds */
	while(!cmaes_TestForTermination(evo))
	{ 						
		/* generate lambda new search points, sample population */
		pop = cmaes_SamplePopulation(evo); /* do not change content of pop */
		/* evaluate the new search points using fitfun from above */ 
		for (i = 0; i < Lambda; ++i) {
			bool needResample = false;
			arFunvals[i] = fitfun(pop[i], Dsize, needResample);
			if (needResample) {
				pop = cmaes_ReSampleSingle(evo, i);
				i--;
			}
		}
		/* update the search distribution used for cmaes_SampleDistribution() */
		cmaes_UpdateDistribution(evo, arFunvals);  
		/* read instructions for printing output or changing termination conditions */ 
		xbest = normolize(cmaes_GetNew(evo, "xbest"));
		fs<<cmaes_Get(evo, "iteration")<<"\t"<<cmaes_Get(evo, "eval")<<"\t"<<cmaes_Get(evo, "fitness")<<"\t"<<cmaes_Get(evo, "fbestever")<<"\t";
		fs.setf(ios::fixed);
		fs.precision(6);
		for(i=0;i<Dsize;i++) 
			fs<<setw(8)<<xbest[i]<<"\t";
		fs.unsetf(ios::fixed);
		fs<<endl;
		
		if(INFO==1) cout<<cmaes_Get(evo, "iteration")<<"\t"<<cmaes_Get(evo, "eval")<<"\t"<<cmaes_Get(evo, "fitness")<<"\t"<<cmaes_Get(evo, "fbestever")<<endl;
		//fflush(stdout); /* useful in MinGW */
		
		delete [] xbest;
	}
	
	printf("Stop:\n%s\n",  cmaes_TestForTermination(evo)); /* print termination reason */
	cmaes_WriteToFile(evo, "all", "allcmaes.dat");         /* write final results */
	
	// recored end time 		
	time(&rawtime);
	timeinfo=localtime( &rawtime );
	fs<<"\nStop time : "<<asctime(timeinfo)<<endl;
	fs.close();	
	
	cmaes_exit(evo); /* release memory */ 		
	delete Tags;
	//delete D;
	delete SD;
	delete Std;	
	delete(evo);
	delete(Rnd);
	system("pause");
	return 0;
}

