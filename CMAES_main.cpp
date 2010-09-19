/*
 *  main_burst_steko.cpp
 *  LT_CMAES
 *
 *  Created by 刁培倫 on 2010/7/21.
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
 */

#include "cmaes_interface.h"
#include "randomc.h"
#include <stdlib.h>
#include <stdio.h>
#include <vector>
#include <iostream>
#include <fstream>
#include <ctime>    
#include <cstdlib>
#include <math.h>
#include <iomanip>
#include <time.h>
#include <climits>
#include <cmath>
#include <string>

#include <omp.h>
#include "LT.h"

// =================================================================
#define K 1000			// K size
#define MaxN 10500		// set Max code word number (set 2*K ,but we only use 1.2*K)
#define Run 12			// how many simulations per fitness 
#define MAXFEC 10000		// set Max function evaluations in CMAES  
#define Lambda 10		// set parameter lambda in CMAES
#define INFO 1			// 1 : show the info during evolution , 0 : don't display
// =================================================================
using namespace std;
using namespace CodeSim;

int g_seed = (int)time(0);
CRandomMersenne RanGen(g_seed);

// =================================================================
// 依照要跑的 degree 去設定 
int 	Dsize = 10;
int		Set_tags[10] = {1, 2, 3, 4, 5, 7, 9, 19, 59, 179};
// =================================================================

int*	Tags;
double D[] = {7.9379E-02,4.0129E-01,1.0121E-01,2.1679E-01,5.0996E-02,
	5.8338E-05,3.9740E-02,7.7470E-02,2.1520E-02,1.1547E-02};
double* SD;
double* Std;


#define epsilonIndex 5
double epsilons[epsilonIndex] = {0.05, 0.07, 0.098, 0.1372, 0.1921},//{0.05, 0.06,0.08, 0.12,0.20},
errorRateBound[epsilonIndex]={4,4,4,4,4},
epsilonBurstBound[epsilonIndex] = {0.5,0.4,0.3,0.2,0.1},
errorDensityBound[epsilonIndex] = {0.33333, 0.26666, 0.2, 0.166666, 0.133333},
areaWeight[epsilonIndex] = {0, 500, 1000, 2000, 4000};

int winSize[epsilonIndex] = {30, 30, 30, 30, 30};
//double failurePenalty[epsilonIndex] = {20, 40, 60, 80, 100};







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

double weighting(double epsilon, double f_rate){
	return epsilon*pow(f_rate,1.5);
	
}


inline double exceed_penalty(double value, double base, double penalty_ratio)
{
	if (value <= base) {
		return 0;
	}
	
	return (value - base) * penalty_ratio;
}

double fitfun(double* Indiv , int dim, bool &needResample){

	normolize(Indiv);
	
	double fit=0, err[epsilonIndex];
	double failureCount[epsilonIndex];
	for (int i=0; i<epsilonIndex; i++) {
		err[i]=0;
		failureCount[i] = 0;
	}
	
	
	#pragma omp parallel for num_threads(4) reduction(+:fit)
	for(int i=0;i<Run;i++){
		//cout << "Run "<< i+1 << endl;
		Codeword<Bit> decodePattern[epsilonIndex];
		
		for (int j=0; j<100; j++) {
			LT_sim<Bit> sim(K, (int) (K*(1+epsilons[epsilonIndex-1])), Dsize, Set_tags, Indiv, RanGen.BRandom());
			
			for (int i=0; i<epsilonIndex; i++) {
				sim.seqReceive(K*(1+epsilons[i])-1);
				sim.decode();
				double temp = sim.failureRate();
				#pragma omp atomic
				err[i] += temp;
				
				if (temp > epsilonBurstBound[i]) {
					#pragma omp atomic
					failureCount[i] += 1;
				}
				Codeword<Bit> t = sim.getResult();
				
				decodePattern[i].insert(decodePattern[i].end(), t.begin(), t.end());
				
			}
		}
		
		
		for (int i=0; i<epsilonIndex; i++) {
			
			int errNO=0, errLen=0;
			for (int p=0; p<winSize[i]; p++) {
				if (decodePattern[i][p].isErased()) {
					errNO ++;
				}
			}
			if(errNO/(double)winSize[i] > errorDensityBound[i])
				errLen=1;
			
			for (int p=winSize[i]; p< 100*K; p++) {
				if (decodePattern[i][p].isErased()) {
					errNO++;
				}
				if (decodePattern[i][p-winSize[i]].isErased()) {
					errNO --;
				}
				
				if(errNO/(double)winSize[i] > errorDensityBound[i])
				{
					errLen ++;
				}
				else {
					if (errLen > 750) {
						//fit +=failurePenalty[i];
						#pragma omp atomic
						failureCount[i] += errLen / 750.0;
					}
					errLen = 0;
				}
				
				
			}	
			
			
		}
		
	}
	
	
	
	for (int i=1; i<epsilonIndex; i++) {
		if (failureCount[i] > 80) {
			//fit +=	200;
			needResample = true;
		}
		if (err[i] > 0) {
			err[i] /= Run*100;
			err[i] = log10(err[i])+4;
			//err[i] += exceed_penalty(err[i], errorRateBound[i], 1000);
			fit += (err[i]+err[i-1]) * (epsilons[i]-epsilons[i-1]) /2 * areaWeight[i];
		}
	}
	
	
	
	//delete [] err;
	return fit;
}

/* the optimization loop */
int main(int argn, char **args) {	
	
	int i; 
	fstream fs;	
	//Rnd = new ran0(0);
  	time_t rawtime;
  	struct tm * timeinfo;
	
	cmaes_t* evo; /* an CMA-ES type struct or "object" */
	double *arFunvals, *const*pop, *xbest;
	
	// open file
	cout << "Enter filename: ";
	string comm;
	getline(cin, comm);
	fs.open(comm.c_str(),fstream::out);
	
	// recored start time 
	time(&rawtime);
	timeinfo=localtime( &rawtime );
	fs<<"Start time : "<<asctime(timeinfo)<<endl;	
	
	cout << "Enter Comment: ";
	getline(cin, comm);
	
	Parameter_init();
	// write Tags and init distribtuion into file
	fs << "Comment: " << comm << '\n';
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
				cout << "R";
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
	//delete(Rnd);
	system("pause");
	return 0;
}


