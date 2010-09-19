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


double Encoder(long int, long int, long int, short*);

double Encoder(long int Num_of_Input, long int Num_of_Output, long int Num_of_Degree, short* bitmap)
{//國光學長記號 ========    k   ================   N   ===============   d_edge   ===========================================  
	
	
	
	//Codeword *R      = new Codeword[Num_of_Output];
	//Codeword *Result = new Codeword[Num_of_Input];
	int zero_counter= 0;
	int one_counter = 0;
	int path_error = 0;
	
	
	//	vector<int>		freq_one(Num_of_Input,0);
	//	vector<int>		freq_zero(Num_of_Input,0);
	//	vector<int>		freq_path(17,0);
	
	vector<double>		temp(Num_of_Degree, 0);
	vector<long int>	d(Num_of_Output, 0);
	vector<long int>	*edge = new vector<long int>[Num_of_Output];
	vector<long int>	R_M(Num_of_Output, 0);
	vector<long int>	erasure(Num_of_Output, 0);
	vector<long int>	DE(Num_of_Input, 0);
	
	//	ofstream fpattern;
	//	ofstream conti_zero;
	//	ofstream conti_one;
	//	ofstream conti_path;
	//	
	//	conti_zero.open("continous_zero.txt",fstream::out|fstream::app);
	//	conti_one.open("continous_one.txt",fstream::out|fstream::app);
	//	conti_path.open("continous_path.txt",fstream::out|fstream::app);
	//	fpattern.open("failure patterns.txt",fstream::out|fstream::app);
	
	//int           Degree_of_Edge[10] = {1, 2, 3, 4, 5, 7, 9, 19, 59, 179};//{1, 2, 3, 4, 5, 8, 20, 50, 100, 120};
	int           *Degree_of_Edge = Set_tags;
	double			Omega[10];// ={0.095874,	0.310842,	0.347227,	0.126544,	0.001954,	0.015687,	0.022324	,0.074265,	0.003683,	0.001599	};//Geometry Weights(Tail2)
	
	for (int i=0; i<10; i++) {
		Omega[i] = D[i];
	}
	//	 double			Omega[10] ={0.058227,	0.472135,	0.119781,	0.190464,	0.039280,	0.015102,	0.050586,	0.024640,	0.027895,	0.001890};//Geometry Weights (with Resampling);
	//	double				Omega[30] ={0.088330251,0.034983829,0.009830305,0.038636479,0.001355699,0.014015336,0.006981217,0.005503599,0.058033691,0.065646432,0.076580531,0.010288606,0.047382095,0.094235979,0.00441917,0.001297473,0.04632022,0.045155864,0.020108309,0.013581224,0.04330058,0.033802575,0.088606895,0.003247101,0.022809331,0.02441986,0.039472951,0.006162829,0.015877843,0.039613726};
	//	 double			Omega[10] = {1, 0, 0, 0, 0, 0, 0, 0, 0, 0};
	//	 	 double				Omega[10]		={0.005586592,0.5,0.166666667,0.083333333,0.05,0.023809524,0.013888889,0.002631579,0.000292,0.153791189};//Soliton
	//	 double				Omega[10] = {0.050728762,0.424347107,0.282085453,0.026033206,0.016145565,0.136520481,0.000641569,0.060746163,0.001547985,0.001203709};//2 seperate region.
	//	 double			Omega[10] = {0.060399311,0.390181264,0.304572196,0.00142503,0.104379941,0.027244498,0.029234722,0.079855937,0.00078359,0.00192351};//Arithmetic weighting
	//	 double			Omega[20] = {0.118489793,0.340930276,0.228353534,0.099050778,0.084985364,0.034303156,0.026062356,0.050101316,0.00058564,0.010044684,0.000734172,0.000461251,0.002187003,0.0000874028,0.000652915,0.00016307,0.00201962,0.000168778,0.000400117,0.000218775};//20tags, Failure Rate@0.05
	//	 double Omega[10] {0.296061	4.43864	2.45673	1.09043	0.470109	0.636914	0.15594	0.826375	0.352833	0.617544}//Mean_Slope
	double				midterm=0, rando= 0;
	long int			i, j, l, flag, s=0, Lflag =0, pattern =1, size, Num_of_Decoding =0, total=0, counter = 0;
	int				rangee =1000;
	
	
	
	vector<unsigned int> vec(Num_of_Output,0);
	
	
	//cout<<"Degree Generating..."<<endl;     
	//===========================Degree Generation==================================     
	for(i=1; i<Num_of_Degree; i++)  Omega[i] = Omega[i-1]+Omega[i];
	
	for(i=0; i<Num_of_Output; i++)
	{
		//   midterm = rand()%rangee;
		
		//    rando=(double)( midterm / ((double)rangee) );
		
		rando = RanGen.Random();
		
		if(rando<Omega[0])
		{
			Lflag=1;
			d[i] = Degree_of_Edge[0];
		}      
		
		else
		{
			for(s=0;s<Num_of_Degree-1;s++)
			{
				if(rando>=Omega[s] && rando<Omega[s+1])
				{
					Lflag=1;
					d[i] = Degree_of_Edge[s+1];
					break;
				}
			}
		}
		
		if(Lflag==0) {d[i] = Degree_of_Edge[s+1];}
		
		for(j=0;j<Num_of_Degree;j++)
		{
			if(d[i]==Degree_of_Edge[j])
			{
				temp[j]=temp[j]+1;
				break;
			}      
		}
	}
	//=====================End of Degree Generation=================================
	
	//cout<<"Connecting Table Generating..."<<endl;
	//====================Connecting Table Generation===============================
	
	
	for (i=0; i<Num_of_Output; i++)
	{
		edge[i].reserve(d[i]); 
		edge[i].assign( d[i], 0 );
	}
	
	for(i=0;i<Num_of_Output;i++)
	{
		for(j=0;j<d[i];j++)
		{
			edge[i][j]=0;
		}
	}
	
	
	for(i=0;i<Num_of_Output;i++)
	{
		for(j=0;j<d[i];j++)
		{
			flag=1;
			while(flag==1)
			{
				flag=0;
				//	edge[i][j]=rand()%Num_of_Input;
				edge[i][j]=RanGen.IRandomX(0,Num_of_Input-1);
				for(l=0;l<j;l++)
				{
					if(edge[i][j]==edge[i][l])
					{
						flag=1;
					}
				}
			}
		}
	}
	
	//=================End of Connecting Table Generation===========================
	
	//cout<<"Encoding ..."<<endl;
	//===========================Encoding===========================================
	
	//	for(i=0;i<Num_of_Output;i++)
	//	{ 
	//		for(j=0;j<d[i];j++)
	//		{  
	//            Mid_Output[i] = Mid_Output[i] + Input[edge[i][j]];
	//		}
	//	}
	
	//======================= End of Encoding ======================================
	
	//cout<<"Receiving Codeword..."<<endl;
	//======================= Receiving Codeword ===================================
	//                  --這個部分暫時沒用到，過渡用-- 
	
	i=0;
	size=0;
	
	
	while(i < Num_of_Output)
	{
		if(erasure[i]==0)
		{
			//R[size] = Mid_Output[i];
			R_M[size]=i;
			size++;
		}
		i++;
	}
	
	
	
	//=================== End of Receiving Codeword ================================
	
	//cout<<"Decoding..."<<endl;
	//========================== Decoding ==========================================
	
	flag=1;
	while( Num_of_Decoding < Num_of_Input && flag==1)
	{
		flag=0;
		for(i=0;i<size;i++)
		{
			if(d[R_M[i]]==1 && DE[edge[R_M[i]][0]]==0)
			{
				DE[edge[R_M[i]][0]]=1;
				//    			Result[edge[R_M[i]][0]] = Mid_Output[R_M[i]];
				Num_of_Decoding++;
				size--;
				R_M[i]=R_M[size];
				i--;
				flag=1;
			}
		}
		
		for(i=0;i<size;i++)
		{
			for(j=0;j<d[R_M[i]];j++)
			{
				if(DE[edge[R_M[i]][j]]==1)
				{
					//		            Mid_Output[R_M[i]] = Mid_Output[R_M[i]] + Result[edge[R_M[i]][j]];
                    d[R_M[i]]--;
					edge[R_M[i]][j]=edge[R_M[i]][d[R_M[i]]];
					j--;
				}
			}
		}
	} 
	
	//====================== End of Decoding =========================
	
	//	for(int de=0; de<Num_of_Input; de++){
	//		//		fpattern<<DE[de]<<'\t';
	//		if(DE[de]==0){
	//			zero_counter++;
	//		}
	//		else{
	//			freq_zero[zero_counter]++;
	//			zero_counter=0;
	//		}
	//		
	//		if(DE[de]==1){
	//			one_counter++;
	//		}
	//		else{
	//			freq_one[one_counter]++;
	//			one_counter=0;
	//		}
	//	}
	//	//	fpattern<<endl;
	//	
	//	for(int de=1; de<Num_of_Input; de++){
	//		//		conti_zero<<freq_zero[de]<<'\t';
	//		//		conti_one<<freq_one[de]<<'\t';
	//	}
	//	
	//	for(int de=0; de<Num_of_Input-15; de++){
	//		path_error=0;
	//		for(int len_of_fundemantal_path=0;len_of_fundemantal_path<16;len_of_fundemantal_path++){
	//			if(DE[de+len_of_fundemantal_path]==0){
	//				path_error++;
	//			}
	//		}
	//		freq_path[path_error]++;
	//		path_error=0;
	//	}
	//	
	//	for(int de=1; de<=16; de++){
	//		//		conti_path<<freq_path[de]<<'\t';
	//	}
	
	//	conti_path<<endl;
	//	conti_zero<<endl;
	//	conti_one<<endl;
	//===================== Calculator =======================
	
	
	//infile<<Num_of_Decoding<<endl;
	total=Num_of_Decoding;
	//*Num_of_Decoded = *Num_of_Decoded + (double)(Num_of_Decoding);
	
	//================================================================
	
	//------------------Memory allocated in Encoding Function-----------------------
	
	delete [] edge;
	//----------------------------------------------------
	//cout<<"Q________Q"<<endl;
	
	for (int i=0; i< Num_of_Input; i++) {
		bitmap[i] = DE[i];
	}
	return 1-(total/(double)Num_of_Input);
	
}

//============================ End of Decoding =================================






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
	int i;		
	double n=0;
	
	
	normolize(Indiv);
	for(i=0;i<Dsize;i++) D[i] = Indiv[i];
	// SD = 累計機率	
	SD[0] = D[0];	
	for(i=1;i<Dsize;i++)	SD[i] = SD[i-1]+D[i];
	
	double fit=0, err[epsilonIndex];
	double failureCount[epsilonIndex];
	//int n=0;
	for (int i=0; i<epsilonIndex; i++) {
		err[i]=0;
		failureCount[i] = 0;
	}
	
	
#pragma omp parallel for num_threads(4) reduction(+:fit)
	for(i=0;i<Run;i++){
		cout << "Run "<< i+1 << endl;
		//for (int i=0; i<15; i++) {
		//			n += weighting( 0.01+0.005*i, Encoder(K, K*(1.01+0.005*i), Dsize));
		//		}
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
			
			//delete [] bitmap;
			
			
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


