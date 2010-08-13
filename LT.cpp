/*
 *  LT.cpp
 *  new_LT
 *
 *  Created by 刁培倫 on 2010/8/6.
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
 */

#include "LT.h"
using namespace LT_class;

LT_sim::LT_sim():Rnd(0){}

LT_sim::LT_sim(int k, int max_n, int tag_size, int *tags, double * omega, int seed):Rnd(seed){
	init(k, max_n, tag_size, tags, omega, seed);
}

void LT_sim::init(int k, int max_n, int tag_size, int *tags, double * omega, int seed){
	Num_of_Input = k;
	Num_of_Output= max_n;
	Num_of_Degree = tag_size;
	
	d.assign(max_n,0);
	edge.assign(max_n,0);
	R_M.assign(max_n,0);
	//erasure.assign(max_n,0);
	
	DE.assign(k,0);
	Rnd.RandomInit(seed);
	Degree_of_Edge.assign(tags, tags+tag_size);
	Omega.assign(omega, omega+tag_size);
	for(int i=1; i<Num_of_Degree; i++)  
		Omega[i] = Omega[i-1]+Omega[i];
	
	Num_of_Decoding = 0;
	ReceivedSize = 0;
	generatedCode = 0;
	codeGen(max_n);
	
	receivedMask.assign(max_n, 0);
}

bool LT_sim::isDecoded(int t){
	return DE[t];
}

void LT_sim::codeGen(int t){
	if (t>= generatedCode){
		for(int i=generatedCode; i<t+1; i++){
			// decide degree
			double rando = Rnd.Random();
			int s, flag;
			
			if(rando<Omega[0])
			{
				d[i] = Degree_of_Edge[0];
			}      
			
			else
			{
				for(s=0;s<Num_of_Degree-1;s++)
				{
					if(rando>=Omega[s] && rando<Omega[s+1])
					{
						d[i] = Degree_of_Edge[s+1];
						break;
					}
				}
			}

			// decide connection
			edge[i].assign( d[i], 0 );
			for(int j=0;j<d[i];j++)
			{
				flag=1;
				while(flag==1)
				{
					flag=0;
					//	edge[i][j]=rand()%Num_of_Input;
					edge[i][j]=Rnd.IRandomX(0,Num_of_Input-1);
					for(int lllll=0;lllll<j;lllll++)
					{
						if(edge[i][j]==edge[i][lllll])
						{
							flag=1;
						}
					}
				}
			}
			
			
		}
		
		generatedCode = t;
	}
}

inline void LT_sim::receive(int t){

	if(receivedMask[t]==0)
	{
		R_M[ReceivedSize]=t;
		ReceivedSize++;
		receivedMask[t]=1;
	}

}

void LT_sim::seqReceive(int t){
	for(int i=0; i <= t; i++)
	{
		receive(i);
	}
}

void LT_sim::decode(){
	int flag=1;
	while( Num_of_Decoding < Num_of_Input && flag==1)
	{
		flag=0;
		for(int i=0;i<ReceivedSize;i++)
		{
			if(d[R_M[i]]==1 && DE[ edge[ R_M[i] ][0] ]==0)
			{
				DE[ edge[ R_M[i] ][0] ]=1;
				//    			Result[edge[R_M[i]][0]] = Mid_Output[R_M[i]];
				Num_of_Decoding++;
				ReceivedSize--;
				R_M[i]=R_M[ReceivedSize];
				i--;
				flag=1;
			}
		}
		
		for(int i=0;i<ReceivedSize;i++)
		{
			for(int j=0;j<d[R_M[i]];j++)
			{
				if(DE[ edge[ R_M[i] ][j] ]==1)
				{
					//		            Mid_Output[R_M[i]] = Mid_Output[R_M[i]] + Result[edge[R_M[i]][j]];
                    d[R_M[i]]--;
					edge[ R_M[i] ][j]=edge[ R_M[i] ][ d[R_M[i]] ];
					j--;
				}
			}
		}
	} 
}

vector<char> LT_sim::encode(vector<char> a){
	return a;
}

vector<char> LT_sim::decode(vector<char> a){
	return a;
}

void LT_sim::reset(){
	
}

double LT_sim::run(){
//	long int           Degree_of_Edge[10] ;//= {1, 2, 3, 4, 5, 7, 9, 19, 59, 179};//{1, 2, 3, 4, 5, 8, 20, 50, 100, 120};
//	double			Omega[10];// ={0.095874,	0.310842,	0.347227,	0.126544,	0.001954,	0.015687,	0.022324	,0.074265,	0.003683,	0.001599	};//Geometry Weights(Tail2)
//	
//	for (int i=0; i<10; i++) {
//		Degree_of_Edge[i]=Degree[i];
//		Omega[i] = D[i];
//	}
	//	 double			Omega[10] ={0.058227,	0.472135,	0.119781,	0.190464,	0.039280,	0.015102,	0.050586,	0.024640,	0.027895,	0.001890};//Geometry Weights (with Resampling);
	//	double				Omega[30] ={0.088330251,0.034983829,0.009830305,0.038636479,0.001355699,0.014015336,0.006981217,0.005503599,0.058033691,0.065646432,0.076580531,0.010288606,0.047382095,0.094235979,0.00441917,0.001297473,0.04632022,0.045155864,0.020108309,0.013581224,0.04330058,0.033802575,0.088606895,0.003247101,0.022809331,0.02441986,0.039472951,0.006162829,0.015877843,0.039613726};
	//	 double			Omega[10] = {1, 0, 0, 0, 0, 0, 0, 0, 0, 0};
	//	 	 double				Omega[10]		={0.005586592,0.5,0.166666667,0.083333333,0.05,0.023809524,0.013888889,0.002631579,0.000292,0.153791189};//Soliton
	//	 double				Omega[10] = {0.050728762,0.424347107,0.282085453,0.026033206,0.016145565,0.136520481,0.000641569,0.060746163,0.001547985,0.001203709};//2 seperate region.
	//	 double			Omega[10] = {0.060399311,0.390181264,0.304572196,0.00142503,0.104379941,0.027244498,0.029234722,0.079855937,0.00078359,0.00192351};//Arithmetic weighting
	//	 double			Omega[20] = {0.118489793,0.340930276,0.228353534,0.099050778,0.084985364,0.034303156,0.026062356,0.050101316,0.00058564,0.010044684,0.000734172,0.000461251,0.002187003,0.0000874028,0.000652915,0.00016307,0.00201962,0.000168778,0.000400117,0.000218775};//20tags, Failure Rate@0.05
	//	 double Omega[10] {0.296061	4.43864	2.45673	1.09043	0.470109	0.636914	0.15594	0.826375	0.352833	0.617544}//Mean_Slope
	//double				/*midterm=0,*/ rando= 0;
//	long int			//i, j, lllll,
//				flag,/* Lflag =0, pattern =1,*/ size;//, Num_of_Decoding =0, total=0 ;//, counter = 0;
	//int				rangee =1000;
	
	
	
	//	vector<unsigned int> vec(Num_of_Output,0);
	
	
	//cout<<"Degree Generating..."<<endl;     
	//===========================Degree Generation==================================     
	
	
//	for(int i=0; i<Num_of_Output; i++)
//	{
//		//   midterm = rand()%rangee;
//		
//		//    rando=(double)( midterm / ((double)rangee) );
//		
//		double rando = Rnd.Random();
//		int s;
//		
//		if(rando<Omega[0])
//		{
//			d[i] = Degree_of_Edge[0];
//		}      
//		
//		else
//		{
//			for(s=0;s<Num_of_Degree-1;s++)
//			{
//				if(rando>=Omega[s] && rando<Omega[s+1])
//				{
//					d[i] = Degree_of_Edge[s+1];
//					break;
//				}
//			}
//		}
//		
////		if(Lflag==0) {
////			d[i] = Degree_of_Edge[s+1];
////		}
//		
//	}
	//=====================End of Degree Generation=================================
	
	//cout<<"Connecting Table Generating..."<<endl;
	//====================Connecting Table Generation===============================
	
	
//	for (int i=0; i<Num_of_Output; i++)
//	{
//		//edge[i].reserve(d[i]); 
//		edge[i].assign( d[i], 0 );
//	}
	
//	for(int i=0;i<Num_of_Output;i++)
//	{
//		for(int j=0;j<d[i];j++)
//		{
//			edge[i][j]=0;
//		}
//	}
	
	
//	for(int i=0;i<Num_of_Output;i++)
//	{
//		for(int j=0;j<d[i];j++)
//		{
//			flag=1;
//			while(flag==1)
//			{
//				flag=0;
//				//	edge[i][j]=rand()%Num_of_Input;
//				edge[i][j]=Rnd.IRandomX(0,Num_of_Input-1);
//				for(int lllll=0;lllll<j;lllll++)
//				{
//					if(edge[i][j]==edge[i][lllll])
//					{
//						flag=1;
//					}
//				}
//			}
//		}
//	}
	
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
	
	//i=0;
//	int ReceivedSize=0;
//	
//	
//	for(int i=0; i < Num_of_Output; i++)
//	{
//		if(erasure[i]==0)
//		{
//			//R[ReceivedSize] = Mid_Output[i];
//			R_M[ReceivedSize]=i;
//			ReceivedSize++;
//		}
//	}
	
	
	
	//=================== End of Receiving Codeword ================================
	
	//cout<<"Decoding..."<<endl;
	//========================== Decoding ==========================================
	
//	int flag=1;
//	while( Num_of_Decoding < Num_of_Input && flag==1)
//	{
//		flag=0;
//		for(int i=0;i<ReceivedSize;i++)
//		{
//			if(d[R_M[i]]==1 && DE[ edge[ R_M[i] ][0] ]==0)
//			{
//				DE[ edge[ R_M[i] ][0] ]=1;
//				//    			Result[edge[R_M[i]][0]] = Mid_Output[R_M[i]];
//				Num_of_Decoding++;
//				ReceivedSize--;
//				R_M[i]=R_M[ReceivedSize];
//				i--;
//				flag=1;
//			}
//		}
//		
//		for(int i=0;i<ReceivedSize;i++)
//		{
//			for(int j=0;j<d[R_M[i]];j++)
//			{
//				if(DE[ edge[ R_M[i] ][j] ]==1)
//				{
//					//		            Mid_Output[R_M[i]] = Mid_Output[R_M[i]] + Result[edge[R_M[i]][j]];
//                    d[R_M[i]]--;
//					edge[ R_M[i] ][j]=edge[ R_M[i] ][ d[R_M[i]] ];
//					j--;
//				}
//			}
//		}
//	} 
	
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
	//	fpattern<<endl;
	
//	for(int de=1; de<Num_of_Input; de++){
//		//		conti_zero<<freq_zero[de]<<'\t';
//		//		conti_one<<freq_one[de]<<'\t';
//	}
	
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
	//total=Num_of_Decoding;
	//*Num_of_Decoded = *Num_of_Decoded + (double)(Num_of_Decoding);
	
	//================================================================
	
	//------------------Memory allocated in Encoding Function-----------------------
	
	//delete [] edge;
	//----------------------------------------------------
	//cout<<"Q________Q"<<endl;
	
	return 1-(Num_of_Decoding/(double)Num_of_Input);//failure rate
	//return total;
	
}