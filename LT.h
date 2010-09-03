/*
 *  LT.h
 *  new_LT
 *
 *  Created by 刁培倫 on 2010/8/6.
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
 */

#include <vector>
#include "randomc.h"
#include "CodeSim.h"
using namespace std;

namespace CodeSim{

	template<class S>
	class LT_sim : public CodingBlock<S,S>{
	public:
		LT_sim();
		LT_sim(int k, int max_n, int tag_size, int *tags, double * omega, int seed);
		void init(int k, int max_n, int tag_size, int *tags, double * omega, int seed);
		void seqReceive(int t);
		void receive(int t);
		void decode();
		bool isDecoded(int t);
		double run();
		Codeword<S> encode(Codeword<S> a);
		Codeword<S> decode(Codeword<S> a);
		void reset();
	private:
		void codeGen(int t);
		
		int Num_of_Input, Num_of_Output, Num_of_Degree, Num_of_Decoding, 
				ReceivedSize, generatedCode;
		int seed;
		CRandomMersenne Rnd;
		vector<long int>	d;//(Num_of_Output, 0);
		vector<vector<long int> >	edge;// = new vector<long int>[Num_of_Output];
		vector<long int>	R_M;//(Num_of_Output, 0);
		//vector<long int>	erasure;//(Num_of_Output, 0);
		Codeword<S>	DE;//(Num_of_Input, 0);
		vector<int>			Degree_of_Edge;
		vector<double>		Omega;
		vector<int>			receivedMask;
	};
	
	template<class S>
	LT_sim<S>::LT_sim():Rnd(0){}
	
	template<class S>
	LT_sim<S>::LT_sim(int k, int max_n, int tag_size, int *tags, double * omega, int seed):Rnd(seed){
		init(k, max_n, tag_size, tags, omega, seed);
	}
	template<class S>
	void LT_sim<S>::init(int k, int max_n, int tag_size, int *tags, double * omega, int seed){
		Num_of_Input = k;
		Num_of_Output= max_n;
		Num_of_Degree = tag_size;
		
		d.assign(max_n,0);
		edge.assign(max_n,0);
		R_M.assign(max_n,0);
		//erasure.assign(max_n,0);
		
		DE.assign(k,0);
		this->seed = seed;
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
	
	template<class S>
	bool LT_sim<S>::isDecoded(int t){
		return DE[t];
	}
	
	template<class S>
	void LT_sim<S>::codeGen(int t){
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
	
	template<class S>
	inline void LT_sim<S>::receive(int t){
		
		if(receivedMask[t]==0)
		{
			R_M[ReceivedSize]=t;
			ReceivedSize++;
			receivedMask[t]=1;
		}
		
	}
	template<class S>
	void LT_sim<S>::seqReceive(int t){
		for(int i=0; i <= t; i++)
		{
			receive(i);
		}
	}
	
	template<class S>
	void LT_sim<S>::decode(){
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
	
	template<class S>
	Codeword<S> LT_sim<S>::encode(Codeword<S> a){
		return Codeword<S>(Num_of_Output, 0);
	}
	
	template<class S>
	Codeword<S> LT_sim<S>::decode(Codeword<S> a){
		for (int i =0; i< a.size() && i< Num_of_Output; i++) {
			if (!a[i].isErased() ) {
				receive(i);
			}
		}
		decode();
		
		Codeword<S> t = DE;
		for (int i=0; i< t.size(); i++) {
			if (t[i] == 0) {
				t[i] = -1;
			}
		}
		return t;
	}
	
	template<class S>
	void LT_sim<S>::reset(){
		//	Num_of_Input = k;
		//	Num_of_Output= max_n;
		//	Num_of_Degree = tag_size;
		
		d.assign(Num_of_Output,0);
		edge.assign(Num_of_Output,0);
		R_M.assign(Num_of_Output,0);
		//erasure.assign(max_n,0);
		
		DE.assign(Num_of_Input,0);
		Rnd.RandomInit(seed);
		
		
		Num_of_Decoding = 0;
		ReceivedSize = 0;
		generatedCode = 0;
		codeGen(Num_of_Output);
		
		receivedMask.assign(Num_of_Output, 0);
	}
	
	template<class S>
	double LT_sim<S>::run(){
		
		return 1-(Num_of_Decoding/(double)Num_of_Input);//failure rate
		
	}
}
