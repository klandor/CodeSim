/*
 *  ConvoCode.cpp
 *  new_LT
 *
 *  Created by 刁培倫 on 2010/9/7.
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
 */

#include "ConvoCode.h"
#include <list>

namespace CodeSim {

	int octToDec(int t){
		if(t <= 0)
			return 0;
		else {
			return ((t%10)&7) ^ (octToDec(t/10) << 3);
		}
		
	}
	
	int setBit(int i, int position, bool value){
		i = i & (-1 - (1 << position) );
		int t = value;
		t <<= position;
		return i+t;
	}
	
	unsigned int countBit(unsigned int i){
		if(i == 0)
			return 0;
		if(i&1)
			return 1 + countBit(i>>1);
		else {
			return countBit(i>>1);
		}
		
	}
	
	
	
	/**
	 *	@brief Construct ConvoCode from file.
	 *	@param filename The file to be read.
	 *
	 *
	 *	The file should consist 2+k lines.
	 *	1st line: N and K
	 *	2nd line: Forney Indices(K numbers)
	 *	3rd to 2+K line: generator matrix (K lines of N numbers, int octal)
	 */
	
	ConvoCode::ConvoCode(string filename){
		ifstream fin(filename.c_str());
		if(fin.fail()){
			cerr << "File opening Error: " << filename << '\n';
		}
		
		fin>>n>>k;
		forneyIndices.assign(k, 0 );
		m=0;
		max_forney = 0;
		for (int i=0; i<k; i++) {
			fin >> forneyIndices[i];
			m += forneyIndices[i];
			if(forneyIndices[i] > max_forney)
				max_forney = forneyIndices[i];
		}
		
		G.assign(k,0);
		for (int i=0; i<k; i++) {
			G[i].assign(n,0);
			for(int j=0; j<n; j++){
				int t;
				fin >> t;
				G[i][j] = octToDec(t);
			}
		}
		
		generateTrellis();
		
	}
	
	void ConvoCode::showInfo(){
		cout << "N: " << n << "\tK: " << k 
		<< "\nForney Indices: ";
		for (int i = 0; i< forneyIndices.size(); i++)
			cout << '\t' <<forneyIndices[i];
		cout << "\nG:\n";
		for(int i = 0; i< k; i++){
			for (int j=0; j< n; j++){
				cout << G[i][j] << ' ';
			}
			cout << '\n';
		}
		cout << "Trellis:\n";
		for(int i=0; i< trellis.size(); i++)
		{
			for(int j=0;j<trellis[i].size(); j++){
				cout << trellis[i][j].out << ' ';
			}
			cout << '\n';
		}
		
	}
	
	void ConvoCode::generateTrellis(){
		trellis.assign(1<<m, 0);
		for(int i=0; i< (1<<m); i++){
			trellis[i].assign(1<<k, state());
			for(int j=0; j< (1<<k); j++){
				// determine next state
				int t = i << 1, p=0, input = j; ;
				for(int d=0; d<k; d++){
					if(forneyIndices[k-1-d] > 0){ 
						t = setBit(t, p, input & 1);
						p += forneyIndices[k-1-d];
					}
					input >>= 1;
				}
				trellis[i][j].next = t & ((1<<m) -1);
				
				// determine output
				t = 0, p=0, input = j;
				for(int d=0; d<k; d++){
					if(forneyIndices[k-1-d] > 0){ 
						t += (( i >> (p-d) ) & ( (1<<forneyIndices[k-1-d]) -1)) << (p+1); // take out the memory of (k-1-d)th input
					}
					t = setBit(t, p, input & 1);
					p += forneyIndices[k-1-d]+1;
					input >>= 1;
				}
				
				
				for(int column = 0; column < n; column++){
					// get column of G
					int g = 0;
					trellis[i][j].out <<= 1;
					for(int row=0; row < k; row++){
						g <<= forneyIndices[row]+1;
						g += G[row][column];
					}
					trellis[i][j].out += countBit(g & t) & 1;
				}
			}
			
		}
	}
	
	
	
	Codeword<Bit> ConvoCode::encode(Codeword<Bit> a){
		Codeword<Bit> output;
		stack<int> s = a.getMessageStack();
		s.push(a.size());
		
		int i, *o = new int[n], m = 0;
		
		for(i=0; i+k<a.size(); i+=k){
			int t=0;
			for(int j=0; j<k; j++){
				t <<= 1;
				t += a[i+j].getValue();
			}
			
			intToArray(o, trellis[m][t].out, n);
			
			output.insert(output.end(), o, o+n);
			
			m = trellis[m][t].next;
			
		}
		// last bits
		int t=0;
		for(int j=0; j< k; j++){
			t <<= 1;
			if(i+j < a.size())
			{
				t += a[i+j].getValue();
			}
			
		}
		
		intToArray(o, trellis[m][t].out, n);
		
		output.insert(output.end(), o, o+n);
		
		m = trellis[m][t].next;
		// end last bits
		
		// todo closing
		for(i = 0; i< max_forney; i++){
			int t=0;
						
			intToArray(o, trellis[m][t].out, n);
			
			output.insert(output.end(), o, o+n);
			
			m = trellis[m][t].next;
		}
		
		delete [] o;
		output.setMessageStack(s);
		return output;
	}
	
	Codeword<Bit> ConvoCode::decode(Codeword<Bit> a){
		Codeword<Bit> output;
		stack<int> s = a.getMessageStack();
		
		while (a.size() % n != 0) {
			a.push_back(0);
		}
		vector< vector<unsigned int> > preState, cost, preIn;
		
		preState.push_back( vector<unsigned int>(1<<m, 0) );
		
		cost.push_back( vector<unsigned int>(1<<m, 999999) );
		preIn.push_back( vector<unsigned int>(1<<m, 0) );
		
		cost[0][0] = 0;
		
		int *o = new int[n];
		for(int i=0; i< a.size()-max_forney*n ; i+=n){
			preState.push_back( vector<unsigned int>(1<<m, 0) );
			cost.push_back( vector<unsigned int>(1<<m, 99999999) );
			preIn.push_back( vector<unsigned int>(1<<m, 0) );
			
			int p = i/n;
			for(int fromState=0; fromState< (1 << m); fromState++){
				for(int in=0; in < (1 << k); in++){
					intToArray(o, trellis[fromState][in].out, n);
					int tempCost = 0;
					for(int s=0; s<n; s++){
						if (a[i+s].isErased()) {
							tempCost += 1;
						}
						else if(a[i+s].getValue() != o[s]){
							tempCost += 9999;
						}
					}
					
					if (cost[p+1][ trellis[fromState][in].next ] > cost[p][fromState] + tempCost) {
						cost[p+1][ trellis[fromState][in].next ] = cost[p][fromState] + tempCost;
						preState[p+1][ trellis[fromState][in].next ] = fromState;
						preIn[p+1][ trellis[fromState][in].next ] = in;
					}
				}
			}
			
			unsigned int t_min = -1;
			for (int j=0; j< (1<<m); j++) {
				if(cost[p+1][j] < t_min){
					t_min = cost[p+1][j];
				}
			}
			for (int j=0; j< (1<<m); j++) {
				cost[p+1][j] -= t_min;
			}
			
		}
		
		// ending stage
		
		for(int i=a.size()-max_forney*n; i< a.size(); i+=n){
			preState.push_back( vector<unsigned int>(1<<m, 0) );
			cost.push_back( vector<unsigned int>(1<<m, 99999999) );
			preIn.push_back( vector<unsigned int>(1<<m, 0) );
			
			int p = i/n;
			for(int fromState=0; fromState< (1 << m); fromState++){
				int in=0; 
				{
					intToArray(o, trellis[fromState][in].out, n);
					int tempCost = 0;
					for(int s=0; s<n; s++){
						if (a[i+s].isErased()) {
							tempCost += 1;
						}
						else if(a[i+s].getValue() != o[s]){
							tempCost += 9999;
						}
					}
					
					if (cost[p+1][ trellis[fromState][in].next ] > cost[p][fromState] + tempCost) {
						cost[p+1][ trellis[fromState][in].next ] = cost[p][fromState] + tempCost;
						preState[p+1][ trellis[fromState][in].next ] = fromState;
						preIn[p+1][ trellis[fromState][in].next ] = in;
					}
				}
			}
			
			unsigned int t_min = -1;
			for (int j=0; j< (1<<m); j++) {
				if(cost[p+1][j] < t_min){
					t_min = cost[p+1][j];
				}
			}
			for (int j=0; j< (1<<m); j++) {
				cost[p+1][j] -= t_min;
			}
			
		}
		
		
		// go through trellis
		int state = 0;
		list<unsigned int> out_temp;
		for (int i = preState.size() - 1; i>0; i--) {
			out_temp.push_front(preIn[i][state]);
			state = preState[i][state];
		}
		
		for (; !out_temp.empty(); out_temp.pop_front()) {
			intToArray(o, out_temp.front(), k);
			output.insert(output.end(), o, o+k);
		}
		
		
		delete [] o;
		output.trim(s.top());
		s.pop();
		output.setMessageStack(s);
		return output;
	}
}