/*
 *  ConvoCode.cpp
 *  new_LT
 *
 *  Created by 刁培倫 on 2010/9/7.
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
 */

#include "ConvoCode.h"

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
		
		
		output.trim(s.top());
		s.pop();
		output.setMessageStack(s);
		return output;
	}
}