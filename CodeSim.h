/*
 *  CodeSim.h
 *  new_LT
 *
 *  Created by 刁培倫 on 2010/8/25.
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
 */

#include<vector>
#include<string>
#include<stack>
using namespace std;

namespace CodeSim {
	
	template<class T>
	class Symbol{
	public:
		Symbol();
		
		bool isErased();		
		void setErased(bool t);
	protected:
		T value;
	private:
		bool erased;
	};
	
	
	
	class Bit : public Symbol<bool>{
	public:
		Bit();
		
		
		Bit(int t);
		

		Bit operator+(Bit t);		
		Bit operator*(Bit t);		
		string toString();
	};
	
	template<class S>
	class Codeword : public vector<S> {
	public:
		stack<unsigned int>  getMessageSize(){
			return messageSize;
		}
		
		void setMessageSize(stack<unsigned int> s){
			messageSize = s;
		}
	private:
		stack<unsigned int> messageSize;
	};
	
	template<class S1, class S2>
	class CodingBlock {
	public:
		virtual Codeword<S2> forward(Codeword<S1> c) = 0;
		virtual Codeword<S1> backword(Codeword<S2> c) = 0;
		
	private:
	};
	
	template<class T>
	class Permutator : public CodingBlock<T,T> {
	public:
		Permutator(){
			blockSize = 1;
			permutationTable.assign(1,1);
			depermutationTable.assign(1,1);
		}
		Permutator(string filename, bool inverseOrder);
		Permutator(string filename){
			Permutator(filename, false);
		}
		Codeword<T> permutate(Codeword<T>);
		Codeword<T> depermutate(Codeword<T>);
		Codeword<T> forward(Codeword<T> c){
			return permutate(c);
		}
		
		Codeword<T> backward(Codeword<T> c){
			return depermutate(c);
		}
		
	private:
		unsigned int blockSize;
		vector<unsigned int> permutationTable, depermutationTable;
	};
}

