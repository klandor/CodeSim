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
#include<fstream>
#include<iostream>
using namespace std;

namespace CodeSim {
	
	template<class T>
	class Symbol{
	public:
		Symbol(){
			erased = false;
		}
		
		
		bool isErased(){
			return erased;
		}
		
		void setErased(bool t){
			erased = t;
		}
		
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
		Codeword(){};
		Codeword(int size, S s){
			vector<S>(size,s);
		}
		stack<int>  getMessageSize(){
			return messageSize;
		}
		
		void setMessageSize(stack<int> s){
			messageSize = s;
		}
		
		void trim(int m){
			if(m < this->size())
			{
				this->erase(this->begin()+m, this->end());
			}
		}
	private:
		stack<int> messageSize;
	};
	
	template<class S1, class S2>
	class CodingBlock {
	public:
		
//		virtual Codeword<S2> forward(Codeword<S1> c) = 0;
//		virtual Codeword<S1> backward(Codeword<S2> c) = 0;
		
	private:
	};
	
	template<class T>
	class Permutator : public CodingBlock<T,T> {
	public:
		
		// default construct, blockSize = 1 (no permutation).
		Permutator(){
			blockSize = 1;
			permutationTable.assign(1,1);
			depermutationTable.assign(1,1);
		}
		
		/* 
			read permutation table from file "filename"
			inverseOrder = true would inverse the table (permutate() become depermutate(), and vice versa)
		*/
		Permutator(const string filename, bool inverseOrder){
			ifstream in(filename.c_str());
			if(in.fail())
			{
				cerr << "Permutator initial error: \""<< filename << "\" can't be opened." ;
				Permutator();
				return;
			}
			
			
			depermutationTable.clear();
			
			
			permutationTable.clear();
			unsigned int t;
			while (in >> t) {
				permutationTable.push_back(t);
			}
			
			blockSize = permutationTable.size();
			depermutationTable.assign(blockSize, 0);
			
			for (int i = 0; i< permutationTable.size(); i++) {
				t = permutationTable[i];
				depermutationTable[t] = i;
			}
			
			
			if(inverseOrder){
				permutationTable.swap(depermutationTable);
			}
			
			
		}
		
		/** 
		 * @brief read permutation table from file "filename"
		 * @param
		 */
		
		Permutator(const string filename){
			Permutator(filename, false);
		}
		
		
		Codeword<T> permutate(Codeword<T> c){
			Codeword<T> output;
			for (int i = 0; i< c.size(); i+=blockSize) {
				Codeword<T> temp;
				temp.assign(blockSize, 0);
				for (int j=0; j< blockSize && i+j < c.size(); j++) {
					temp[permutationTable[j]] = c[i+j];
				}
				
				output.insert(output.end(), temp.begin(), temp.end());
			}
			stack<int> s = c.getMessageSize();
			s.push(c.size());
			output.setMessageSize(s);
			return output;
		}
		
		Codeword<T> depermutate(Codeword<T> c){
			Codeword<T> output;
			for (int i = 0; i< c.size(); i+=blockSize) {
				Codeword<T> temp;
				temp.assign(blockSize, 0);
				for (int j=0; j< blockSize && i+j < c.size(); j++) {
					temp[depermutationTable[j]] = c[i+j];
				}
				
				output.insert(output.end(), temp.begin(), temp.end());
			}
			stack<int> s = c.getMessageSize();
			output.trim(s.top());
			s.pop();
			output.setMessageSize(s);
			return output;
		}
		
		Codeword<T> forward(Codeword<T> c){
			return permutate(c);
		}
		
		Codeword<T> backward(Codeword<T> c){
			return depermutate(c);
		}
		
		
		void print(){
			for(int i = 0; i< blockSize; i++){
				cout << depermutationTable[i] << endl;
			}
		}
	private:
		int blockSize;
		vector<unsigned int> permutationTable, depermutationTable;
	};
}

