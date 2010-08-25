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
using namespace std;

namespace CodeSim {
	
	template<class T>
	class Symbol{
	public:
		bool inErased(){
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
		Bit(){
			value = 0;
		}
		
		template<class T>
		Bit(T t){
			value = t;
			
		}
		
		template<class T>
		Bit operator=(T t){
			value = t;
			return *this;
		}
		Bit operator+(Bit t){
			return Bit(this->value ^ t.value);
		}
		
		Bit operator*(Bit t){
			return Bit(this->value & t.value);
		}
		
		string toString()
		{
			if (value) {
				return "1";
			}
			else {
				return "0";
			}
		}

		
	};
	
	template<class S1, class S2>
	class CodingBlock {
	public:
		virtual S2 forward(S1 c) = 0;
		virtual S1 backword(S2 c) = 0;
		
	private:
	};
}

