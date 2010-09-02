/*
 *  Symbol.cpp
 *  new_LT
 *
 *  Created by 刁培倫 on 2010/8/27.
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
 */


#include<vector>
#include<string>
#include"CodeSim.h"
using namespace std;

namespace CodeSim {

	
	
	
	Bit::Bit(){
			value = 0;
	}
		
	//template<class T>
	Bit::Bit(int t){
		value = t;
	}
		
//	template<class T>
//	Bit Bit::operator=(T t){
//		value = t;
//		return *this;
//	}
	Bit Bit::operator+(Bit t){
		return Bit(this->value ^ t.value);
	}
	
	Bit Bit::operator*(Bit t){
		return Bit(this->value & t.value);
	}
	
	string Bit::toString()
	{
		if (isErased()) {
			return "-1";
		}
		if (value) {
			return "1";
		}
		else {
			return "0";
		}
	}
		
	
}