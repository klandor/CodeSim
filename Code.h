/*
 *  Code.h
 *  new_LT
 *
 *  Created by 刁培倫 on 2010/8/13.
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
 */

#include <vector>
using namespace std;

namespace CodeSim{
	
	template<class T>
	class Code {
		virtual vector<T> encode(vector<T> a) = 0;
		virtual vector<T> decode(vector<T> a, vector<bool> erasure) = 0;
		virtual void reset() = 0;
	};
}