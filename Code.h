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
	
	class Code {
		virtual vector<char> encode(vector<char> a) = 0;
		virtual vector<char> decode(vector<char> a) = 0;
		virtual void reset() = 0;
	};
}