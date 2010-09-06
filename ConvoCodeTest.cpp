/*
 *  ConvoCode.cpp
 *  new_LT
 *
 *  Created by 刁培倫 on 2010/9/6.
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
 */


#include "ConvoCode.h"
#include "LTCode.h"

using namespace std;
using namespace CodeSim;
int main(){
	cout << octToDec(19) << "\n";
	ConvoCode cc("convo.txt");
	cc.showInfo();
	
	return 0;
}