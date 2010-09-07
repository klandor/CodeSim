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
	
	
	Codeword<Bit> a;
	char c;
	while (cin >> c) {
		if (c == '0') {
			a.push_back(0);
		}
		else if (c == '1'){
			a.push_back(1);
		}
		else {
			break;
		}

	}
	
	Codeword<Bit> b = cc.encode(a);
	
	for (int i = 0; i<b.size(); i++) {
		cout << b[i].toString();
	}
	cout << '\n';
	
	
	
	return 0;
}