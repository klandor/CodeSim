/*
 *  ConvoCodeSim.cpp
 *  new_LT
 *
 *  Created by 刁培倫 on 2011/5/7.
 *  Copyright 2011 __MyCompanyName__. All rights reserved.
 *
 */

#include "ConvoCode.h"
#include "LTCode.h"
#include "randomc.h"
#include <ctime>
#include <cmath>
#define L 100000
#define Layer 4
#include <omp.h>

using namespace CodeSim;
using namespace std;

int main(int argn, char **args){
	
	
	CRandomMersenne r(time(0));
	ConvoCode cc("convo-2-4-6-6.txt");
	cout << "ChannelErasureRate";
	for (int i=0; i<Layer; i++) {
		cout << "\tLayer" << i+1;
	}
	cout << "\tTotalBits\n";
	
	double errRate = 0.1;
	while (1) {
		unsigned long err[Layer], total=0;
		
		for (int i=0; i<Layer; i++) {
			err[i]=0;
		}
		
		while (1) {
			#pragma omp parallel for num_threads(6)
			for (int i=0; i<6; i++) {
				Codeword<Bit> a;
				a.reserve(Layer*L);
				for (int i=0; i<Layer*L; i++) {
					a.push_back(r.IRandomX(0, 1));
				}
				Codeword<Bit> b = cc.encode(a);
				
				for (int i=0; i< b.size(); i++) {
					if (r.Random() < errRate) {
						b[i].setErased(true);
					}
				}
				
				Codeword<Bit> c = cc.decode(b);
				
				for (int i=0; i<c.size(); i++) {
					if (!(a[i] == c[i])) {
						#pragma omp atomic
						err[i%4]++;
					}
				}
				#pragma omp atomic
				total += L;
			}
			// check stop
			if (total >= 1000*L) {
				break;
			}
			
			bool enough = true;
			for (int i=0; i<Layer; i++) {
				if (err[i] < 100) {
					enough = false;
				}
			}
			
			if (enough) {
				break;
			}
		}
		
		// enough # of errors
		
		cout << errRate;
		
		for (int i=0; i<Layer; i++) {
			cout << '\t' << err[i]/(double)total;
		}
		
		cout << '\t' << total << endl;
		
		errRate *= pow(10, -1.0/5);
		
	}
	
}