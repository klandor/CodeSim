/*
 *  fitting_data_set_generator.cpp
 *  new_LT
 *
 *  Created by 刁培倫 on 2010/11/14.
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
 */


#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <algorithm>
using namespace std;



int main(){
	string s;
	vector< vector<double> > data;

	getline(cin, s);
	{
		istringstream iss(s);
		double d;
		int i=0;
		while (iss >> d) {
			data.push_back(vector<double>());
			data[i].push_back(d);
			i++;
		}
	}
	
	while ( getline(cin, s) ) {
		istringstream iss(s);
		double d;
		for (int i=0; iss >> d; i++) {
			data[i].push_back(d);
		}
	}
	
	for (vector< vector<double> >::iterator i=data.begin(); i!=data.end(); i++) {
		
		sort(i->begin(), i->end());
		
		int size = i->size();
		for (int s=0; s<100; s++) {
			double sum=0;
			int n=0;
			for (int x=s*(size/100); x< (s+1)*(size/100); x++) {
				sum += (*i)[x];
				n++;
			}
			cout << sum/n << '\t';
		}
		cout << '\n';
	}
	return 0;
}