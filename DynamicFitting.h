/*
 *  RegressionLib.h
 *  new_LT
 *
 *  Created by 刁培倫 on 2010/12/16.
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
 */

#include "nr.h"
#include "randomc.h"
#include <vector>
#include <utility>
#include <cmath>
#include <ctime>

#ifndef DYNAMICFITTING_H
#define DYNAMICFITTING_H

void sigmoid(const DP x, Vec_I_DP &a, DP &y, Vec_O_DP &dyda);
void sigmoidRange(vector<double> &x,
				  vector<double> &y,
				  vector< pair<double,double> >& a);
int regression(vector<double> &data_x, vector<double> &data_y, vector<double> &guess_a
			   ,void funcs(const DP, Vec_I_DP &, DP &, Vec_O_DP &) 
			   ,double &ss);
void dynamicFitting(vector<double> &data_x, 
					vector<double> &data_y, 
					void funcs(const DP, Vec_I_DP &, DP &, Vec_O_DP &),
					void argRange(vector<double> &,
								  vector<double> &,
								  vector< pair<double,double> >&),
					vector<double> &a,
					double &ss);

void x25_50_75(vector<double> &data_x, vector<double> &data_y, vector<double> & x25_50_75);

double find_x(vector<double> &data_x, vector<double> &data_y, double target);

double min(vector<double> &data);
double max(vector<double> &data);




void sigmoid(const DP x, Vec_I_DP &a, DP &y, Vec_O_DP &dyda)
{
	DP ex;
	
	if(a.size()!=5)
		cerr << "sigmoid: wrong # of parameters\n";
	ex = pow(1+exp((-x+a[3])/a[1]),(-a[2]));
	y=a[4]+a[0]*ex;
	dyda[0] = ex;
	dyda[2] = -a[0]*log(1+exp((-x+a[3])/a[1]))*ex;
	ex /= 1+exp((-x+a[3])/a[1]);
	dyda[1] = a[0]*a[2]*(-x+a[3])*exp((-x+a[3])/a[1])/a[1]/a[1]*ex;
	dyda[3] = -a[0]*a[2]*exp((-x+a[3])/a[1])/a[1]*ex;
	dyda[4] = 1;
}

void sigmoidRange(vector<double> &x,
				  vector<double> &y,
				  vector< pair<double,double> >& a){


	a.clear();
	a.reserve(5);
	pair<double,double>  t;
	double y_min = min(y), y_max = max(y);
	
//	a = max(y)-min(y) ''Auto {{previous: 3.81787}}	
	t.first = y_max-y_min;
	t.second = -t.first;
	t.first *= 3;
	a.push_back(t);
	
//	b = xwtr(x,y-min(y),.5)/4 ''Auto {{previous: -0.0300599}}
	t.first = find_x(x, y, y_min + (y_max-y_min)*0.75) 
						- find_x(x, y, y_min + (y_max-y_min)*0.25);
	t.first /= 4;
	t.second = -t.first;
	t.first *= 3;
	a.push_back(t);
//	c = 1 ' {{previous: 1.13287}}
	t.first = 3;
	t.second = 0;
	a.push_back(t);
//	x0 = x50(x,y-min(y),.5) ''Auto {{previous: 0.136189}}
	t.first = find_x(x, y, (y_max + y_min)/2 );
	
	t.second = -t.first;
	t.first *= 3;
	a.push_back(t);
//	y0 = min(y) ''Auto {{previous: -3.98481}}	
	t.first = y_min;
	t.second = -t.first;
	t.first *= 3;
	a.push_back(t);
}

int regression(vector<double> &data_x, vector<double> &data_y, vector<double> &guess_a
			   ,void funcs(const DP, Vec_I_DP &, DP &, Vec_O_DP &) 
			   ,double &ss){
	if (data_x.size() != data_y.size()) {
		cerr << "regression: data size not matched."<<endl;
		return -1;
	}
	if (data_x.size() == 0) {
		cerr << "regression: data size is zero."<<endl;
		return -1;
	}
	int NPT = data_x.size(), MA = guess_a.size();
	Vec_DP x(NPT),y(NPT),sig(NPT);
	Vec_BOOL ia(MA);
	Vec_DP a(MA);
	Mat_DP covar(MA,MA),alpha(MA,MA);
	DP alamda,chisq,ochisq;
	
	for (int i=0;i<NPT;i++) {
		x[i]=data_x[i];
		y[i]=data_y[i];
		//		  for (j=0;j<MA;j+=3) {
		//			y[i] += a[j]*exp(-SQR((x[i]-a[j+1])/a[j+2]));
		//		  }
		//		  y[i] *= (1.0+SPREAD*NR::gasdev(idum));
		sig[i]=1;//SPREAD;
	}
	for (int i=0;i<MA;i++) ia[i]=true;
	for (int i=0;i<MA;i++) a[i]=guess_a[i];
	
	alamda = -1;
	NR::mrqmin(x,y,sig,a,ia,covar,alpha,chisq,sigmoid,alamda);
	int k=0;
	int itst=0;
	while (itst < 4 && k<200) {
		//		cout << endl << "Iteration #" << setw(3) << k;
		//		cout << setw(18) << "chi-squared:" << setw(13) << chisq;
		//		cout << setw(11) << "alamda:" << setw(10) << alamda << endl;
		//		cout << setw(8) << "a[0]" << setw(9) << "a[1]";
		//		cout << setw(9) << "a[2]" << setw(9) << "a[3]";
		//		cout << setw(9) << "a[4]" << setw(9) << "a[5]" << endl;
		//		cout << fixed << setprecision(6);
		//		for (int i=0;i<MA;i++) cout << ' ' << a[i];
		//		cout << endl;
		k++;
		ochisq=chisq;
		NR::mrqmin(x,y,sig,a,ia,covar,alpha,chisq,sigmoid,alamda);
		fabs(ochisq-chisq) < 0.0000000001 ? itst++ : itst=0;
		//if (itst > 4 ||  ) continue;
		
	}
	
	ss = chisq;
	for (int i=0;i<MA;i++) guess_a[i] = a[i];
	alamda=0.0;
	//	NR::mrqmin(x,y,sig,a,ia,covar,alpha,chisq,sigmoid,alamda);
	//	cout << endl << "Uncertainties:" << endl;
	//	for (int i=0;i<MA;i++) cout << setw(9) << sqrt(covar[i][i]);
	//	cout << endl;
	//	cout << endl << "Expected results:" << endl;
	//	cout << setw(9) << 5.0 << setw(9) << 2.0 << setw(9) << 3.0;
	//	cout << setw(9) << 2.0 << setw(9) << 5.0 << setw(9) << 3.0 << endl;
	return k;
}

double min(vector<double> &data){
	if (data.size() == 0) {
		return 0;
	}
	
	double min = data[0];
	for (vector<double>::iterator i = data.begin()+1; i!=data.end(); i++) {
		if (*i< min) {
			min = *i;
		}
	}
	
	return min;
}
double max(vector<double> &data){
	if (data.size() == 0) {
		return 0;
	}
	
	double max = data[0];
	for (vector<double>::iterator i = data.begin()+1; i!=data.end(); i++) {
		if (*i> max) {
			max = *i;
		}
	}
	
	return max;
}

double find_x(vector<double> &data_x, vector<double> &data_y, double target){
	if (data_x.size() != data_y.size()) {
		cerr << "find_x: data size not matched." << endl;
	}
	if (data_x.size() == 0) {
		cerr << "find_x: empty data." << endl;
	}
	
	double diff = abs(target - data_y[0]);
	double x = data_x[0];
	
	for (int i=1; i<data_x.size(); i++) {
		if (abs(target - data_y[i]) < diff) {
			diff = abs(target - data_y[i]);
			x = data_x[i];
		}
	}
	
	return x;
}

void dynamicFitting(vector<double> &data_x, 
					vector<double> &data_y, 
					void funcs(const DP, Vec_I_DP &, DP &, Vec_O_DP &),
					void argRange(vector<double> &,
								  vector<double> &,
								  vector< pair<double,double> >&),
					vector<double> &a,
					double &ss){
	
	vector< pair<double,double> > aRange;
	
	argRange(data_x,data_y,aRange);
	
	CRandomMersenne ran(time(0));
	
	vector<double> ran_a(aRange.size(), 0);
	vector< vector<double> > a_list;
	for (int i=0; i<aRange.size(); i++) {
		a_list.push_back(vector<double>());
		a_list[i].assign(500, 0);
		for (int j=0; j<500; j++) {
			a_list[i][j] = aRange[i].first + j*(aRange[i].second - aRange[i].first)/500.0;
		}
		
		for (int j=0; j<500; j++){
			int b = ran.IRandom(j,500);
			swap(a_list[i][j], a_list[i][b]);
		}
		
	}
	
	
	double min_ss = (1 << 10);
	for (int fit=0; fit<500; fit++) {
		for (int i=0; i<aRange.size(); i++) {
			ran_a[i] = a_list[i][fit];
			cout << "\t" << ran_a[i];
		}
		cout << endl;
		double ss_t;
		regression(data_x, data_y, ran_a, funcs,ss_t);
		if (ss_t < min_ss) {
			min_ss == ss_t;
			a = ran_a;
			ss = ss_t;
		}
	}
	
	
}

#endif
