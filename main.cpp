#include <iostream>
#include <limits>
#include "LT.h"
using namespace std;
using namespace CodeSim;
int main (int argc, char * const argv[]) {
    
	//LT_sim sim;
	CRandomMersenne r(5565);
	
	
	int Degree_of_Edge[10] = {1, 2, 3, 4, 5, 7, 9, 19, 59, 179};
	double Omega[10] = {7.9379E-02, 4.0129E-01, 1.0121E-01, 2.1679E-01, 5.0996E-02, 
				5.8338E-05, 3.9740E-02, 7.7470E-02, 2.1520E-02, 1.1547E-02};
	double sum[16];
	for (int i = 0; i< 16; i++)
		sum[i] = 0;
	
	int K = 1000;
	
	for (int i=0; i<1000; i++) {
		cout << "Run: " << i << '\n';
		LT_sim<int> sim;
		sim.init(K, K*1.2, 10, Degree_of_Edge, Omega, r.IRandom(0, INT_MAX));
		for (int i = 0; i< 16; i++) {
//			double t = Encoder(K, K*(1.05+0.01*i), Dsize);
//			ErrorCount[i][(int)(t*16)]++;
//			BER[i]+=t;
		
			sim.seqReceive(K*(1.05+0.01*i) );
			sim.decode();
			sum[i] += sim.run();
		}
	}
	for (int i = 0; i< 16; i++)
		cout << sum[i] / 1000 << '\t';
	
	
    return 0;
}
