"""
LT_BER input
1000 1000
20
1	2	3	4	5	7	8	11	16	19	23	32	41	53	64	83	101	128	163	199	
0.0928815	    0.368117	    0.188153	   0.0521204	    0.129776	  0.00526887	    0.027561	   0.0194776	   0.0108784	  0.00835065	   0.0213981	   0.0363515	   0.0012552	  0.00026333	 0.000510834	  0.00186076	 0.000139882	  0.00165875	 0.000187628	   0.0337885	
16
0	0.03	0.06	0.09	0.12	0.15	0.18	0.21	0.24	0.27	0.3	0.33	0.36	0.39	0.42	0.45
"""

def fitness(parameters):
	p=parameters[0:];
	for i in range(len(p)):
		if p[i] < 0:
			p[i] = -p[i]
	sum=sum(p)
	for i in range(len(p)):
		p[i] /= sum
	return p
a = [1, 0, 3]
print fitness(a)

