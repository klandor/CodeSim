import sys
from scipy import array
import math
import pybrain
from pybrain.optimization import ExactNES

# if(len(sys.argv) == 1):
# 	print "Usage: NES_fit.py filename"
# 	sys.exit()

def parseInt(a):
	b=[]
	for i in range(len(a)):
		b.append(int(a[i]))
	return b

def parseFloat(a):
	b=[]
	for i in range(len(a)):
		b.append( float(a[i]) )
	return b

def toStringList(a):
	b=[]
	for i in range(len(a)):
		b.append( `a[i]` )
	return b


def normalize(p):
	for i in range(len(p)):
		if p[i] < 0:
			p[i] = -p[i]
	from scipy import sum
	sum_p=float(sum(p))
	for i in range(len(p)):
		p[i] /= sum_p

def mygetline(file):
	while(True):
		s=file.readline()
		if(s == ''):
			return s
		if(s[0] != '#'):
			return s

def log10(x):
	return math.log(x, 10)

def abs(x):
	return math.fabs(x)


class LT_exp:
	def __init__(self, filename):
		f=open(filename,'r')
		self.K = int(mygetline(f))
		self.Dsize = int(mygetline(f))
		self.Tags = parseInt(mygetline(f).split())
		if(self.Dsize != len(self.Tags)):
			print "Error: Dsize and Tags doesn't match"
			sys.exit()
		self.D = parseFloat(mygetline(f).split())
		if(self.Dsize != len(self.D)):
			print "Error: Dsize and initial distribution doesn't match"
			sys.exit()
		tmp=mygetline(f).split()
		self.STEPS = int(tmp[0])
		self.Delta = float(tmp[1])
		mygetline(f) #cubic fitting
		mygetline(f) #error exponent
		self.epsilons = parseFloat(mygetline(f).split())
		self.epsilons_w = parseFloat(mygetline(f).split())
		self.targetErrorRate = parseFloat(mygetline(f).split())
		self.targetErrorRate_w = parseFloat(mygetline(f).split())


	def fitness_LT(self, dist):
		LT_BER_format = "{K} {Run}\n{Dsize}\n{tag_list}\n{distribution}\n{STEPS}\n{epsilon_list}"
		STEPS = 61
		Delta = 0.005
		Run = 10000
		MaxEpsilon = (Delta*(STEPS-1))
		P_e_min = log10(1/float(self.K*Run))
		
		epsilon_list = array(parseFloat(range(STEPS)))*Delta
		str_epsilon_list = []
		for i in range(STEPS):
			str_epsilon_list.append('%.3f' % epsilon_list[i])
		
		d=dist[0:]
		normalize(d)
		from subprocess import Popen
		from subprocess import PIPE
		
		p = Popen('./LT_BER.out', stdin=PIPE, stdout=PIPE)
		input = LT_BER_format.format(K=self.K, Run=Run, Dsize=self.Dsize, tag_list=' '.join(toStringList(self.Tags)), distribution=' '.join(toStringList(d)), STEPS=STEPS, epsilon_list= ' '.join(str_epsilon_list))
		#K='`self.K`', Run='100', Dsize='`self.Dsize`', tag_list=' '.join(self.Tags), distribution='`self.D`', STEPS='16', epsilon_list= "")
		#print input
		err =  p.communicate(input)[0]
		err = err.split()
		err = parseFloat(err)
		for i in range(len(err)):
			if(err[i]>0):
				err[i] = log10(err[i])
			else:
				err[i] = P_e_min
		
		
		fit = 0
		
		for i in range(len(self.epsilons)) :
			fit += (err[int(self.epsilons[i]/Delta)] - P_e_min)  / abs( P_e_min ) * self.epsilons_w[i];
			#parameters[i] = (err[int(self.epsilons[i]/Delta)] - P_e_min)  / abs( P_e_min );
		
		for i in range(len(self.targetErrorRate)):	
			min_diff=999
			for j in range(STEPS):
				if ( abs(err[j]-log10(self.targetErrorRate[i])) < min_diff):
					min_diff = abs(err[j]-log10(self.targetErrorRate[i]))
					min_i = j
			fit += min_i/float(STEPS-1) * self.targetErrorRate_w[i];

		
		##print err
		return fit

lt = LT_exp(sys.argv[1])
# d = lt.D
# d = d.split()
# d =parseFloat(d)
#print lt.fitness_LT(lt.D)
nes = ExactNES(lt.fitness_LT, lt.D, minimize=True, maxEvaluations=10000, verbose=True)
print nes.learn()
# p = pow(2).power
# from pybrain.optimization import ExactNES
# print ExactNES(pow(32.60).power, [0], maxEvaluations = 10000, minimize=True).learn()[0][0]
