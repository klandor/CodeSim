import sys
from scipy import array
import math
import pybrain
from pybrain.optimization import ExactNES
import os
os.nice(20)
if(len(sys.argv) == 1):
	print "Usage: NES_fit.py filename"
	sys.exit()

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
		#read from config file
		f=open(filename,'r')
		tmp=mygetline(f).split()
		self.K = int(tmp[0])
		self.Run = int(tmp[1])
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
		
		# prepare LT_BER process
		LT_BER_format = "{K} {Run}\n{Dsize}\n{tag_list}\n{STEPS}\n{epsilon_list}\n"
		self.STEPS = 61
		self.Delta = 0.005
		self.MaxEpsilon = (self.Delta*(self.STEPS-1))
		self.P_e_min = log10(1/float(self.K*self.Run))
		
		self.epsilon_list = array(parseFloat(range(self.STEPS)))*self.Delta
		str_epsilon_list = []
		for i in range(self.STEPS):
			str_epsilon_list.append('%.3f' % self.epsilon_list[i])
		from subprocess import Popen
		from subprocess import PIPE
		
		self.p = Popen('../../LT_BER.out', stdin=PIPE, stdout=PIPE)
		input = LT_BER_format.format(K=self.K, Run=self.Run, Dsize=self.Dsize, tag_list=' '.join(toStringList(self.Tags)), STEPS=self.STEPS, epsilon_list= ' '.join(str_epsilon_list))
		self.p.stdin.write(input)
		
		#prepare listener for ExactNES to write out results
		self.file_result = open('result_%s' % filename, 'w')
		self.result_line = 0
		self.times_of_eval = 0
	def fitness_LT(self, dist):
		self.times_of_eval += 1
		
		d=dist[0:]
		normalize(d)
		input=' '.join(toStringList(d))
		#K='`self.K`', Run='100', Dsize='`self.Dsize`', tag_list=' '.join(self.Tags), distribution='`self.D`', STEPS='16', epsilon_list= "")
		#print input
		self.p.stdin.write('%s\n' % input)
		err = self.p.stdout.readline()
		err = err.split()
		err = parseFloat(err)
		for i in range(len(err)):
			if(err[i]>0):
				err[i] = log10(err[i])
			else:
				err[i] = self.P_e_min
		
		
		fit = 0
		
		for i in range(len(self.epsilons)) :
			fit += (err[int(self.epsilons[i]/self.Delta)] - self.P_e_min)  / abs( self.P_e_min ) * self.epsilons_w[i];
			#parameters[i] = (err[int(self.epsilons[i]/Delta)] - P_e_min)  / abs( P_e_min );
		
		for i in range(len(self.targetErrorRate)):	
			min_diff=999
			for j in range(self.STEPS):
				if ( abs(err[j]-log10(self.targetErrorRate[i])) < min_diff):
					min_diff = abs(err[j]-log10(self.targetErrorRate[i]))
					min_i = j
			fit += min_i/float(self.STEPS-1) * self.targetErrorRate_w[i];
		
		#print fit
		
		return fit
	def printResult(self, parameter_list, eval):
		self.result_line += 1
		p_list = parameter_list[0:]
		normalize(p_list)
		self.file_result.write(`self.result_line`+'\t'+ `self.times_of_eval`+'\teval\t'+`eval`+'\tpara\t' + ' '.join(toStringList(p_list))+'\n')
		self.file_result.flush()


# class ResultWriter:
# 	def __init__(self, filename):
# 		self.filename = filename
# 		self.file = open(filename, 'w')
# 	def printLine(self, line):
# 		self.file.write(line+'\n')
# 		self.file.flush()
# 	def close(self):
# 		self.file.close()

"""main()"""
filename = sys.argv[1]
lt = LT_exp(filename)
# out = dup(sys.stdout, open('result_%s' % filename, 'w'))
# sys.stdout = out
# writer = ResultWriter('dist_%s' % filename)
nes = ExactNES(lt.fitness_LT, lt.D, minimize=True, maxEvaluations=10000, verbose=True, listener=lt.printResult)
nes_result = nes.learn()
final = nes_result[0]
normalize(final)
final = toStringList(final)
print ' '.join(final)
# dist = open('dist_%s' % filename, 'w')
# dist.write(' '.join(final) + '\n')
# dist.close()
# p = pow(2).power
# from pybrain.optimization import ExactNES
# print ExactNES(pow(32.60).power, [0], maxEvaluations = 10000, minimize=True).learn()[0][0]
