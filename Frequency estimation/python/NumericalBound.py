#!/usr/bin/env python3
#import computeamplification as CA
import computeamplification_rev as CA
import math
import sys
import csv

################################# Parameters ##################################
#sys.argv = ["NumericalBound.py", "107614", "."]
#sys.argv = ["NumericalBound.py", "896308", "."]
if len(sys.argv) < 3:
    print("Usage:",sys.argv[0],"[n] [ResDir]" )
    sys.exit(0)

#number of local reports
n = int(sys.argv[1])
#result directory
ResDir = sys.argv[2]

#desired delta of final shuffled privacy guarantee
delta = 10**(-8)
#number of iterations of binary search. The higher T is, the more accurate the result
num_iterations = 10
#This is a parameter of the empirical analysis computation that can be tuned for efficiency. The larger step is, the less accurate the result, but more efficient the algorithm.
step = 100

#local epsilon (max)
epsl_max = math.log(n / (16 * math.log(2 / delta)))
#local epsilon (min)
#epsl_min = 0.001
epsl_min = 0.01
#local epsilon (shift)
#epsl_shift = 0.001
epsl_shift = 0.01

outfile = ResDir + "/numerical-bound_n" + str(n) + "_d10-8.csv"
f = open(outfile, "w")

epsl = epsl_max
print("epsL,eps_lower,eps_upper", file=f)
writer = csv.writer(f, lineterminator="\n")
while epsl > epsl_min:
    #There are 2 main functions, empiricalanalysis and theoryanalysis.
    #empiricalanalysis computes the shuffled privacy guarantee empirically. The 1 and 0 correspond to returning either an upper bound on the privacy guarantee, or a lower bound.
    numerical_upperbound = CA.numericalanalysis(n, epsl, delta, num_iterations, step, True)
    numerical_lowerbound = CA.numericalanalysis(n, epsl, delta, num_iterations, step, False)
    print(epsl, numerical_lowerbound, numerical_upperbound)
    lst = [epsl, numerical_lowerbound, numerical_upperbound]
    writer.writerow(lst)
    epsl = epsl - epsl_shift
f.close()
