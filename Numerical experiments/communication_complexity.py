import numpy as np
import matplotlib.pyplot as plt
import math

n = 7
t = 2
q = 2**61-1
logq = math.ceil(math.log2(q))
lmd = 128
logNmin = 4
logNmax = 9
N = [2**i for i in range(logNmin, logNmin + logNmax)]

# semihonest
comm_semihonest = np.zeros(logNmax)
for i in range(logNmax): # i = 0 ... logNmax-1
    comm_semihonest[i] = n**2 * logq * N[i]

# baseline
comm_baseline = np.zeros(logNmax)
for i in range(logNmax): # i = 0 ... logNmax-1
    k = n - t
    comm_baseline[i] = lmd * k**2 * math.comb(n,t) + n**2 * logq * N[i]

# allsets
comm_allsets = np.zeros(logNmax)
for i in range(logNmax): # i = 0 ... logNmax-1
    k = t + 1
    comm_allsets[i] = lmd * k**2 * math.comb(n,t+1) + n**2 * logq * N[i]

# partition
comm_partition = np.zeros(logNmax)
for i in range(logNmax): # i = 0 ... logNmax-1
    k = t + 1
    ell = int(math.lcm(n,t+1)/(t+1))
    comm_partition[i] = lmd * k**2 * ell + n**2 * logq * N[i]

# binomial
comm_binomial = np.zeros(logNmax)
for i in range(logNmax): # i = 0 ... logNmax-1
    k = n - t
    comm_binomial[i] = lmd * k**2 * math.comb(n,t) + n**2 * logq * N[i]
