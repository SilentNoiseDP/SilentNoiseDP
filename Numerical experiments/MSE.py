import numpy as np
import matplotlib.pyplot as plt
import math

n = 3
t = 1
D = 1
eps = 0.1
delta = 10**(-8)
p = 1 - np.exp(-eps/D) # 1 - p = alpha**(-1)
N = 2**10 # repetition

# semihonest
noise_semihonest = np.zeros(N)
for h in range(N):
    m = n - t
    for j in range(n):
        noise_semihonest[h] = noise_semihonest[h] + (np.random.negative_binomial(1/m, p) - np.random.negative_binomial(1/m, p))
var_semihonest = np.var(noise_semihonest)

# baseline
noise_baseline = np.zeros(N)
for h in range(N):
    for j in range(math.comb(n, t)):
        noise_baseline[h] = noise_baseline[h] + (np.random.geometric(p) - np.random.geometric(p))
var_baseline = np.var(noise_baseline)

# allsets
noise_allsets = np.zeros(N)
for h in range(N):
    m = math.comb(n-t, t+1)
    for j in range(math.comb(n, t+1)):
        noise_allsets[h] = noise_allsets[h] + (np.random.negative_binomial(1/m, p) - np.random.negative_binomial(1/m, p))
var_allsets = np.var(noise_allsets)

# partition
noise_partition = np.zeros(N)
for h in range(N):
    ell = int(math.lcm(n,t+1)/(t+1))
    m = math.ceil(ell*(1-t*(t+1)/n))
    for j in range(ell):
        noise_partition[h] = noise_partition[h] + (np.random.negative_binomial(1/m, p) - np.random.negative_binomial(1/m, p))
var_partition = np.var(noise_partition)

# binomial
d = 1
D1 = 1 
D2 = np.sqrt(1)
Di = 1
c1 = 2*D2*np.sqrt(2*np.log(1.25/delta))
c2 = 4/(1-delta/10) * (D2*7*np.sqrt(2)/4 * np.sqrt(np.log(10/delta)) + D1*1/3) \
    + 4 * (Di*(2/3)*np.log(1.25/delta) + Di*(2/3)*np.log(20*d/delta)*np.log(10/delta))
K = 10**6
M =  math.floor(eps / (c1 / np.sqrt(K) + c2 / K))
noise_binomial = np.zeros(N)
for h in range(N):
    for j in range(math.comb(n,t)):
        noise_binomial[h] = noise_binomial[h] + (1/M)*(np.random.binomial(K, 1/2) - K/2)
var_binomial = np.var(noise_binomial)

