import numpy as np
from sympy import *
import math
import itertools
import time

def GenM(n, k, Z):
    M = np.zeros((n, len(Z)))    
    for j in range(len(Z)):
        A = Z[j]
        for i in range(k):
            M[A[i]][j] = 1
    return M

def GenP(n, k, Z, q):
    s = n - k + 1
    P = np.zeros((s, len(Z)))
    for j in range(len(Z)):
        A = Z[j]
        a = np.zeros(s); i = 1
        for h in range(n):
            if h not in A:
                a[i] = h + 1
                i = i + 1
        vlist = []
        for i in range(s):
            row =[int(a[i]**j) for j in range(s)]
            vlist.append(row)
        V = Matrix(vlist)
        V_inv = matrix2numpy(V.inv_mod(q))
        for i in range(s):
            P[i][j] = V_inv[i][0]
    Pval = np.zeros((n, len(Z)))
    for j in range(len(Z)):
        for i in range(n):
            for h in range(s):
                Pval[i][j] = (Pval[i][j] + (i+1)**h * P[h][j]) % q
    return Pval

def GenHV(n, k, s, q):
    v1list = []
    for i in range(s):
        row =[int((j+1)**i) for j in range(s)]
        v1list.append(row)
    V1 = Matrix(v1list)
    V1_inv = matrix2numpy(V1.inv_mod(q))
    v2list = []
    for i in range(s):
        row =[int((j+1)**i) for j in range(s, n)]
        v2list.append(row)
    V2 = matrix2numpy(Matrix(v2list))
    P = np.matmul(V1_inv, V2) % q
    G = np.zeros((s, n))
    for i in range(s):
        for j in range(n):
            if j < s:
                if i == j: G[i][j] = 1 % q
            else:
                G[i][j] = P[i][j-s] % q
    # G = (I, P) is the generator matrix of the (n,s)-Reed-Solomon code in the standard form
    H = np.zeros((n-s, n))
    for i in range(n-s):
        for j in range(n):
            if j < s:
                H[i][j] = - P[j][i] % q
            else:
                if i == j-s: H[i][j] = 1 % q
    # H = (-P^T, I) is the parity-check matrix of the (n,s)-Reed-Solomon code
    return H, V1_inv
    
n = 3
t = 1
ID = 0

D = 1; eps = 0.1; delta = 10**(-8); p = 1 - np.exp(-eps/D)

q = 2**61-1
logq = math.ceil(math.log2(q))

lmdlist = [64, 128, 256]

Nmin = 4
Nmax = 10
logNmax = Nmax - Nmin + 1 
Nlist = [2**i for i in range(Nmin, Nmax + 1)]

# WAN
print("WAN")
latency = 0.02 # sec
bandwidth = 50 * 10**6 # bps

# # LAN
# print("LAN")
# latency = 0 # sec
# bandwidth = 1000 * 10**6 # bps

timePRF = 3.5 * 10**(-9)

time_semihonest = np.zeros((len(lmdlist), len(Nlist)))
time_baseline = np.zeros((len(lmdlist), len(Nlist)))
time_allsets = np.zeros((len(lmdlist), len(Nlist)))
time_partition = np.zeros((len(lmdlist), len(Nlist)))
time_binomial = np.zeros((len(lmdlist), len(Nlist)))

for lmd_id in range(len(lmdlist)):
    lmd = lmdlist[lmd_id]
    print("security parameter:"); print(lmd)
    for N_id in range(len(Nlist)):
        N = Nlist[N_id]
        
        msg = np.random.randint(20, size=(N, n))

        # semi-honest
        # noise generation
        m = n - t
        start = time.perf_counter()
        for rpt in range(N):
            noise = np.random.negative_binomial(1/m, p) - np.random.negative_binomial(1/m, p)
        end = time.perf_counter()
        time_semihonest[lmd_id][N_id] += end - start

        # reconstruction
        start = time.perf_counter()
        for rpt in range(N):
            result = 0
            for i in range(n):
                result = (result + msg[rpt][i]) % q
        end = time.perf_counter()
        time_semihonest[lmd_id][N_id] += end - start

        # baseline
        # key distribution
        k = n - t
        time_baseline[lmd_id][N_id] += 2*latency + lmd * k**2 * math.comb(n,t) / bandwidth

        # noise generation
        Z_baseline = list(itertools.combinations(range(n), k))
        M_baseline = GenM(n, k, Z_baseline)
        P_baseline = GenP(n, k, Z_baseline, q)
        start = time.perf_counter()
        for rpt in range(N):
            share = 0
            for j in range(len(Z_baseline)):
                if M_baseline[ID][j] == 1:
                    noise = np.random.geometric(p) - np.random.geometric(p)
                    share = (share + noise * P_baseline[ID][j]) % q 
        end = time.perf_counter()
        time_baseline[lmd_id][N_id] += N * len(Z_baseline) * timePRF + (end - start)

        # reconstruction
        s = n - k + 1
        H_baseline, V_baseline = GenHV(n, k, s, q)
        start = time.perf_counter()
        for rpt in range(N):
            result = 0
            parity = np.matmul(H_baseline, msg[rpt]) % q
            if all([parity[i]==0 for i in range(n-s)]):
                result = np.matmul(msg[rpt][0:s], V_baseline)[0] % q
        end = time.perf_counter()
        time_baseline[lmd_id][N_id] += latency + (n**2 * logq * N / bandwidth + (end - start))

        # allsets
        # key distribution
        k = t + 1
        time_allsets[lmd_id][N_id] += 2*latency + lmd * k**2 * math.comb(n,t+1) / bandwidth

        # noise generation
        m = math.comb(n-t, t+1)
        Z_allsets = list(itertools.combinations(range(n), k))
        M_allsets = GenM(n, k, Z_allsets)
        P_allsets = GenP(n, k, Z_allsets, q)
        start = time.perf_counter()
        for rpt in range(N):
            share = 0
            for j in range(len(Z_allsets)):
                if M_allsets[ID][j] == 1:
                    noise = np.random.negative_binomial(1/m, p) - np.random.negative_binomial(1/m, p)
                    share = (share + noise * P_allsets[ID][j]) % q 
        end = time.perf_counter()
        time_allsets[lmd_id][N_id] += N * len(Z_allsets) * timePRF + (end - start)

        # reconstruction
        s = n - k + 1
        H_allsets, V_allsets = GenHV(n, k, s, q)
        start = time.perf_counter()
        for rpt in range(N):
            result = 0
            parity = np.matmul(H_allsets, msg[rpt]) % q
            if all([parity[i]==0 for i in range(n-s)]):
                result = np.matmul(msg[rpt][0:s], V_allsets)[0] % q
        end = time.perf_counter()
        time_allsets[lmd_id][N_id] += latency + (n**2 * logq * N / bandwidth + (end - start))

        # partition
        # key distribution
        k = t + 1
        ell = int(math.lcm(n,t+1)/(t+1))
        time_partition[lmd_id][N_id] += 2*latency + lmd * k**2 * ell / bandwidth

        # noise generation
        m = math.ceil(ell*(1-t*(t+1)/n))
        Z_partition = []
        for j in range(ell):
            Z_partition.append([(j*k + i) % n for i in range(k)])
        M_partition = GenM(n, k, Z_partition)
        P_partition = GenP(n, k, Z_partition, q)
        start = time.perf_counter()
        for rpt in range(N):
            share = 0
            for j in range(len(Z_partition)):
                if M_partition[ID][j] == 1:
                    noise = np.random.negative_binomial(1/m, p) - np.random.negative_binomial(1/m, p)
                    share = (share + noise * P_partition[ID][j]) % q 
        end = time.perf_counter()
        time_partition[lmd_id][N_id] += N * len(Z_partition) * timePRF + (end - start)

        # reconstruction
        s = n - k + 1
        H_partition, V_partition = GenHV(n, k, s, q)
        start = time.perf_counter()
        for rpt in range(N):
            result = 0
            parity = np.matmul(H_partition, msg[rpt]) % q
            if all([parity[i]==0 for i in range(n-s)]):
                result = np.matmul(msg[rpt][0:s], V_partition)[0] % q
        end = time.perf_counter()
        time_partition[lmd_id][N_id] += latency + n**2 * logq * N / bandwidth + (end - start)
        
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

        # key distribution
        k = n - t
        time_binomial[lmd_id][N_id] += 2*latency + lmd * k**2 * math.comb(n,t) / bandwidth

        # noise generation
        Z_binomial = list(itertools.combinations(range(n), k))
        M_binomial = GenM(n, k, Z_binomial)
        P_binomial = GenP(n, k, Z_binomial, q)
        start = time.perf_counter()
        for rpt in range(N):
            share = 0
            for j in range(len(Z_binomial)):
                if M_binomial[ID][j] == 1:
                    noise = (1/M)*(np.random.binomial(K, 1/2) - K/2)
                    share = (share + noise * P_binomial[ID][j]) % q 
        end = time.perf_counter()
        time_binomial[lmd_id][N_id] += N * len(Z_binomial) * timePRF + (end - start)

        # reconstruction
        s = n - k + 1
        H_binomial, V_binomial = GenHV(n, k, s, q)
        start = time.perf_counter()
        for rpt in range(N):
            result = 0
            parity = np.matmul(H_binomial, msg[rpt]) % q
            if all([parity[i]==0 for i in range(n-s)]):
                result = np.matmul(msg[rpt][0:s], V_binomial)[0] % q
        end = time.perf_counter()
        time_binomial[lmd_id][N_id] += latency + n**2 * logq * N / bandwidth + (end - start)
    
    print("semi-honest:")
    for N_id in range(len(Nlist)):
        print(time_semihonest[lmd_id][N_id])
    print("baseline:")
    for N_id in range(len(Nlist)):
        print(time_baseline[lmd_id][N_id])
    print("allsets:")
    for N_id in range(len(Nlist)):
        print(time_allsets[lmd_id][N_id])
    print("partition:")
    for N_id in range(len(Nlist)):
        print(time_partition[lmd_id][N_id])
    print("binomial:")
    for N_id in range(len(Nlist)):
        print(time_binomial[lmd_id][N_id])
