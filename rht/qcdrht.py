# -*- coding: utf-8 -*-


# The program below is extracted from the following reference:

# Lepage, G. Peter. 1998. Lattice QCD for novices. In <em>Strong Interactions
# at Low and Intermediate Energies: Proceedings of the 13th Annual HUGS
# at CEBAF, Jefferson Laboratory</em> (edited by  J. L. Goity), pp. 49–90.
# River Edge, N.J.: World Scientific.
# Preprint available at arxiv.org/abs/hep-lat/0506036

# To understand what the code does and how to interpret the results, please
# see the Lepage article.

# The only changes made in the published program are minor modifications
# needed to make it compatible with more-recent versions of Python. In
# particular, the older Numeric package has been replaced with the modern
# NumPy package.

# tweaked by rht Apr 26 2012


from pylab import *
from time import time
#from numba import double, int_
#from numba import autojit
#from numba.decorators import jit

tic = time()

#@jit(arg_types=[double[:], int_])
#@autojit(backend="ast")


def S(j,x):         # harm. osc. S 
    jp = (j+1) % N  # next site 
    jm = (j-1) % N  # previous site 
    return a*x[j]**2/2 + x[j]*(x[j]-x[jp]-x[jm])/a 

def compute_G(x): 
    g = [sum([x[j] * x[(j + n) % N] for j in range(N)])/N for n in range(N)]
    return g 

def MCaverage(x,G): 
    def update(x,repeat): 
        for i in range(repeat):
            #update once
            for j in range(1,N):
                old_x = x[j]                        # save original value 
                old_Sj = S(j,x) 
                x[j] += uniform(-eps,eps)     # update x[j] 
                if dS>0 and (exp(-(S(j,x) - old_Sj)) < random()):
                    x[j] = old_x                    # restore old value
 
    update(x, 5 * N_cor)   # thermalize x
    for alpha in range(N_cf):     # loop on random paths 
        update(x, N_cor) 
        G[alpha] = compute_G(x)
    for n in range(N):            # compute MC averages 
        avg_G = sum(G[alpha][n] for alpha in range(N_cf)) / N_cf
        #print("G(%d) = %g" % (n,avg_G) )
    return x, avg_G


def bootstrap(G):
    N_cf = len(G) 
    # choose random config from G ensemble
    G_bootstrap = [G[randint(N_cf)] for i in range(N_cf)]
    #G_bootstrap = array(G)[randint(N_cf, size=N_cf)]
    return G_bootstrap 


def bin(G,binsize): 
    G_binned = []                       # binned ensemble 
    for _ in range(len(G),binsize):   # loop on bins 
        G_avg = sum([G[_ + __] for __ in range(binsize)])  # loop on bin elements
        G_binned.append(G_avg/binsize)  # keep bin avg 
    return G_binned 


# set parameters: 

N = 20 
N_cor = 50 
N_cf = 100
a = 0.5 
eps = 1.4 

# create arrays:
x = zeros(N) 
G = zeros((N_cf,N)) 
# do the simulation: 

x, avg_G = MCaverage(x,G) 

# To test the binning and bootstrap codes add the following
# to the the file: 

print('avg G\n', G.mean(axis=0) )
print('avg G (binned)\n', mean(bin(G,4),axis=0) )
print('avg G (bootstrap)\n', mean(bootstrap(G),axis=0) )

# The average of the binned copy of G should be the same as the 
# average of G itself; the average of the bootstrap copy should 
# be different by an amount of order the Monte Carlo error. 
# Compute averages for several bootstrap copies to get a good 
# feel for the errors. Finally one wants to extract energies. 
# This is done by adding code to compute ∆E(t): 

def deltaE(G):      # Delta E(t)
    avgG = mean(G, axis=0)
    adE = log(absolute(avgG[:-1]/avgG[1:]))
    return adE/a 

print('Delta E\n',deltaE(G) )
# print 'Delta E (bootstrap)\n',deltaE(bootstrap(G)) 

# Again repeating the evaluation for 50 or 100 bootstrap 
# copies of G gives an estimate of the statistical errors 
# in the energies. Additional code can be added to evaluate 
# standard deviations from these copies: 

def bootstrap_deltaE(G,nbstrap=100):    # Delta E + errors 
    avgE = deltaE(G)                # avg deltaE 
    bsE = [] 
    for i in range(nbstrap):    # bs copies of deltaE 
        g = bootstrap(G) 
        bsE.append(deltaE(g)) 
    bsE = array(bsE) 
    sdevE = bsE.std(axis=0)               # spread of deltaE's 
    print("\n%2s %10s %10s" % ("t","Delta E(t)","error") )
    print(2 *"-" )
    for i in range(int(len(avgE)/2)): 
        print("%2d %10g %10g" % (i,avgE[i],sdevE[i]) )

bootstrap_deltaE(G) 

print("time %f" %(time()-tic))
