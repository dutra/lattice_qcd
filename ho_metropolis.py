from pylab import *
from time import time

tic = time()

#====================================================================================================
def update(x,repeat=1):
    for i in range(repeat):
        for j in range(1,N):
            old_x = x[j]                        # save original value
            old_Sj = S(j,x)
            x[j] += uniform(-eps,eps)     # update x[j]
            dS = S(j,x) - old_Sj                # change in action
            if dS > 0 and exp(-dS) < uniform(0,1):
                x[j] = old_x                    # restore old value


def S(j,x):         # harm. osc. S
    jp = (j+1) % N  # next site
    jm = (j-1) % N  # previous site
    return a * ( (1/2) * ((x[j+1]-x[j])/a)**2 + (1/2)*x[j]**2)


def compute_O(x):
    g = sum(S(j,x) for j in range(N))
    return g

def MCaverage(x,G):

    update(x, 10 * nrepeat)   # thermalize x

    for alpha in range(npaths):     # loop on random paths
        update(x, nrepeat)           # metropolis generate paths
        G[alpha] = compute_O(x)      # find total S for x

    avg = sum(exp(-G[alpha]) for alpha in range(npaths))
    return avg




N = 20 # number of time intervals
nrepeat = 10 # number of iterations for each path
npaths = 5000 # number of paths
a = 10# T/N
eps = 1.4



def start_x(u) :
    x = zeros(N+1, dtype=float)
    x[0]=x[-1]=u
    return x

poss = []
# do the simulation:
for i in arange(0,1,.1) :
    G = zeros(npaths, dtype=float)
    x = start_x(i)
    poss.append(MCaverage(x,G))
    print(poss[-1])

print(time() - tic)
plot(poss)
savefig("pic%s.png" % time.strftime("%m%d-%H%M"))
#show()
