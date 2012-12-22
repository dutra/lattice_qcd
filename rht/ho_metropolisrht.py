from pylab import *
from time import time, strftime

tic = time()
####
#parameters
#periodic boundary condition
N = 20 # number of time intervals
nrepeat = 10 # number of iterations for each path
npaths = 50 # number of paths
a = 10# T/N
eps = 1.4
####

# observables
def propagator(x):
    return exp(-sum(S(j,x) for j in range(N)))


def twopointcorrelator(x):
    return array([sum(x[j] * x[(j+n)%N] for j in range(N)) for n in range(N)])


def S(j,x):         # harm. osc. local action
    jn = (j+1) % N  # next site
    jp = (j-1) % N  # previous site
    return a * ( .5 * ((x[jn]-x[j])/a)**2 + .5*x[j]**2)
    return a * ( .5 * ((x[jn]-x[jp])/a)**2 + .5*x[j]**2)
############


def generate_thermalized_paths(x):
    accept = 0
    total = 0
    def update(x,repeat=1):
        for i in range(repeat):
            for j in range(1,N):
                old_x = x[j]                        # save original value
                old_Sj = S(j,x)
                x[j] += uniform(-eps,eps)     # update x[j]
                dS = S(j,x) - old_Sj                # change in action
                total += 1
                accept += 1
                if dS > 0 and exp(-dS) < uniform(0,1):
                    x[j] = old_x                    # restore old value
                    acceptance -= 1
    paths = []
    update(x, 10 * nrepeat)   # thermalize x
    for alpha in range(npaths):     # loop on random paths
        update(x, nrepeat)           # metropolis generate paths
        paths.append(x)
    acceptance = accept * 100. / total
    print("acceptance %f" %acceptance)
    return paths


def thermal_average(observable, paths):
    # here, axis=0 is used so that the observable can be a tuple in general
    return array([observable(x) for x in paths]).mean(axis=0)


def start_x(u) :
    x = zeros(N+1, dtype=float)
    x[0]=x[-1]=u
    return x


# do the simulation:
# numerical experiments
#1 propagator
#1.1 change a, while keeping T and N fixed.
X = arange(-1, 1, .1)
#for a in arange(.2,10,2):
poss = [thermal_average(propagator, generate_thermalized_paths(start_x(i)))
        for i in X]
plot(poss)
savefig("pic%s.png" % strftime("%m%d-%H%M"))

#2 two point correlation function
paths = generate_thermalized_paths(start_x(0))
corr = thermal_average(twopointcorrelator, paths)
def bootstrap(x):
    N = len(x)
    return array(x)[randint(N, size=N)]
def deltaE(corr):
    return log(abs(corr[:-1] / corr[1:])) / a
#def bootstrap_deltaE(nbstrap=100):
    #deltaEs = delta
print(deltaE(corr))
    

print("total time spent %f" % (time() - tic))
#show()
