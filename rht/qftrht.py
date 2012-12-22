from pylab import *
from time import time, strftime
from pyqm import K1D
from montecarloalgorithms import hmc



tic = time()
####
#parameters
#periodic boundary condition
N = 20 # number of time intervals
nrepeat = 10 # number of iterations for each path
npaths = 20 # number of paths
a = 10# T/N
eps = 1.4
####

def propagator(x):
    return exp(-sum(S(j,x) for j in range(N)))


def twopointcorrelator(x):
    return array([sum(x[j] * x[(j+n)%N] for j in range(N)) for n in range(N)])


def S(x):
    k = K1D(a, N)
    return a * .5 * (dot(-x, ravel(dot(k, x))) + dot(x,x))


def generate_thermalized_paths(x):
    def gradH(q):
        "harmonic oscillator gradient"
        return q
    paths = [x]
    for i in range(50):
        print paths[-1]
        paths.append(hmc(S, gradH, paths[-1], a, L=nrepeat))
    return paths


def thermal_average(observable, paths):
    # here, axis=0 is used so that the observable can be a composite object in general
    return array([observable(x) for x in paths]).mean(axis=0)


def start_x(u) :
    x = zeros(N, dtype=float)
    x[0] = x[-1] = u
    return x


# do the simulation:
# numerical experiments
#1 propagator
#1.1 change a, while keeping T and N fixed.
X = r_[-1:1:-1]
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
