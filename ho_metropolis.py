from pylab import *
from time import *
from math import pi

tic = time()
acceptance = 0
total = 0

#====================================================================================================
def update(x,repeat=1):
    global acceptance
    global total
    for i in range(repeat):
        for j in range(1,N):
            total+=1
            old_x = x[j]                        # save original value
            old_Sj = S(j,x)
            x[j] += uniform(-eps,eps)     # update x[j]
            dS = S(j,x) - old_Sj                # change in action
            acceptance+=1
            if dS > 0 and exp(-dS) < uniform(0,1):
                x[j] = old_x                    # restore old value
                acceptance-=1


def S(j,x):         # harm. osc. S
    jp = (j+1) % N  # next site
    jm = (j-1) % N  # previous site
    return a * ( (1/2) * ((x[j+1]-x[j])/a)**2 + (1/2)*x[j]**2)


def compute_O(x):
    g = exp(-sum(S(j,x) for j in range(N)))
    return g

def MCaverage(x,G):

    update(x, 10 * nrepeat)   # thermalize x

    for alpha in range(npaths):     # loop on random paths
        update(x, nrepeat)           # metropolis generate paths
        G[alpha] = compute_O(x)      # find total S for x


    avg = sum(G[alpha] for alpha in range(npaths))
    return avg




N = 10 # number of time intervals ~ a
nrepeat = 10 # number of iterations for each path #N_cor # 1/sqrt(a)
npaths = 1000 # number of paths
a = 10 # T/N 1/nrepeat^2
eps = 0.5 # 1/sqrt(a)



def start_x(u) :
    x = zeros(N+1, dtype=float)
    x[0]=x[-1]=u
    return x

def avg(G) :
    return sum(G)/len(G)
def sdev(G) :
    v = asarray(G)
    return absolute(avg(v**2)-avg(v)**2)**0.5

poss = []


# do the simulation:
for i in arange(0,1,0.1) :
    aposs = []
    for j in range(5) :

        G = zeros(npaths, dtype=float)
        x = start_x(i)
        aposs.append(MCaverage(x,G))
        print(i,aposs[-1],acceptance/total*100)
        total = 0
        acceptance = 0

    print("-----------------")
    poss.append(avg(aposs))
    print("i: {0}, pos: {1}, sdev: {2}".format(i,poss[-1],sdev(aposs)))
    print("-----------------")


print(time() - tic)
plot(poss,'bo')
savefig("pic{0}-N{1}_nr{2}_np{3}_a{4}_e{5:.1f}.png".format(strftime("%m%d-%H%M"),N,nrepeat,npaths,a,eps))
show()
