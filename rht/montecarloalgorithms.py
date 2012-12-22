import scipy
import pylab
from pyqm import checktime
from integrator import leapfrog
t = checktime()

# list of MC algorithms
# metropolis-hasting
# hamiltonian monte carlo / hybrid monte carlo
# langevin
# heat bath
# overrelaxation


def generate_paths(path,npaths):
    paths = [path]
    U = lambda x: x*x
    gradU = lambda x: 2*x
    for i in range(npaths):
        newpath = hmc(U, gradU, paths[-1], .5, 200)
        paths.append(newpath)
    return paths


def hmc(U, gradU, current_q, epsilon, L=5):
    while 1:
        q = hmconce(U, gradU, current_q, epsilon, L=5)
        if isinstance(q, scipy.ndarray): # to ensure that a new path is generated
            return q

def hmconce(U, gradU, current_q, epsilon, L=5):
    # based on
    # http://www.cs.toronto.edu/~radford/ham-mcmc.abstract.html
    # http://www.cs.toronto.edu/~radford/ham-mcmc-simple
    # more, see http://deeplearning.net/tutorial/hmc.html
    """
    U         as Neal said "a function to evaluate minus the log of the density of the
    distribution to be sampled, plus any constant -> the potential energy"
    gradU     function, gradient of U
    epsilon   leapfrog stepsize
    L         number of leapfrog steps
    current_q current states (position variables only)
    """
    q = current_q
    length = len(q) if (isinstance(q, list) or isinstance(q,scipy.ndarray)) else 1
    p = scipy.randn(length)
    current_p = p


    p, q = leapfrog(gradU, p, q, epsilon, L)

    # negate momentum at end of trajectory to make the proposal symmetric
    p = -p

    # evaluate potential and kinetic at start and end of trajectory

    current_U = U(current_q)
    current_K = scipy.sum(current_p**2) / 2
    proposed_U = U(q)
    proposed_K = scipy.sum(p**2) / 2
    dH = proposed_U-current_U + proposed_K-current_K  

    # accept or reject the state at end of trajectory, returning either
    # the position at the end of the trajectory or the initial position

    if scipy.random.uniform() < scipy.exp(-dH):
        return scipy.copy(q)  # accept
    else:
        return False
        #return current_q  # reject

if __name__ == '__main__':
    paths = generate_paths(5,10)
    print paths


