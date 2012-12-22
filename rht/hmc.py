from pylab import *
# based on https://github.com/koepsell/pyhmc
# http://prb.aps.org/abstract/PRB/v45/i2/p679_1

"""Hybrid Monte Carlo Sampling

This package is a straight-forward port of the functions hmc2.m and
hmc2_opt.m from the MCMCstuff matlab toolbox written by Aki Vehtari
<http://www.lce.hut.fi/research/mm/mcmcstuff/>.
   
The code is originally based on the functions hmc.m from the netlab toolbox
written by Ian T Nabney <http://www.ncrg.aston.ac.uk/netlab/index.php>.

The portion of algorithm involving "windows" is derived from the C code for
this function included in the Software for Flexible Bayesian Modeling
written by Radford Neal <http://www.cs.toronto.edu/~radford/fbm.software.html>.
   
This software is distributed under the BSD License (see LICENSE file).

Authors
-------
- Kilian Koepsell <kilian@berkeley.edu>
"""



#-----------------------------------------------------------------------------
# Functions
#-----------------------------------------------------------------------------
#1. scipy.optimize.check_grad(func, grad, x0, *args)
"""Check the correctness of a gradient function by comparing it against a
    (forward) finite-difference approximation of the gradient."""


def hmc(f, x, gradf,
        args=(),
        steps=1,
        nsamples=1,
        nomit=0,
        stepadj=0.2,
        window=1,
        ):
    """Hybrid Monte Carlo sampling.

    Uses a hybrid Monte Carlo algorithm to sample from the distribution P ~
    exp(-f), where f is the first argument to hmc. The Markov chain starts
    at the point x, and the function gradf is the gradient of the `energy'
    function f.

    Parameters
    ----------

    f : function
      Energy function.
    x : 1-d array
      Starting point for the sampling Markov chain.
    gradf : function
      Gradient of the energy function f.

    Optional Parameters
    -------------------

    steps : int
      Defines the trajectory length (i.e. the number of leapfrog
      steps at each iteration).
    nsamples : int
      The number of samples retained from the Markov chain.
    nomit : int
      The number of samples omitted from the start of the chain
    stepadj : float
      The step adjustment used in leap-frogs. Default: 0.2
    args : tuple
      additional arguments to be passed to f() and gradf().

    Returns
    -------

    samples : array
      Array with data samples in rows.

    """
    if window > steps:
        window = steps
        print "setting window size to step size %d" % window

    nparams = len(x)
    epsilon = stepadj

    # Initialize matrix of returned samples
    samples = zeros((nsamples, nparams))


    # Initialise momenta from stored state
    p = randn(nparams)
        
    # Main loop.
    all_args = [f,
                x,
                gradf,
                args,
                p,
                samples,
                nsamples,
                nomit,
                window,
                steps,
                epsilon]

    nreject = hmc_main_loop(*all_args)
    print '\nFraction of samples rejected:  %g\n'%(nreject/float(nsamples))
    return samples
    

def hmc_main_loop(f, x, gradf, args, p, samples,
                  nsamples, nomit, window, steps,
                  epsilon):
    nparams = len(x)
    nreject = 0              # number of rejected samples
    window_offset = 0        # window offset initialised to zero
    k = -nomit       # nomit samples are omitted, so we store
    
    # Evaluate starting energy.
    E = f(x, args)

    while k < nsamples:  # samples from k >= 0
        # Store starting position and momenta
        xold = x
        pold = p
        # Recalculate Hamiltonian as momenta have changed
        Eold = E
        Hold = E + 0.5*(p**2).sum()

        # Decide on window offset, if windowed HMC is used
        if window > 1:
            window_offset = int(window*rand())

        have_rej = 0
        have_acc = 0
        n = window_offset
        direction = -1 # the default value for direction 
                       # assumes that windowing is used

        while direction == -1 or n != steps:
            # if windowing is not used or we have allready taken
            # window_offset steps backwards...
            if direction == -1 and n==0:
                # Restore, next state should be original start state.
                if window_offset > 0:
                    x = xold
                    p = pold
                    n = window_offset

                # set direction for forward steps
                E = Eold
                H = Hold
                direction = 1
                stps = direction
            else:
                if n*direction+1<window or n > (steps-window):
                    # State in the accept and/or reject window.
                    stps = direction
                else:
                    # State not in the accept and/or reject window. 
                    stps = steps-2*(window-1)

                # First half-step of leapfrog.
                p = p - direction*0.5*epsilon*gradf(x, *args)
                x = x + direction*epsilon*p
                
                # Full leapfrog steps.
                for m in range(abs(stps)-1):
                    p = p - direction*epsilon*gradf(x, *args)
                    x = x + direction*epsilon*p

                # Final half-step of leapfrog.
                p = p - direction*0.5*epsilon*gradf(x, *args)

                E = f(x, *args)
                H = E + 0.5*(p**2).sum()

                n += stps

            if window != steps+1 and n < window:
                # Account for state in reject window.  Reject window can be
                # ignored if windows consist of the entire trajectory.
                if not have_rej:
                    rej_free_energy = H
                else:
                    rej_free_energy = -addlogs(-rej_free_energy, -H)

                if not have_rej or rand() < exp(rej_free_energy-H):
                    E_rej = E
                    x_rej = x
                    p_rej = p
                    have_rej = 1

            if n > (steps-window):
                # Account for state in the accept window.
                if not have_acc:
                    acc_free_energy = H
                else:
                    acc_free_energy = -addlogs(-acc_free_energy, -H)

                if not have_acc or  rand() < exp(acc_free_energy-H):
                    E_acc = E
                    x_acc = x
                    p_acc = p
                    have_acc = 1
  
        # Acceptance threshold.
        a = exp(rej_free_energy - acc_free_energy)


        print 'New position is\n',x

        # Take new state from the appropriate window.
        if a > rand():
            # Accept 
            E = E_acc
            x = x_acc
            p = -p_acc # Reverse momenta
            print 'Finished step %4d  Threshold: %g\n'%(k,a)
        else:
            # Reject
            if k >= 0:
                nreject = nreject + 1

            E = E_rej
            x = x_rej
            p = p_rej
            print '  Sample rejected %4d.  Threshold: %g\n'%(k,a)

        if k >= 0:
            # Store sample
            samples[k,:] = x;

        # Set momenta for next iteration
        # Replace all momenta
        p = randn(nparams)

        k += 1

    return nreject


def addlogs(a,b):
    """Add numbers represented by their logarithms.
    
            Description
            Add numbers represented by their logarithms.
            Computes log(exp(a)+exp(b)) in such a fashion that it 
            works even when a and b have large magnitude.
    """
    
    if a>b:
        return a + log(1+exp(b-a))
    else:
        return b + log(1+exp(a-b))


if __name__ == '__main__':
    from time import time as now
    
    def f_multivariate_normal(x,M):
        """Energy function for multivariate normal distribution
        """
        return .5*dot(dot(x,M),x)

    def g_multivariate_normal(x,M):
        """Energy gradient for multivariate normal distribution
        """
        return .5*dot(x,M+M.T)

    # sample from 5-dimensional unit variance gaussian
    dim = 5
    M = eye(dim)
    x0 = zeros(dim)    
    t0 = now()
    samples = hmc(f_multivariate_normal, x0, g_multivariate_normal, args=(M,),
                  nsamples=10**3, nomit=10**3, steps=100, stepadj=.05, )
    dt = now()-t0

    # check covariance matrix of sampled data
    C = cov(samples,rowvar=0,bias=1)

    print "mean squared error of sample covariance (expect .001): %f" % ((C-linalg.inv(M))**2).mean()
    print "time (python): ",dt


