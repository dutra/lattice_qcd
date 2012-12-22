import scipy
import scipy.integrate as integrate

#list of odeint implementations
# http://hplgit.github.com/odespy/doc/api/odespy.html
# http://github.com/olivierverdier/odelab
# http://github.com/jtksai/PyCOOL
# http://github.com/simon-r/PyParticles

#about higher order symplectic integrator
#http://iopscience.iop.org/1538-3881/119/1/425/fulltext/990231.text.html

#state-of-the-art of ode as of dec 8 2012
# http://lh3lh3.users.sourceforge.net/solveode.shtml
# nonstiff classic: RK4 good overall, easy to implement
# stiff classic: rosenbrock (easy), bulirsch-stoer-bader-deuflhard
# best known ode solver: LSODE (evolves into ODEPACK) and VODE


def leapfrog(gradU, p, q, epsilon, L):
    # L: number of leapfrog steps
    # epsilon: leapfrog stepsize

    # make a half step at the beginning
    p -= epsilon * gradU(q) / 2
    q += epsilon * p
    
    # alternate full steps for position and momentum
    for i in range(L-1):
        p += -epsilon * gradU(q)
        q += epsilon * p

    # make a half step for momentum at the end
    p -= epsilon * gradU(q) / 2

    return p, q

def leapfrogarray(gradU, p, q, epsilon, L):
    # make a half step at the beginning
    p -= epsilon * gradU(q) / 2
    ps = [p]
    qs = [q]
    qs.append(q + epsilon * p)
    
    # alternate full steps for position and momentum
    for i in range(L-1):
        ps.append(ps[-1] - epsilon * gradU(qs[-1]))
        qs.append(qs[-1] + epsilon * ps[-1])
    # make a half step for momentum at the end
    ps.append( -epsilon * gradU(qs[-1]) / 2)

    return ps, qs

def normaleulerintegratortest(gradU, p, q, epsilon, L):
    # straightforward
    ps = [p]
    qs = [q]
    for i in range(L):
        dp = - epsilon * gradU(qs[-1])
        dq = epsilon * ps[-1]
        ps.append(ps[-1] + dp)
        qs.append(qs[-1] + dq)
    return ps, qs


def semiimpliciteulerintegratortest(gradU, p, q, dt, L):
    ps = [p]
    qs = [q]
    for i in range(L):
        ps.append(ps[-1] - dt * gradU(qs[-1]))
        # here qs uses the updated ps
        qs.append(qs[-1] + dt * ps[-1])
    return ps, qs


def odeinttest(gradU, p, q, dt, L):
    # LSODA, automatically switch between stiff and nonstiff
    def evolve(pq,t=0):
        p, q = scipy.split(pq,2)
        dp = -dt * gradU(q)
        dq = dt * p
        return scipy.concatenate((dp, dq))
   
    ps, qs = scipy.integrate.odeint(evolve, [p,q], range(L)).T

    return ps, qs


def symplecticintegratortest(gradU, p, q, dt, L):
    # http://en.wikipedia.org/wiki/Symplectic_integrator
    ps, qs = [p], [q]

    #ruth 3rd order
    #ruth3c = [2./3, -2./3, 1]
    #ruth3d = [7./24, 3./4, -1./24]

    #ruth 4th order
    ruth4common = 2 - 2**1./3
    c1 = c4 = 1. / 2 / ruth4common
    c2 = c3 = (1.- 2**1./3) / 2 / ruth4common
    d1 = d3 = 1. / ruth4common
    d2 = -d1 * 2**1./3
    d4 = 0
    ruth4c = [c1,c2,c3,c4]
    ruth4d = [d1,d2,d3,d4]

    #semi-implicit euler
    #dp = -dt * gradU(q)
    #dq = dt * (p+dp)

    #verlet 2nd order
    #dp1 = -dt * gradU(q) * .5
    #dq1 = dt * (p+dp1)
    #dq2 = -dt * gradU(q+dq1) *.5
    #dp = dp1
    #dq = dq1+dq2

    
    for i in range(L):
        p, q = ps[-1], qs[-1]

        #ruth 3rd order
        #dp1 = -dt * gradU(q) * ruth3c[0]
        #dq1 = dt * (p+dp1) * ruth3d[0]
        #dp2 = -dt * gradU(q+dq1) * ruth3c[1]
        #dq2 = dt * (p+dp1+dp2) * ruth3d[1]
        #dp3 = -dt * gradU(q+dq1+dq2) * ruth3c[2]
        #dq3 = dt * (p+dp1+dp2+dp3) * ruth3d[2]
        #dp = dp1+dp2+dp3
        #dq = dq1+dq2+dq3

        #ruth 4th order
        dp1 = -dt * gradU(q) * ruth4c[0]
        dq1 = dt * (p+dp1) * ruth4d[0]
        dp2 = -dt * gradU(q+dq1) * ruth4c[1]
        dq2 = dt * (p+dp1+dp2) * ruth4d[1]
        dp3 = -dt * gradU(q+dq1+dq2) * ruth4c[2]
        dq3 = dt * (p+dp1+dp2+dp3) * ruth4d[2]
        dp4 = -dt * gradU(q+dq1+dq2+dq3) * ruth4c[3]
        dq4 = dt * (p+dp1+dp2+dp3+dp4) * ruth4d[3]
        dp = dp1+dp2+dp3+dp4
        dq = dq1+dq2+dq3+dq4
        ps.append(p+dp)
        qs.append(q+dq)
    return ps, qs


# for testing
if __name__ == '__main__':
    # testing for the sanity of leapfrog method and other methods
    p, q = scipy.randn(2)
    pylab.ion()

    #verified to be symplectic
    #ps,qs = leapfrogintegratortest(lambda x: x, p, q, .1, 200)
    #pylab.plot(qs)

    #not symplectic, as expected
    #psnormal,qsnormal = normaleulerintegratortest(lambda x: x, p, q, .1, 200)
    #pylab.plot(qsnormal)

    # symplectic
    pssemi,qssemi = semiimpliciteulerintegratortest(lambda x: x, p, q, .1, 200)
    pylab.plot(qssemi)
    t.next()

    psodeint,qsodeint = odeinttest(lambda x: x, p, q, .1, 200)
    pylab.plot(qsodeint)
    t.next()

    # supposed to be symplectic
    pssym,qssym = symplecticintegratortest(lambda x: x, p, q, .1, 200)
    pylab.plot(qssym)
    t.next()
    #pylab.legend(['leapfrog', 'normaleuler', 'semiimplicit','symplectic'])
    raw_input()

