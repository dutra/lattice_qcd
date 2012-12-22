from pylab import *
from scipy.optimize import curve_fit
ion()

#f = open("data/pion.dat","r")
times, c = loadtxt("data/px0py0pz0_pi_Nsrc174_Ncfg739_32x256_um0p0840_sm0p0743_P.dat",
        unpack = 1, skiprows=1)

num = len(times)/256
times = times.reshape((num,256))
c = c.reshape((num,256))
for i in range(4):
    figure()
    #plot(times[i], log(c[i] / roll(c[i],-1)))
    #plot(times[i], log(c[i] / c[i+1]))
    plot(times[i], log(c[i]))

cmean = c.mean(axis=0)
cstd = c.std(axis=0) / sqrt(len(c))
time = times[0]

#curve fitting
# see http://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.curve_fit.html
NT = time[-1]
def func(x, a,E):
    return 2*a * exp(-NT*E/2) * cosh((NT/2-x)*E)

popt, pcov = curve_fit(func, time, cmean, sigma=cstd)



#plotting
figure()
errorbar(time, cmean, cstd)
plot(time, func(time,*popt))

#calculating chi square
cfitted = func(time, *popt)
reducedchi2 = sum(((cmean - cfitted)/cstd)**2 / (cmean.size - 2))
print reducedchi2

raw_input()
