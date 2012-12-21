from pylab import *
ion()

#f = open("data/pion.dat","r")
times, c = loadtxt("data/px0py0pz0_pi_Nsrc174_Ncfg739_32x256_um0p0840_sm0p0743_P.dat",
        unpack = 1, skiprows=1)

num = len(times)/256
times = times.reshape((num,256))
c = c.reshape((num,256))
#for i in range(4):
    #figure()
    ##plot(times[i], log(c[i] / roll(c[i],-1)))
    #plot(times[i], log(c[i] / c[i+1]))
    ##plot(times[i], log(c[i]))

cmean = c.mean(axis=0)
cstd = c.std(axis=0) / sqrt(len(c))

#print shape(cmean)
#print shape(times[0])
errorbar(times[0], cmean, cstd)


raw_input()
