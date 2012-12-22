from pylab import *

data = loadtxt("wf.txt", unpack=1, delimiter=',')
plot(*data)
savefig("wf.png")
