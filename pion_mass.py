
from pylab import *

f = open("data/pion.dat","r")

c = zeros((739,256))

for i in arange(0,739) :
    for j in arange(0,256) :
        c[i][j] = f.readline().split(" ")[1]

m = zeros((739,256))
for i in arange(0,739):
    for j in arange(0,255):
        m[i][j] = log(c[i][j]/c[i][j+1])

avg = zeros(256)
for i in arange(0,256) :
    avg[i] = sum(m[j][i] for j in arange(0,256))/256

plot(avg)
show()
