from materials import *

from matplotlib.pyplot import *
from numpy import arange, log10

w = Material("Water")
we = w.electron
wee = we.elastic
web = we.brem


mu_c = wee.mu_c

E = 40e6
F = lambda s: wee.F(mu_c(E), E, s)




def U(a, b, x):
    return 1/(b-a) if a <= x <= b else 0

def dist(a, b, x):
    return a*U(0, b, x) + b*U(b, 1, x)
X = arange(0, 1, .01)
for s in arange(0, 50, 1):
    a, b = wee.F(mu_c(E), E, s)
    #print(a, b)
    
    Y = [dist(a, b, x) for x in X]

    plot(X, Y)
show()









##E = arange(0, 10e6, 10e3)
##IMFP = array([wee.mu_c(en) for en in E])
##
##wee_ = [1/wee(x) for x in E]
##plot(E, wee_)
##plot(E, 1/web(E))
##show()
