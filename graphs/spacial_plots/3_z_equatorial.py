import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from scipy.stats import norm, ks_2samp
##import functions
import sys
sys.path.insert(0, 'C:\Python34\masters\graphs')
import functions as fun

##RA = [0], dec = [1], dist = [2], disterr = [3], redshifts = [4]
lowZ = fun.lowZ 
midZ = fun.midZ
highZ = fun.highZ
degtorad = 0.0174533

def stdev(x, mu, N):
        sumsq = 0
        i=0
        while i < N:
                sumsq = sumsq + (x[i]-mu)**2
                i+=1
        return np.sqrt((1/(N-1))*sumsq)

##Test against H0
lowZdeviation = []
midZdeviation = []
highZdeviation = []

i = 0
while i < len(lowZ):
        lowZdeviation.append(( fun.D(fun.HubbleIntegrate( fun.extract(lowZ, 4)[i] )) - fun.extract(lowZ, 5)[i])/fun.D(fun.HubbleIntegrate( fun.extract(lowZ, 4)[i] )))
        i+=1

i = 0
while i < len(midZ):
        midZdeviation.append(( fun.D(fun.HubbleIntegrate( fun.extract(midZ, 4)[i] )) - fun.extract(midZ, 5)[i])/fun.D(fun.HubbleIntegrate( fun.extract(midZ, 4)[i] )))
        i+=1 ##(D_L(z) - D_L) / D_L(z)

i = 0
while i < len(highZ):
        highZdeviation.append(( fun.D(fun.HubbleIntegrate( fun.extract(highZ, 4)[i] )) - fun.extract(highZ, 5)[i])/fun.D(fun.HubbleIntegrate( fun.extract(highZ, 4)[i] )))
        i+=1 ##(D_L(z) - D_L) / D_L(z)

'''
fig = plt.figure()
##x, y, size of marker, color gradients from deviation of H0, min gradient, max gradient, color scheme.
plt.scatter([(x*degtorad)-2*np.pi if x*degtorad > np.pi else x*degtorad for x in fun.extract(lowZ, 0)],
           [x*degtorad for x in fun.extract(lowZ, 1)],
           s=10, c=lowZdeviation, vmin=-0.05, vmax=0.05, cmap=plt.cm.seismic)
plt.colorbar()
plt.grid(True)
plt.xlabel('RA')
plt.ylabel('DEC')
plt.title('Deviation from expected distance modulus for z<0.1')
plt.show()


fig = plt.figure()
ax = fig.add_subplot(111, projection="mollweide")
##x, y, size of marker, color gradients from deviation of H0, min gradient, max gradient, color scheme.
ax.scatter([(x*degtorad)-2*np.pi if x*degtorad > np.pi else x*degtorad for x in fun.extract(midZ, 0)],
           [x*degtorad for x in fun.extract(midZ, 1)],
           s=10, c=midZdeviation, vmin=-0.05, vmax=0.05, cmap=plt.cm.seismic)
ax.grid(True)
plt.xlabel('RA')
plt.ylabel('DEC')
plt.title('Deviation from expected distance modulus for 0.1<z<0.5')
plt.show()


fig = plt.figure()
ax = fig.add_subplot(111, projection="mollweide")
##x, y, size of marker, color gradients from deviation of H0, min gradient, max gradient, color scheme.
ax.scatter([(x*degtorad)-2*np.pi if x*degtorad > np.pi else x*degtorad for x in fun.extract(highZ, 0)],
           [x*degtorad for x in fun.extract(highZ, 1)],
           s=10, c=highZdeviation, vmin=-0.05, vmax=0.05, cmap=plt.cm.seismic)
ax.grid(True)
plt.xlabel('RA')
plt.ylabel('DEC')
plt.title('Deviation from expected distance modulus for 0.5<z')
plt.show()
'''

lowZhist = plt.hist(lowZdeviation, bins=40, range=[-0.5, 0.5])
(lowZmu, lowZsigma) = norm.fit(lowZdeviation)
lowZSE = lowZsigma/np.sqrt(len(lowZdeviation))
emptySGC = plt.hist([], range=[-0.5, 0.5], alpha=0.5, label='$N=%i$\n$\mu=%.3f \pm %.3f$\n$\sigma=%.3f$'%(len(lowZdeviation), lowZmu, lowZSE, lowZsigma), color='white')
plt.xlabel('$\delta d_L$')
plt.ylabel('Frequency')
plt.title('Histogram of deviation from expected distance modulus for z<0.1')
plt.legend()
plt.show()


midZhist = plt.hist(midZdeviation, bins=40, range=[-0.5, 0.5])
(midZmu, midZsigma) = norm.fit(midZdeviation)
midZSE = midZsigma/np.sqrt(len(midZdeviation))
emptySGC = plt.hist([], range=[-0.5, 0.5], alpha=0.5, label='$N=%i$\n$\mu=%.3f \pm %.3f$\n$\sigma=%.3f$'%(len(midZdeviation), midZmu, midZSE, midZsigma), color='white')
plt.xlabel('$\delta d_L$')
plt.ylabel('Frequency')
plt.title('Histogram of deviation from expected distance modulus for 0.1<z<0.5')
plt.legend()
plt.show()


highZhist = plt.hist(highZdeviation, bins=20, range=[-0.5, 0.5])
(highZmu, highZsigma) = norm.fit(highZdeviation)
highZSE = highZsigma/np.sqrt(len(highZdeviation))
emptySGC = plt.hist([], range=[-0.5, 0.5], alpha=0.5, label='$N=%i$\n$\mu=%.3f \pm %.3f$\n$\sigma=%.3f$'%(len(highZdeviation), highZmu, midZSE, highZsigma), color='white')
plt.xlabel('$\delta d_L$')
plt.ylabel('Frequency')
plt.title('Histogram of deviation from expected distance modulus for 0.5<z')
plt.legend()
plt.show()
