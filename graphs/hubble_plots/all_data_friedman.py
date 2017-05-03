##import graphing
import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import quad
##import functions
import sys
sys.path.insert(0, 'C:\Python34\masters\graphs')
import functions as fun

redshifts = fun.redshifts
distmod = fun.distmod
dist = fun.dist
disterror = fun.disterror


##Plot Hubble
z = np.arange(min(redshifts), max(redshifts), 0.001)
expD = [fun.expectedD(z) for z in z] ##mx + c
expH = [fun.HubbleIntegrate(z) for z in z]

##Graph for redshift vs distance
plt.figure()
plt.errorbar(redshifts, dist, yerr=disterror, fmt='o', markersize=1, elinewidth=1)
plt.plot(z, expH, ls='-', color='black')


plt.xlabel('Redshift')
plt.ylabel('Distance Modulus')
plt.title('Best fit for Union2.1 data')
plt.show()
