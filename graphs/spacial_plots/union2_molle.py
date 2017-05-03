import matplotlib.pyplot as plt
import numpy as np
from astropy import units as u
from astropy.coordinates import SkyCoord
import scipy.stats
##import functions
import sys
sys.path.insert(0, 'C:\Python34\masters\graphs')
import functions as fun

degtorad = 0.0174533

clusterRA = [x*degtorad for x in fun.RA]
clusterdec = [x*degtorad for x in fun.dec]
clusterRA = [x-2*np.pi if x > np.pi else x for x in clusterRA]
z = fun.redshifts
name = fun.name

i = 0
while i < len(z):
        if z[i] > 0.2:
                del clusterRA[i]
                del clusterdec[i]
                del z[i]
                del name[i]
        else:
                i+=1


#########EQUITORIAL COORDINATE

fig = plt.figure()#figsize=(10, 5))
ax = fig.add_subplot(111, projection="mollweide")
ax.grid(True)

ax.plot(clusterRA, clusterdec, '.', color='blue')
plt.xlabel('RA')
plt.ylabel('DEC')
plt.title('Equitorial Coordinates of Planck cluster objects')
plt.legend()
plt.show()


#####GALACTIC COORDINATE TRANSFORMATION#######
NgalRA = []
Ngaldec = []
i = 0
while i < len(clusterRA):
        x = SkyCoord(ra=clusterRA[i]*u.radian, dec=clusterdec[i]*u.radian, frame='icrs')
        l = x.galactic.l.radian
        if l > np.pi:
                l = l-2*np.pi
        NgalRA.append(l)
        Ngaldec.append(x.galactic.b.radian)
        i+=1



fig = plt.figure()#figsize=(10, 5))
ax = fig.add_subplot(111, projection="mollweide")
ax.grid(True)

i=0
ax.plot(NgalRA, Ngaldec, '.', color='blue')
while i < len(NgalRA):
    ax.annotate(name[i], xy=(NgalRA[i], Ngaldec[i]), textcoords='data', color='red')
    i+=1
plt.xlabel('Galactic Longitude')
plt.ylabel('Galactic Latitude')
plt.title('Galactic Coordinates of Union2.1 Supernovae')
plt.legend()
plt.show()
