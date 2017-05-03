import numpy as np
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
from matplotlib.patches import Rectangle

from scipy.stats import norm, ks_2samp
from scipy.integrate import quad
##import functions
import sys
sys.path.insert(0, 'C:\Python34\masters\graphs')
import functions as fun

degtorad = 0.0174533

###South Cluster
sclusterRA = [x*degtorad for x in fun.clustersouthRA]
sclusterdec = [x*degtorad for x in fun.clustersouthdec]
sclusterz = fun.clustersouthz
sclusterRA = [x-2*np.pi if x > np.pi else x for x in sclusterRA]

##North Cluster
nclusterRA = [float(x)*degtorad for x in fun.north_sc_RA]
nclusterdec = [float(x)*degtorad for x in fun.north_sc_dec]
nclusterz = [float(x) for x in fun.north_sc_z]
nclusterRA = [x-2*np.pi if x > np.pi else x for x in nclusterRA]

##Planck Cluster
pclusterRA = [float(x)*degtorad for x in fun.newclusterRA]
pclusterdec = [float(x)*degtorad for x in fun.newclusterdec]
pclusterz = fun.newclusterz
pclusterRA = [x-2*np.pi if x > np.pi else x for x in pclusterRA]

##Abell Cluster
aclusterRA = [float(x)*degtorad for x in fun.abelRA]
aclusterdec = [float(x)*degtorad for x in fun.abeldec]
aclusterz = fun.abelz
aclusterRA = [x-2*np.pi if x > np.pi else x for x in aclusterRA]


fig = plt.figure()#figsize=(10, 5))
ax = fig.add_subplot(111, projection="mollweide")
ax.grid(True)

i=0
ax.plot(sclusterRA, sclusterdec, '.', color='blue', label='ROSAT south')
ax.plot(nclusterRA, nclusterdec, '.', color='red', label='ROSAT north')
ax.plot(pclusterRA, pclusterdec, '.', color='grey', label='Planck SZ')
ax.plot(aclusterRA, aclusterdec, '.', color='green', label='Abell')
plt.xlabel('Right Ascenscion')
plt.ylabel('Declination')
plt.legend()
plt.show()
