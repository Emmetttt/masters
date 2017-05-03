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

##REFLEX II
scRA = [float(x)*degtorad for x in fun.scRA]    ##Southern Superclusters
scdec = [float(x)*degtorad for x in fun.scdec]
scz = [float(x) for x in fun.scz]
scRA = [x-2*np.pi if x > np.pi else x for x in scRA]


##Abell
ascRA = [float(x)*degtorad for x in fun.abelscRA]
ascdec = [float(x)*degtorad for x in fun.abelscdec]
ascz = fun.abelscz
ascRA = [x-2*np.pi if x > np.pi else x for x in ascRA]



fig = plt.figure()#figsize=(10, 5))
ax = fig.add_subplot(111, projection="mollweide")
ax.grid(True)

i=0
ax.plot(scRA, scdec, '.', color='blue', label='REFLEX II')
ax.plot(ascRA, ascdec, '.', color='grey', label='Abell')
plt.xlabel('Right Ascenscion')
plt.ylabel('Declination')
plt.legend()
plt.show()
