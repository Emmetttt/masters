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

def deltaH(SNered, SNedist):
    return (fun.HubbleIntegrate(float(SNered)) - (float(SNedist)))  / fun.HubbleIntegrate(float(SNered))*100

##SUPERNOVA##
SNeRA = fun.RA
SNedec = fun.dec
SNedist = fun.dist
SNedisterror = fun.disterror
SNedistMpc = fun.distMpc
SNered = fun.redshifts
degtorad = 0.0174533

##SUPERCLUSTER##
scRA = [float(x) for x in fun.scRA]    ##Southern Superclusters
scdec = [float(x) for x in fun.scdec]
scz = [float(x) for x in fun.scz]
scName = fun.scName
scRadius = [float(x) for x in fun.R]

##ABEL SUPERCLUSTER##
ascRA = fun.abelscRA
ascdec = fun.abelscdec
ascz = fun.abelscz
ascRadius = fun.abelscradii

##SDSS SUPERCLUSTER##
SDSSRA = fun.SDSSRA
SDSSdec = fun.SDSSdec
SDSSz = fun.SDSSz
SDSSradius = fun.SDSSradius

##Combine Clusters
clusterRA = SDSSRA + scRA + ascRA
clusterdec = SDSSdec + scdec + ascdec
clusterz = SDSSz + scz + ascz
clusterRadii = SDSSradius + scRadius + ascRadius

##Voids
voidRA = fun.minvoidRA + fun.isovoidRA
voiddec = fun.minvoiddec + fun.isovoiddec
voidz = fun.minvoidz + fun.isovoidz
voidRadii = fun.minvoidradius + fun.isovoidradius

distancescluster=len(SNered) * [1000000]
distancevoid = len(SNered) * [1000000]
SNSC = len(SNered) * [0]
SNVoid = len(SNered) * [0]
i=0
j=0
k=0
while k < len(SNered): ##Find the nearest sc for each k Supernova
    while i < len(clusterz):
        distancesNew = fun.D(fun.HubbleIntegrate(fun.polar_dist_comp(SNeRA[k], clusterRA[i], SNedec[k], clusterdec[i], SNered[k], clusterz[i]))) ##Keep SN constant, loop over supercluster
        if distancesNew < clusterRadii[i]:
            SNSC[k] = 1
        if distancesNew < distancescluster[k]:
            distancescluster[k] = distancesNew
        i+=1
    while j < len(voidRA):
        distancesNew = fun.D(fun.HubbleIntegrate(fun.polar_dist_comp(SNeRA[k], voidRA[j], SNedec[k], voiddec[j], SNered[k], voidz[j]))) ##Keep SN constant, loop over supercluster
        if distancesNew < voidRA[j]:
            SNVoid[k] = 1
        if distancesNew < distancevoid[k]:
            distancevoid[k] = distancesNew
        j+=1
    j=0
    i=0
    k+=1
    
##Calculate dH for each supernovae
dH = []
i=0
while i < len(SNeRA):
    dH.append(deltaH(SNered[i], fun.D(SNedist[i])))
    i+=1


plt.figure(1)
plt.subplot(211)
plt.plot([np.log10(x) for x in distancescluster], dH, '.')

plt.subplot(212)
plt.plot([np.log10(x) for x in distancevoid], dH, '.')

plt.xlabel('Logarithmic distance to nearest supercluster in Mpc')
plt.ylabel('$\delta d_L$')
plt.legend()
plt.grid(True)
plt.show()
