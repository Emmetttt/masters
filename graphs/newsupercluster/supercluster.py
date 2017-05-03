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
    return (fun.D(fun.HubbleIntegrate(float(SNered))) - (float(SNedist)))  / fun.D(fun.HubbleIntegrate(float(SNered)))



##USER INPUT
binno = 25

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


confirmedSC = [] ##[0] is SNe number, [1] is corresponding SC
unconfirmedSC = []

i=0
z=0
k=0
while k < len(SNered): ##Find the nearest sc for each k Supernova
    if SNered[k] < max(scz)+0.02 and SNedec[k] < max(scdec)+5:
        while i < len(scz):
            distancesNew = fun.polar_dist_comp(SNeRA[k], scRA[i], SNedec[k], scdec[i], SNered[k], scz[i]) ##Keep SN constant, loop over supercluster
            if fun.D(fun.HubbleIntegrate(distancesNew)) < scRadius[i]:
                print(SNeRA[k], scRA[i], SNedec[k], scdec[i], SNered[k], scz[i],i,k)
                confirmedSC.append(k)
                i=len(scz)
            if i == len(scz)-1:
                unconfirmedSC.append(k)
            i+=1
        i=0
        k+=1
    else:
        if SNered[k] < max(scz)+0.02 and SNedec[k] > max(scdec)+5:
            z+=1
        del SNered[k]
        del SNeRA[k]
        del SNedec[k]
        del SNedist[k]
        del SNedisterror[k]
        del SNedistMpc[k]

##Loop to determine the deviation from hubble constant
##for each Supernova
within_dist_dH = []
outside_dist_dH = []
m=0

while m < len(SNered):
    if m in confirmedSC:
        within_dist_dH.append(deltaH(SNered[m], SNedistMpc[m]))
    elif m in unconfirmedSC:
        outside_dist_dH.append(deltaH(SNered[m], SNedistMpc[m]))
    m+=1


print(ks_2samp(within_dist_dH, outside_dist_dH))

binwidth = 0.07

##Histograms
plt.hist(within_dist_dH, bins=np.arange(-1, 1, binwidth), alpha=0.5, label='Cluster origin', color='blue')
##Best fit for data
(LSSmu, LSSsigma) = norm.fit(within_dist_dH)
emptyLSS = plt.hist([], range=[-1,1], alpha=0.5, label='$N=%i, \mu=%.2f$, $\sigma=%.2f$'%(len(within_dist_dH), LSSmu, LSSsigma), color='white')

##Plot of the histogram of the NLSS

NLSShist = plt.hist(outside_dist_dH, bins=np.arange(-1, 1, binwidth), alpha=0.5, label='Not cluster origin', color='green')
(NLSSmu, NLSSsigma) = norm.fit(outside_dist_dH)
emptyNLSS = plt.hist([], range=[-1,1], alpha=0.5, label='$N=%i, \mu=%.2f$, $\sigma=%.2f$'%(len(outside_dist_dH), NLSSmu, NLSSsigma), color='white')
plt.xlabel('% Deviation from expected distance')
plt.ylabel('Frequency')
plt.title('Histogram of luminosity distance diviation for supernovae of \ncluster origin and not of cluster origin')
plt.legend()
plt.show()

within_dist_dH.sort()
outside_dist_dH.sort()

plt.plot(within_dist_dH, fun.cumfreq(within_dist_dH), label="Cluster Origin, (N=%i)"%len(within_dist_dH), color='blue')
plt.plot(outside_dist_dH, fun.cumfreq(outside_dist_dH), label="Not cluster origin (N=%i)"%len(outside_dist_dH), color='green')
emptyNLSS = plt.hist([], range=[-1,1], alpha=0.5, label='Stat=%.3f \np-val=%.4f'%(ks_2samp(within_dist_dH, outside_dist_dH)[0], ks_2samp(within_dist_dH, outside_dist_dH)[1]), color='white')
plt.legend()
#plt.axis([-0.5, 0.5, 0, 1])
plt.xlabel('% Deviation from expected distance')
plt.ylabel('Cumulative frequency')
plt.title('Cumulative frequency graph for SN of cluster and non-cluster origin')
plt.show()
        







