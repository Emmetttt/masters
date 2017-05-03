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

def taylorexpanded(z, H):
    c = fun.c
    DE = 0.7205
    M = 0.2795
    return (c/H)*(z + (0.5*(1 + DE - (M/2)))*z**2)

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

##Combine Clusters
clusterRA = scRA + ascRA
clusterdec = scdec + ascdec
clusterz = scz + ascz
clusterRadii = scRadius + ascRadius

scin = []
scout = []
i=0
k=0
while k < len(SNered): ##Find the nearest sc for each k Supernova
    while i < len(clusterz):
        distancesNew = fun.polar_dist_comp(SNeRA[k], clusterRA[i], SNedec[k], clusterdec[i], SNered[k], clusterz[i]) ##Keep SN constant, loop over supercluster
        if fun.D(fun.HubbleIntegrate(distancesNew)) < clusterRadii[i]:
            scin.append([SNedist[k], SNedisterror[k], SNered[k], clusterRadii[i], fun.D(fun.HubbleIntegrate(distancesNew)) - clusterRadii[i]])
            i = len(clusterz)+1
        if i == len(clusterz) - 1:
            scout.append([SNedist[k], SNedisterror[k], SNered[k]])
        i+=1
    i=0
    k+=1
    
##Calculate dH for each supernovae
dHin = []
i=0
while i < len(scin):
    dHin.append(deltaH(scin[i][2], scin[i][0]))
    i+=1

##Calculate dH for each supernovae
dHout = []
i=0
while i < len(scout):
    dHout.append(deltaH(scout[i][2], scout[i][0]))
    i+=1

dists = [fun.D(x[0]) for x in scin]
distserr = [fun.errorD(x[0],x[1]) for x in scin]
reds = [x[2] for x in scin]
redshiftrange = np.linspace(min(reds), max(reds), 10)
bestfit = [taylorexpanded(x,71.2) for x in redshiftrange]
avgfit = [taylorexpanded(x,69.79) for x in redshiftrange]

##Graph for redshift vs distance
plt.figure()
plt.errorbar(reds, dists, yerr=distserr, fmt='.', markersize=5, color="blue", label='SNe from superclusters')
plt.plot(redshiftrange, bestfit, '-', color="blue", label='Best fit H0 = 71.2')
plt.plot(redshiftrange, avgfit, '-', color="grey", label='Global Avg H0 = 69.79')

plt.legend(loc="upper left")
plt.xlabel('Redshift')
plt.ylabel('Distance (Mpc)')
plt.title('Hubble Diagram for supernovae originating in superclusters')
plt.axis([0, max(reds)*1.1, 0, max(dists)*1.1])
plt.legend()
plt.show()




###KS TEST
print(ks_2samp(dHin, dHout))

'''
####H values
Hin = [fun.HubbleLCDM(x[2], fun.D(x[0])) for x in SNein]
Hinmean = sum(Hin)/len(Hin)
Hinerror = fun.stdev(Hin, Hinmean)/np.sqrt(len(Hin))
print('Hubble Constant within distance to cluster: ', Hinmean, ' pm ', Hinerror)
Hout = [fun.HubbleLCDM(x[2], fun.D(x[0])) for x in SNeout]
Houtmean = sum(Hout)/len(Hout)
Houterror = fun.stdev(Hout, Houtmean)/np.sqrt(len(Hout))
print('Hubble Constant within distance to cluster: ', Houtmean, ' pm ', Houterror)
'''

####CUMFREQ
dHin.sort()
dHout.sort()

plt.plot(dHin, fun.cumfreq(dHin), label="Cluster Origin, (N=%i)"%len(dHin), color='blue')
plt.plot(dHout, fun.cumfreq(dHout), label="Not cluster origin (N=%i)"%len(dHout), color='green')
emptyNLSS = plt.hist([], range=[-1,1], alpha=0.5, label='Stat=%.3f \np-val=%.4f'%(ks_2samp(dHin, dHout)[0], ks_2samp(dHin, dHout)[1]), color='white')
plt.legend()
#plt.axis([-0.5, 0.5, 0, 1])
plt.xlabel('% Deviation from expected distance')
plt.ylabel('Cumulative frequency')
plt.title('Cumulative frequency graph for SN of cluster and non-cluster origin')
plt.show()



