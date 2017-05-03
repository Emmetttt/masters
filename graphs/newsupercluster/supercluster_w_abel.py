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


##Combine Clusters
clusterRA = scRA + ascRA
clusterdec = scdec + ascdec
clusterz = scz + ascz


distancescluster=len(SNered) * [[1000000]]
i=0
k=0
while k < len(SNered): ##Find the nearest sc for each k Supernova
    while i < len(clusterz):
        distancesNew = fun.polar_dist_comp(SNeRA[k], clusterRA[i], SNedec[k], clusterdec[i], SNered[k], clusterz[i]) ##Keep SN constant, loop over supercluster
        if distancesNew < distancescluster[k]:
            distancescluster[k] = fun.polar_dist_comp(SNeRA[k], clusterRA[i], SNedec[k], clusterdec[i], SNered[k], clusterz[i])
        i+=1
    i=0
    k+=1

distancesMpc = []
i=0
while i < len(distancescluster):
    distancesMpc.append(fun.D(fun.HubbleIntegrate(distancescluster[i])))
    i+=1

##Calculate dH for each supernovae
dH = []
i=0
while i < len(SNeRA):
    dH.append(deltaH(SNered[i], SNedist[i]))
    i+=1

plt.plot([np.log10(x) for x in distancesMpc], dH, '.')
#plt.plot([np.percentile(distances, bestpercentile), np.percentile(distances, bestpercentile)], [min(dH), max(dH)], '-.')
plt.xlabel('Logarithmic distance to nearest supercluster in Mpc')
plt.ylabel('% deviation luminsity distance')
plt.title('Deviation from H0 for supernovae as a function of\ndistance to nearest supercluster')
plt.legend()
plt.grid(True)
plt.show()

###KS TEST
'''
bestKS = 10
besti = 0
i = 3
while i < 1000:
    k=0
    cin = []
    cout = []
    while k < len(distancesMpc):
        if distancesMpc[k] < i:
            cin.append(dH[k])
        else:
            cout.append(dH[k])
        k+=1
    stat = ks_2samp(cin, cout)[1]
    if stat < bestKS:
        bestKS = stat
        besti = i
        print(besti, stat)
    i+=1
'''
cin = []
cout = []
SNein = []
SNeout = []
i=0
while i < len(distancesMpc):
    if distancesMpc[i] < 20:
        cin.append(dH[i])
        SNein.append([SNedist[i], SNedisterror[i], SNered[i]])
    else:
        cout.append(dH[i])
        SNeout.append([SNedist[i], SNedisterror[i], SNered[i]])
    i+=1

print(ks_2samp(cin, cout))

####H values
Hin = [fun.HubbleLCDM(x[2], fun.D(x[0])) for x in SNein]
Hinmean = sum(Hin)/len(Hin)
Hinerror = fun.stdev(Hin, Hinmean)/np.sqrt(len(Hin))
print('Hubble Constant within distance to cluster: ', Hinmean, ' pm ', Hinerror)
Hout = [fun.HubbleLCDM(x[2], fun.D(x[0])) for x in SNeout]
Houtmean = sum(Hout)/len(Hout)
Houterror = fun.stdev(Hout, Houtmean)/np.sqrt(len(Hout))
print('Hubble Constant within distance to cluster: ', Houtmean, ' pm ', Houterror)

####CUMFREQ
cin.sort()
cout.sort()

plt.plot(cin, fun.cumfreq(cin), label="Cluster Origin, (N=%i)"%len(cin), color='blue')
plt.plot(cout, fun.cumfreq(cout), label="Not cluster origin (N=%i)"%len(cout), color='green')
emptyNLSS = plt.hist([], range=[-1,1], alpha=0.5, label='Stat=%.3f \np-val=%.4f'%(ks_2samp(cin, cout)[0], ks_2samp(cin, cout)[1]), color='white')
plt.legend()
#plt.axis([-0.5, 0.5, 0, 1])
plt.xlabel('% Deviation from expected distance')
plt.ylabel('Cumulative frequency')
plt.title('Cumulative frequency graph for SN of cluster and non-cluster origin')
plt.show()



