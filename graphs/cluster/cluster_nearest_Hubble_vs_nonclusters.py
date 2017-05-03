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

distancesMpc = [eval(x) for x in open("C:\\Python34\\masters\\graphs\\newcluster\\SDSS_planck.txt", "r").read().split(", ")]


##Calculate dH for each supernovae
dH = []
i=0
while i < len(SNeRA):
    dH.append(deltaH(SNered[i], SNedist[i]))
    i+=1

    
clusterIN = []
clusterOUT = []

i = 0
while i < len(SNeRA):
    if i in [566, 567, 569, 571, 576]:
        clusterIN.append([SNered[i], SNedist[i], SNedisterror[i], dH[i]])
    elif i in [564, 565, 568, 570, 572, 573, 574, 575, 577]:
        clusterOUT.append([SNered[i], SNedist[i], SNedisterror[i], dH[i]])
    elif distancesMpc[i] < 5:
        clusterIN.append([SNered[i], SNedist[i], SNedisterror[i], dH[i]])
    i+=1

cin = [x[3] for x in clusterIN]
cout = [x[3] for x in clusterOUT]

print(ks_2samp(cin, cout))

'''

redshiftrange = np.linspace(minredshift, maxredshift, 10)
INH = [fun.c*x/largebestH for x in redshiftrange]
OUTH = [fun.c*x/smallbestH for x in redshiftrange]

##Graph for redshift vs distance
plt.figure()
plt.errorbar([x[0] for x in clusterIN], [fun.D(x[1]) for x in clusterIN], yerr=[fun.errorD(x[1],x[2]) for x in clusterIN], fmt='.', markersize=5, color="blue", label='SNe risiding in clusters')
plt.errorbar([x[0] for x in clusterIN], [fun.D(x[1]) for x in clusterIN], yerr=[fun.errorD(x[1],x[2]) for x in clusterIN], fmt='.', markersize=5, color="red", label='SNe outside of clusters')
plt.plot(redshiftrange, smallH, '-', color="red", label='H0 = 73.7 ± 4')
plt.plot(redshiftrange, largeH, '-', color="blue", label='H0 = 64.1 ± 3.5')

plt.legend(loc="upper left")
plt.xlabel('Redshift')
plt.ylabel('Distance (Mpc)')
plt.title('Hubble Diagrams for the maximal and minimal H0')
plt.axis([min(np.append(largered, smallred))*0.9, max(np.append(largered, smallred))*1.1, min(np.append(largedist, smalldist))*0.9, max(np.append(largedist, smalldist))*1.1])
plt.legend()
plt.show()
'''

####CUMFREQ
cin.sort()
cout.sort()

plt.plot(cin, fun.cumfreq(cin), label="Cluster Origin, (N=%i)"%len(cin), color='blue')
plt.plot(cout, fun.cumfreq(cout), label="Not cluster origin (N=%i)"%len(cout), color='green')
emptyNLSS = plt.hist([], range=[-1,1], alpha=0.5, label='Stat=%.3f \np-val=%.4f'%(ks_2samp([x[3] for x in clusterIN], [x[3] for x in clusterOUT])[0], ks_2samp([x[3] for x in clusterIN], [x[3] for x in clusterOUT])[1]), color='white')
plt.legend()
#plt.axis([-0.5, 0.5, 0, 1])
plt.xlabel('% Deviation from expected distance')
plt.ylabel('Cumulative frequency')
plt.title('Cumulative frequency graph for SN of cluster and non-cluster origin')
plt.show()



