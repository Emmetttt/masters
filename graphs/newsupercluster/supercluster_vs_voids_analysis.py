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

def bestfit(x, a):
    return ((3*10**5)/a)*(x + (0.5*(1 + 0.7205 - (0.2795/2)))*x**2 + (1/6)*(1+0.7205 - (0.2795/2) - (3*((0.2795/2)-0.7205)**2+1))*x**3)


##SUPERNOVA##
SNeRA = fun.RA
SNedec = fun.dec
SNedist = fun.dist
SNedisterror = fun.disterror
SNedistMpc = fun.distMpc
SNered = fun.redshifts
degtorad = 0.0174533

SNeSC = eval(open("C:\\Python34\\masters\\graphs\\newsupercluster\\SNSC.txt", "r").read())
SNevoid = eval(open("C:\\Python34\\masters\\graphs\\newsupercluster\\SNvoid.txt", "r").read())

SCdist = [fun.D(x[0]) for x in SNeSC]
SCdisterror = [fun.errorD(x[0], x[1]) for x in SNeSC]
SCred = [x[2] for x in SNeSC]

Voiddist = [fun.D(x[0]) for x in SNevoid]
Voiddisterror = [fun.errorD(x[0], x[1]) for x in SNevoid]
Voidred = [x[2] for x in SNevoid]


##Calculate dH for each supernovae
dHin = []
i=0
while i < len(SCdist):
    dHin.append(deltaH(SCred[i], SCdist[i]))
    i+=1

##Calculate dH for each supernovae
dHout = []
i=0
while i < len(Voiddist):
    dHout.append(deltaH(Voidred[i], Voiddist[i]))
    i+=1



redshiftrange = np.linspace(0, max(Voidred+SCred), 10)
SCfit = [bestfit(x,71.2) for x in redshiftrange]
Voidfit = [bestfit(x,69.79) for x in redshiftrange]

##Graph for redshift vs distance
plt.figure()
plt.errorbar(SCred, SCdist, yerr=SCdisterror, fmt='.', markersize=5, color="blue", label='SNe from superclusters')
plt.errorbar(Voidred, Voiddist, yerr=Voiddisterror, fmt='.', markersize=5, color="red", label='SNe from Voids')
plt.plot(redshiftrange, SCfit, '-', color="blue", label='Supercluster H0 = 71.45 ± 1.67')
plt.plot(redshiftrange, Voidfit, '-', color="red", label='Void H0 = 71.84 ± 1.21')

plt.legend(loc="upper left")
plt.xlabel('Redshift')
plt.ylabel('Distance (Mpc)')
plt.title('Hubble Diagram for supernovae originating in superclusters')
plt.axis([0, max(Voidred+SCred)*1.1, 0, max(Voiddist+SCdist)*1.1])
plt.legend()
plt.show()





###KS TEST
print(ks_2samp(dHin, dHout))

####CUMFREQ
dHin.sort()
dHout.sort()

plt.plot(dHin, fun.cumfreq(dHin), label="Cluster Origin, (N=%i)"%len(dHin), color='blue')
plt.plot(dHout, fun.cumfreq(dHout), label="Void origin (N=%i)"%len(dHout), color='green')
emptyNLSS = plt.hist([], range=[-1,1], alpha=0.5, label='Stat=%.3f \np-val=%.4f'%(ks_2samp(dHin, dHout)[0], ks_2samp(dHin, dHout)[1]), color='white')
plt.legend()
#plt.axis([-0.5, 0.5, 0, 1])
plt.xlabel('$\delta d_L$')
plt.ylabel('Cumulative frequency')
plt.title('Cumulative frequency graph for SN of cluster and non-cluster origin')
plt.show()


