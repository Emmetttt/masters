##import graphing
import matplotlib.pyplot as plt
import numpy as np
##import functions
import sys
sys.path.insert(0, 'C:\Python34\masters\graphs')
import functions as fun

redshifts = fun.redshifts
dist = fun.distMpc
disterror = fun.disterror
RA = fun.RA
dec = fun.dec



###INPUT###
maxredshift = 0.05
minredshift = 0.015
error = 40 ##Size of annulus


i = 0
SGC = [12.85, -27.13]
    
SN_SGC = []
All = []

while i < len(redshifts):
    if RA[i] > (SGC[0] - error) and RA[i] < (SGC[0] + error) and dec[i] > (SGC[1] - error) and dec[i] < (SGC[1] + error) and redshifts[i] > minredshift and redshifts[i] < maxredshift:
        SN_SGC.append([redshifts[i], dist[i], disterror[i]])
        All.append([redshifts[i], dist[i], disterror[i]])
        i+=1
    else:
        i+=1

##Plots
del SN_SGC[-6] ##obvious outlier
del SN_SGC[4] ##obvious outlier
SGCy = np.array([x[0] for x in SN_SGC])
SGCxerror = np.array([x[2] for x in SN_SGC])
SGCx = np.array([x[1] for x in SN_SGC])



##Calculate polynomial for SGC
SGCx = SGCx[:,np.newaxis]
#a, _, _, _ = np.linalg.lstsq(SGCx, SGCy)

SGCx_new = np.linspace(min(SGCx), max(SGCx), 50)

i = 0
a = 60
minchi = 100
bestgrad = 0
while i < 10000:
    a = a + (0.001 * i)
    SGCy_new = (a*SGCx_new)/(3*10**5)
    if fun.ChiSq(SGCy, SGCx, a) < minchi:
        minchi = fun.ChiSq(SGCy, SGCx, a)
        bestgrad = a
    i+=1

SGCy_new = (bestgrad*SGCx_new)/(3*10**5)

print("SGC Hubble Constant Found: ", bestgrad, " kms^-1Mpc^-1")
print("Chi Squared value of: ", fun.ChiSq(SGCy, SGCx, bestgrad))


##Hubble
z = np.arange(0, max(SGCy)+0.01, 0.005)
expD = [z/0.00023333 for z in z] ##0.0002333 is 70 / 3*10^5

##Graph for redshift vs distance
plt.figure()
plt.plot(expD, z, color="red", label='H')
plt.plot(SGCx_new, SGCy_new, color="blue", label='SGC') ##Best fit for SGC
plt.errorbar(SGCx, SGCy, xerr=SGCxerror, fmt='x', markersize=5, color="blue")

plt.legend(loc="upper left")
plt.ylabel('Redshift')
plt.xlabel('Distance (Mpc)')
plt.title('Hubble Diagrams for NGC and SGC')
plt.axis([0, max(expD), 0, max(SGCy)+0.01])
plt.show()
