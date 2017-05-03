##import graphing
from astropy import units as u
from astropy.coordinates import SkyCoord
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


###INPUT###
maxredshift = 0.2
minredshift = 0
error = 50 ##Size of annulus
binno = 7

ShowPosGraph = True

##Data
RA = fun.RA
dec = fun.dec
dist = fun.distMpc
disterror = fun.disterrMpc
redshifts = fun.redshifts
name = fun.name

j=0
while j < len(RA):
        if RA[j] > 180:
                RA[j] = RA[j] - 360
        j+=1
   

i = 0
NGC = [-30, -62] ##MAX
SGC = [150, 62]  ##MIN
    
SN_NGC = []
SN_SGC = []
All = []
degtorad = 0.0174533
'''
while i < len(redshifts):
    if redshifts[i] < maxredshift:
        All.append([redshifts[i], dist[i], disterror[i], float(RA[i])*degtorad, float(dec[i])*degtorad, name[i]])
        if dec[i] > (SGC[1] - error) and RA[i] > (SGC[0] - error) and dec[i] < (SGC[1] + error) and RA[i] < (SGC[0] + error):
            SN_SGC.append([redshifts[i], dist[i], disterror[i], float(RA[i])*degtorad, float(dec[i])*degtorad, name[i]])
        elif dec[i] < (SGC[1] + error - 180) and RA[i] < (SGC[0] + error - 360):
            SN_SGC.append([redshifts[i], dist[i], disterror[i], float(RA[i])*degtorad, float(dec[i])*degtorad, name[i]])
        elif dec[i] < (NGC[1] + error) and RA[i] > (NGC[0] - error) and RA[i] < (NGC[0] + error) and dec[i] > (NGC[1] - error):
            SN_NGC.append([redshifts[i], dist[i], disterror[i], float(RA[i])*degtorad, float(dec[i])*degtorad, name[i]])
        elif RA[i] > (NGC[0] - error) and RA[i] < (NGC[0] + error) and dec[i] > (NGC[1] - error + 180):
            SN_NGC.append([redshifts[i], dist[i], disterror[i], float(RA[i])*degtorad, float(dec[i])*degtorad, name[i]])
    i+=1
'''


while i < len(redshifts):
    if redshifts[i] < maxredshift:
        All.append([redshifts[i], dist[i], disterror[i], float(RA[i])*degtorad, float(dec[i])*degtorad, name[i]])
        if (dec[i] > (SGC[1] - error) and dec[i] < (SGC[1] + error)) or dec[i] < (SGC[1] + error - 180):
                if (RA[i] > (SGC[0] - error) and RA[i] < (SGC[0] + error)) or RA[i] < (SGC[0] + error - 360):
                        SN_SGC.append([redshifts[i], dist[i], disterror[i], float(RA[i])*degtorad, float(dec[i])*degtorad, name[i]])
        if (dec[i] < (NGC[1] + error) and dec[i] > (NGC[1] - error)) or dec[i] > (NGC[1] - error + 180):
                if (RA[i] > (NGC[0] - error) and RA[i] < (NGC[0] + error)):
                        SN_NGC.append([redshifts[i], dist[i], disterror[i], float(RA[i])*degtorad, float(dec[i])*degtorad, name[i]])
    i+=1





##Plots
NGCred = np.array([x[0] for x in SN_NGC])
NGCdisterror = np.array([x[2] for x in SN_NGC])
NGCdist = np.array([x[1] for x in SN_NGC])
NGCRA = [x[3] for x in SN_NGC]
NGCdec = [x[4] for x in SN_NGC]
NGC_H = fun.HMpc(np.array(NGCred), np.array(NGCdist))
NGC_dH = fun.dHMpc(np.array(NGCred), np.array(NGCdist), np.array(NGCdisterror))

SGCred = np.array([x[0] for x in SN_SGC])
SGCdisterror = np.array([x[2] for x in SN_SGC])
SGCdist = np.array([x[1] for x in SN_SGC])
SGCRA = [x[3] for x in SN_SGC]
SGCdec = [x[4] for x in SN_SGC]
SGC_H = fun.HMpc(np.array(SGCred), np.array(SGCdist))
SGC_dH = fun.dHMpc(np.array(SGCred), np.array(SGCdist), np.array(SGCdisterror))

Allred = [x[0] for x in All]
Alldist = [x[1] for x in All]
Alldisterror = [x[2] for x in All]
AllRA = [x[3] for x in All]
Alldec = [x[4] for x in All]

minres = 100000
NGCH = 0
i = 0
while i < 1000:
        h = 65 + i/100
        res = fun.residual(fun.dmod(NGCdist), fun.dmod(NGCdisterror), NGCred, h, 0.2807, 0.7193)
        if res < minres:
                minres = res
                NGCH = h
        i+=1

redshiftrange = np.linspace(minredshift, maxredshift, 10)
NGCH = [fun.c*x/67.9 for x in redshiftrange]
SGCH = [fun.c*x/68.8 for x in redshiftrange]

##Graph for redshift vs distance
plt.figure()
plt.errorbar(NGCred, NGCdist, yerr=NGCdisterror, fmt='.', markersize=5, color="blue", label='Maximal expansion')
plt.errorbar(SGCred, SGCdist, yerr=SGCdisterror, fmt='.', markersize=5, color="red", label='Minimal expansion')
plt.plot(redshiftrange, SGCH, '-', color="red", label='H0 = 68.8 ± 1.5')
plt.plot(redshiftrange, NGCH, '-', color="blue", label='H0 = 67.9 ± 4.0')

plt.legend(loc="upper left")
plt.xlabel('Redshift')
plt.ylabel('Distance (Mpc)')
plt.title('Hubble Diagrams for maximal and minimal expansion rates')
plt.axis([0, max(SGCred)*1.1, 0, max(SGCdist)*1.1])
plt.legend()
plt.show()

##Graph for redshfit vs H0
plt.figure(2)
plt.errorbar(NGCred, NGC_H, yerr=NGC_dH, fmt='.', markersize=5, color="blue", label='Dipole')
plt.errorbar(SGCred, SGC_H, yerr=SGC_dH, fmt='.', markersize=5, color="red", label='SGC')
plt.plot(redshiftrange, 10*[68.4], '-', color="blue", label='H0 = 68.4 ± 2.2')
plt.plot(redshiftrange, 10*[68.4+2.2], ':', color="blue", alpha=0.3)
plt.plot(redshiftrange, 10*[68.4-2.2], ':', color="blue", alpha=0.3)
plt.plot(redshiftrange, 10*[66.9], '-', color="red", label='H0 = 66.9 ± 1.6')
plt.plot(redshiftrange, 10*[66.9+1.6], ':', color="red", alpha=0.3)
plt.plot(redshiftrange, 10*[66.9-1.6], ':', color="red", alpha=0.3)



plt.legend(loc="upper left")
plt.xlabel('Redshift')
plt.ylabel('Hubble constant (km s$^{-1}$ Mpc$^{-1}$)')
plt.title('Hubble values for SGC and its corresponding dipole')
plt.legend()
plt.show()




##Histogram of dH
NGC_dH = []
SGC_dH = []

i = 0
while i < len(NGCred):
        NGC_dH.append(( fun.D(fun.HubbleIntegrate(NGCred[i])) - NGCdist[i])/fun.D(fun.HubbleIntegrate(NGCred[i])))
        i+=1
i = 0
while i < len(SGCred):
        SGC_dH.append((fun.D(fun.HubbleIntegrate(SGCred[i])) - SGCdist[i])/fun.D(fun.HubbleIntegrate(SGCred[i])))
        i+=1

print(ks_2samp(NGC_dH, SGC_dH))

NGC_dH.sort()
SGC_dH.sort()

plt.plot(NGC_dH, fun.cumfreq(NGC_dH), label="Maximal expansion", color='blue')
plt.plot(SGC_dH, fun.cumfreq(SGC_dH), label="Minimal expansion", color='red')
emptyNLSS = plt.hist([], range=[-0.3,0.3], alpha=0.5, label='Stat=%.3f \np-val=%.4f'%(ks_2samp(NGC_dH, SGC_dH)[0], ks_2samp(NGC_dH, SGC_dH)[1]), color='white')
plt.xlabel('$\delta d_L$')
plt.ylabel('Cumulative Frequency')
plt.title('Cumulative frequency for the deviations in the\n maximal and minimal rate of expansion')
plt.legend()
plt.show()
        





#########EQUITORIAL COORDINATE

fig = plt.figure()#figsize=(10, 5))
ax = fig.add_subplot(111, projection="mollweide")
ax.grid(True)

SGCboxxe = np.array([SGC[0]-error, SGC[0]-error, 180, 180, SGC[0]-error])*degtorad
SGCboxye = np.array([SGC[1]-error, 90, 90, SGC[1]-error, SGC[1]-error])*degtorad

SGCboxxw = np.array([-180, SGC[0]+error-360, SGC[0]+error-360, -180, -180])*degtorad
SGCboxyw = np.array([SGC[1]-error, SGC[1]-error, 90, 90, SGC[1]-error])*degtorad

NGCboxx = np.array([NGC[0]-error, NGC[0]-error, NGC[0]+error, NGC[0]+error])*degtorad
NGCboxy = np.array([NGC[1]-error, NGC[1]+error, NGC[1]+error, NGC[1]-error])*degtorad



ax.plot(AllRA, Alldec, '.', color='black', label='All SNe for z<0.2')
ax.plot(NGCRA, NGCdec, '.', color='blue', label='Maximal expansion')
ax.plot(NGCboxx, NGCboxy, '-', color='blue')
ax.plot(SGCRA, SGCdec, '.', color='red', label='Minimal expansion')
ax.plot(SGCboxxe, SGCboxye, '-', color='red')
ax.plot(SGCboxxw, SGCboxyw, '-', color='red')
plt.xlabel('RA')
plt.ylabel('DEC')
plt.title('Positions of dipole')
plt.legend()

if ShowPosGraph:
        plt.show()


