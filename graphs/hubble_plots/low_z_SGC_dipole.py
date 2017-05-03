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
maxredshift = 0.1
minredshift = 0.01
error = 30 ##Size of annulus
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
##https://oup.silverchair-cdn.com/oup/backfile/Content_public/Journal/mnras/437/3/10.1093/mnras/stt2024/3/stt2024.pdf?Expires=1488468567&Signature=hK3NyZEPoBrbIbDt~MZLL4WKXzySHAkXnsLa0NA0BBSGH6P715QSYg~Ssff5H8q50dyztDgbl1RHBkzv49jKl5o0YgjJ9~5WFWKF2~M3O9akCsWLK4iwDMX9tJzAvD0xJbN7JkOaFn0jMzH4efoH6Vo9Gbe7eSorXatPe0ffbHobHDfnEXYQL-bEjru1J~vYQblOOeIi5SeZCG8QA9asKO6-eaWVwVaX8GnGXl0Fs~l0Yex8b21xDrvemEV3eejEvTPtcYDkgbA43uR-vLZ4pnwTcvbEyUnfqw-lpahf5~jQiANQT6kTAEb9zcGLsN2cthGDbK7PvjBk9R8oadZbOw__&Key-Pair-Id=APKAIUCZBIA4LVPAVW3Q
##Table 3
##The regions are located in the SGC (6dF SGC, 3511 deg2) at RA = 0 − 50 and 330 − 360) with DEC = −40 − 0 as well as at RA = 150 − 220o with DEC = −50 − 0 in the NGC (6dF NGC, 2578 deg2)
NGC = [192.85-360, 20.13]
#NGC = [220.85-360, -40.13]
SGC = [12.85, -20.13]
    
SN_NGC = []
SN_SGC = []
All = []
degtorad = 0.0174533

while i < len(redshifts):
    if redshifts[i] < maxredshift and redshifts[i] > minredshift:
        All.append([redshifts[i], dist[i], disterror[i], float(RA[i])*degtorad, float(dec[i])*degtorad, name[i]])
        if dec[i] > (SGC[1] - error) and dec[i] < (SGC[1] + error) and RA[i] > (SGC[0] - error) and RA[i] < (SGC[0] + error):
            SN_SGC.append([redshifts[i], dist[i], disterror[i], float(RA[i])*degtorad, float(dec[i])*degtorad, name[i]])
        elif dec[i] > (NGC[1] - error) and dec[i] < (NGC[1] + error):
            if RA[i] < (NGC[0] + error):
                    SN_NGC.append([redshifts[i], dist[i], disterror[i], float(RA[i])*degtorad, float(dec[i])*degtorad, name[i]])
            elif RA[i] > (NGC[0] - error + 360):
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
NGCH = [fun.c*x/68.41 for x in redshiftrange]
SGCH = [fun.c*x/66.9 for x in redshiftrange]

##Graph for redshift vs distance
plt.figure()
plt.errorbar(NGCred, NGCdist, yerr=NGCdisterror, fmt='.', markersize=5, color="blue", label='Dipole')
plt.errorbar(SGCred, SGCdist, yerr=SGCdisterror, fmt='.', markersize=5, color="red", label='SGC')
plt.plot(redshiftrange, SGCH, '-', color="red", label='H0 = 66.9 ± 1.6')
plt.plot(redshiftrange, NGCH, '-', color="blue", label='H0 = 68.4 ± 2.2')

plt.legend(loc="upper left")
plt.xlabel('Redshift')
plt.ylabel('Distance (Mpc)')
plt.title('Hubble Diagrams for SGC and its corresponding dipole')
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

plt.plot(SGCred, SGC_dH, '.')
plt.plot([min(SGCred), max(SGCred)], [0,0])
plt.show()

NGChist = plt.hist(NGC_dH, bins=binno, alpha=0.5, label='NGC', color='blue')
##Best fit for data
(NGCmu, NGCsigma) = norm.fit(NGC_dH)
emptyNGC = plt.hist([], range=[-1,1], alpha=0.5, label='$N=%i, \mu=%.2f$, $\sigma=%.2f$'%(len(NGC_dH), NGCmu, NGCsigma), color='white')

##Plot of the histogram of the NLSS

SGChist = plt.hist(SGC_dH, bins=binno, alpha=0.5, label='SGC', color='red')
(SGCmu, SGCsigma) = norm.fit(SGC_dH)
emptySGC = plt.hist([], range=[-1,1], alpha=0.5, label='$N=%i, \mu=%.2f$, $\sigma=%.2f$'%(len(SGC_dH), SGCmu, SGCsigma), color='white')
plt.xlabel('% Deviation from expected distance')
plt.ylabel('Frequency')
plt.title('Histogram of deviation from expected distance for SGC and its corresponding dipole supernovae')
plt.legend()
plt.show()


NGC_dH.sort()
SGC_dH.sort()

plt.plot(NGC_dH, fun.cumfreq(NGC_dH), label="NGC", color='blue')
plt.plot(SGC_dH, fun.cumfreq(SGC_dH), label="SGC", color='red')
emptyNLSS = plt.hist([], range=[-0.3,0.3], alpha=0.5, label='Stat=%.3f \np-val=%.4f'%(ks_2samp(NGC_dH, SGC_dH)[0], ks_2samp(NGC_dH, SGC_dH)[1]), color='white')
plt.xlabel('% Deviation from expected distance')
plt.ylabel('Frequency')
plt.title('Cumulative frequency for the deviations in the SGC and its corresponding dipole.')
plt.legend()
plt.show()
        





#########EQUITORIAL COORDINATE

fig = plt.figure()#figsize=(10, 5))
ax = fig.add_subplot(111, projection="mollweide")
ax.grid(True)

SGCboxx = np.array([SGC[0]-error, SGC[0]-error, SGC[0]+error, SGC[0]+error, SGC[0]-error])*degtorad
SGCboxy = np.array([SGC[1]-error, SGC[1]+error, SGC[1]+error, SGC[1]-error, SGC[1]-error])*degtorad

NGCboxxwest = np.array([-180, NGC[0]+error, NGC[0]+error, -180])*degtorad
NGCboxywest = np.array([NGC[1]-error, NGC[1]-error, NGC[1]+error, NGC[1]+error])*degtorad
NGCboxxeast = np.array([180, NGC[0]-error+360, NGC[0]-error+360, 180])*degtorad
NGCboxyeast = np.array([NGC[1]-error, NGC[1]-error, NGC[1]+error, NGC[1]+error])*degtorad


ax.plot(AllRA, Alldec, '.', color='black', label='All SNe for z<0.1')
ax.plot(NGCRA, NGCdec, '.', color='blue', label='NGC')
ax.plot(NGCboxxwest, NGCboxywest, '-', color='blue')
ax.plot(NGCboxxeast, NGCboxyeast, '-', color='blue')
ax.plot(SGCRA, SGCdec, '.', color='red', label='SGC')
ax.plot(SGCboxx, SGCboxy, '-', color='red')
plt.xlabel('RA')
plt.ylabel('DEC')
plt.title('Positions for z<0.1 supernovae in equatorial coordinates')
plt.legend()

if ShowPosGraph:
        plt.show()




#####GALACTIC COORDINATE TRANSFORMATION#######

NgalRA = []
Ngaldec = []
i = 0
while i < len(NGCRA):
        x = SkyCoord(ra=NGCRA[i]*u.radian, dec=NGCdec[i]*u.radian, frame='icrs')
        l = x.galactic.l.radian
        if l > np.pi:
                l = l-2*np.pi
        NgalRA.append(l)
        Ngaldec.append(x.galactic.b.radian)
        i+=1

SgalRA = []
Sgaldec = []
i = 0
while i < len(SGCRA):
        x = SkyCoord(ra=SGCRA[i]*u.radian, dec=SGCdec[i]*u.radian, frame='icrs')
        l = x.galactic.l.radian
        SgalRA.append(l)
        Sgaldec.append(x.galactic.b.radian)
        i+=1



fig = plt.figure()#figsize=(10, 5))
ax = fig.add_subplot(111, projection="mollweide")
ax.grid(True)

ax.plot(NgalRA, Ngaldec, '.', color='blue', label='NGC')
ax.plot(SgalRA, Sgaldec, '.', color='red', label='SGC')
plt.xlabel('RA')
plt.ylabel('DEC')
plt.title('Galactic Coordinates')
plt.legend()

if ShowPosGraph:
        plt.show()

