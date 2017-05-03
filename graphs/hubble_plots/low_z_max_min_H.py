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
##The regions are located in the small (6dF small, 3511 deg2) at RA = 0 − 50 and 330 − 360) with DEC = −40 − 0 as well as at RA = 150 − 220o with DEC = −50 − 0 in the large (6dF large, 2578 deg2)
large = [78.14028977641314  , 30.554766905652883]
small = [5.870396949370161  , 54.1304526638786]
    
SN_large = []
SN_small = []
All = []
degtorad = 0.0174533

while i < len(redshifts):
    if redshifts[i] < maxredshift and redshifts[i] > minredshift:
        All.append([redshifts[i], dist[i], disterror[i], float(RA[i])*degtorad, float(dec[i])*degtorad, name[i]])
        if dec[i] > (small[1] - error) and dec[i] < (small[1] + error) and RA[i] > (small[0] - error) and RA[i] < (small[0] + error):
            SN_small.append([redshifts[i], dist[i], disterror[i], float(RA[i])*degtorad, float(dec[i])*degtorad, name[i]])
        if dec[i] > (large[1] - error) and dec[i] < (large[1] + error) and RA[i] > (large[0] - error) and RA[i] < (large[0] + error):
            SN_large.append([redshifts[i], dist[i], disterror[i], float(RA[i])*degtorad, float(dec[i])*degtorad, name[i]])
    i+=1


##Plots
largered = np.array([x[0] for x in SN_large])
largedisterror = np.array([x[2] for x in SN_large])
largedist = np.array([x[1] for x in SN_large])
largeRA = [x[3] for x in SN_large]
largedec = [x[4] for x in SN_large]
large_H = fun.HMpc(np.array(largered), np.array(largedist))
large_dH = fun.dHMpc(np.array(largered), np.array(largedist), np.array(largedisterror))

smallred = np.array([x[0] for x in SN_small])
smalldisterror = np.array([x[2] for x in SN_small])
smalldist = np.array([x[1] for x in SN_small])
smallRA = [x[3] for x in SN_small]
smalldec = [x[4] for x in SN_small]
small_H = fun.HMpc(np.array(smallred), np.array(smalldist))
small_dH = fun.dHMpc(np.array(smallred), np.array(smalldist), np.array(smalldisterror))


Allred = [x[0] for x in All]
Alldist = [x[1] for x in All]
Alldisterror = [x[2] for x in All]
AllRA = [x[3] for x in All]
Alldec = [x[4] for x in All]


minres = 100000
largeH = 0
i = 0
while i < 1000:
        h = 65 + i/100
        res = fun.residual(fun.dmod(largedist), fun.dmod(largedisterror), largered, h, 0.2807, 0.7193)
        if res < minres:
                minres = res
                largeH = h
        i+=1

        
largebestH = 73.7
smallbestH = 64.1

redshiftrange = np.linspace(minredshift, maxredshift, 10)
largeH = [fun.c*x/largebestH for x in redshiftrange]
smallH = [fun.c*x/smallbestH for x in redshiftrange]



##Graph for redshift vs distance
plt.figure()
plt.errorbar(largered, largedist, yerr=largedisterror, fmt='.', markersize=5, color="blue", label='Maximal H0')
plt.errorbar(smallred, smalldist, yerr=smalldisterror, fmt='.', markersize=5, color="red", label='Minimal H0')
plt.plot(redshiftrange, smallH, '-', color="red", label='H0 = 73.7 ± 4.5')
plt.plot(redshiftrange, largeH, '-', color="blue", label='H0 = 64.1 ± 3.5')

plt.legend(loc="upper left")
plt.xlabel('Redshift')
plt.ylabel('Distance (Mpc)')
plt.title('Hubble Diagrams for the maximal and minimal H0')
plt.axis([min(np.append(largered, smallred))*0.9, max(np.append(largered, smallred))*1.1, min(np.append(largedist, smalldist))*0.9, max(np.append(largedist, smalldist))*1.1])
plt.legend()
plt.show()



##Graph for redshfit vs H0
plt.figure(2)
plt.errorbar(largered, large_H, yerr=large_dH, fmt='.', markersize=5, color="blue", label='Maximal H0')
plt.errorbar(smallred, small_H, yerr=small_dH, fmt='.', markersize=5, color="red", label='Minimal H0')
plt.plot(redshiftrange, 10*[73.7], '-', color="blue", label='H0 = 73.7 ± 4.5')
plt.plot(redshiftrange, 10*[73.7+4.5], ':', color="blue", alpha=0.3)
plt.plot(redshiftrange, 10*[73.7-4.5], ':', color="blue", alpha=0.3)
plt.plot(redshiftrange, 10*[64.1], '-', color="red", label='H0 = 64.1 ± 3.5')
plt.plot(redshiftrange, 10*[64.1+3.5], ':', color="red", alpha=0.3)
plt.plot(redshiftrange, 10*[64.1-3.5], ':', color="red", alpha=0.3)



plt.legend(loc="upper left")
plt.xlabel('Redshift')
plt.ylabel('Hubble constant (km s$^{-1}$ Mpc$^{-1}$)')
plt.title('Hubble values for small and its corresponding dipole')
plt.axis([min(np.append(largered, smallred))*0.9, max(np.append(largered, smallred))*1.1, min(np.append(large_H, small_H))*0.9, max(np.append(large_H, small_H))*1.1])
plt.legend()
plt.show()




##Histogram of dH
large_dH = []
small_dH = []

i = 0
while i < len(largered):
        large_dH.append(( fun.D(fun.HubbleIntegrate(largered[i])) - largedist[i])/fun.D(fun.HubbleIntegrate(largered[i])))
        i+=1
i = 0
while i < len(smallred):
        small_dH.append((fun.D(fun.HubbleIntegrate(smallred[i])) - smalldist[i])/fun.D(fun.HubbleIntegrate(smallred[i])))
        i+=1

print(ks_2samp(large_dH, small_dH))

largehist = plt.hist(large_dH, bins=binno, alpha=0.5, label='large', color='blue')
##Best fit for data
(largemu, largesigma) = norm.fit(large_dH)
emptylarge = plt.hist([], range=[-1,1], alpha=0.5, label='$N=%i, \mu=%.2f$, $\sigma=%.2f$'%(len(large_dH), largemu, largesigma), color='white')

##Plot of the histogram of the NLSS

smallhist = plt.hist(small_dH, bins=binno, alpha=0.5, label='small', color='red')
(smallmu, smallsigma) = norm.fit(small_dH)
emptysmall = plt.hist([], range=[-1,1], alpha=0.5, label='$N=%i, \mu=%.2f$, $\sigma=%.2f$'%(len(small_dH), smallmu, smallsigma), color='white')
plt.xlabel('% Deviation from expected distance')
plt.ylabel('Frequency')
plt.title('Histogram of deviation from expected distance for small and its corresponding dipole supernovae')
plt.legend()
plt.show()


large_dH.sort()
small_dH.sort()

plt.plot(large_dH, fun.cumfreq(large_dH), label="Maximal H0", color='blue')
plt.plot(small_dH, fun.cumfreq(small_dH), label="Minimal H0", color='red')
emptyNLSS = plt.hist([], range=[-0.3,0.3], alpha=0.5, label='Stat=%.3f \np-val=%.4f'%(ks_2samp(large_dH, small_dH)[0], ks_2samp(large_dH, small_dH)[1]), color='white')
plt.xlabel('$\delta d_L$')
plt.ylabel('Frequency')
plt.title('Cumulative frequency for the deviations \nof SNe found in the area of maximal and minimal expansion.')
plt.legend()
plt.show()
        





#########EQUITORIAL COORDINATE

fig = plt.figure()#figsize=(10, 5))
ax = fig.add_subplot(111, projection="mollweide")
ax.grid(True)

smallboxx = np.array([small[0]-error, small[0]-error, small[0]+error, small[0]+error, small[0]-error])*degtorad
smallboxy = np.array([small[1]-error, small[1]+error, small[1]+error, small[1]-error, small[1]-error])*degtorad

largeboxx = np.array([large[0]-error, large[0]-error, large[0]+error, large[0]+error, large[0]-error])*degtorad
largeboxy = np.array([large[1]-error, large[1]+error, large[1]+error, large[1]-error, large[1]-error])*degtorad

ax.plot(AllRA, Alldec, '.', color='black', label='All SNe for z<0.1')
ax.plot(largeRA, largedec, '.', color='blue', label='Maximal H0')
ax.plot(largeboxx, largeboxy, '-', color='blue')
ax.plot(smallRA, smalldec, '.', color='red', label='Minimal H0')
ax.plot(smallboxx, smallboxy, '-', color='red')
plt.xlabel('RA')
plt.ylabel('DEC')
plt.title('Positions for z<0.1 supernovae in equatorial coordinates')
plt.legend()

if ShowPosGraph:
        plt.show()



'''
#####GALACTIC COORDINATE TRANSFORMATION#######

NgalRA = []
Ngaldec = []
i = 0
while i < len(largeRA):
        x = SkyCoord(ra=largeRA[i]*u.radian, dec=largedec[i]*u.radian, frame='icrs')
        l = x.galactic.l.radian
        if l > np.pi:
                l = l-2*np.pi
        NgalRA.append(l)
        Ngaldec.append(x.galactic.b.radian)
        i+=1

SgalRA = []
Sgaldec = []
i = 0
while i < len(smallRA):
        x = SkyCoord(ra=smallRA[i]*u.radian, dec=smalldec[i]*u.radian, frame='icrs')
        l = x.galactic.l.radian
        SgalRA.append(l)
        Sgaldec.append(x.galactic.b.radian)
        i+=1



fig = plt.figure()#figsize=(10, 5))
ax = fig.add_subplot(111, projection="mollweide")
ax.grid(True)

ax.plot(NgalRA, Ngaldec, '.', color='blue', label='large')
ax.plot(SgalRA, Sgaldec, '.', color='red', label='small')
plt.xlabel('RA')
plt.ylabel('DEC')
plt.title('Galactic Coordinates')
plt.legend()

if ShowPosGraph:
        plt.show()

'''
