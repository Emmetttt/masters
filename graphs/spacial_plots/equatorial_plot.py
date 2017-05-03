import numpy as np
import matplotlib.pyplot as plt
from astropy import units as u
from astropy.coordinates import SkyCoord
##import functions
import sys
import astropy.coordinates as coord
from astropy.io import ascii
sys.path.insert(0, 'C:\Python34\masters\graphs')
import functions as fun

RA = fun.RA
dec = fun.dec
dist = fun.distMpc
disterror = fun.disterrMpc
redshifts = fun.redshifts

##Loop to convert the RA into -180<RA<180, so they're smoother
j=0
while j < len(RA):
        if RA[j] > 180:
                RA[j] = RA[j] - 360
                #RA[j] = -RA[j]
        j+=1

##Loop to find bad coords which are just on dec~0
BadCoord = []
i=0

while i < len(dec): ##Attempt to find all  that lay on ~~0 range
	if dec[i] > -5 and dec[i] < 5:
		BadCoord.append([i, dec[i]])
	i+=1

##Loop to split into low z < 0.1, mid z < 0.5, large z > 0.5
lowZ = []
midZ = []
highZ = []
degtorad = 0.0174533

k = 0

while k < len(dec):
        if redshifts[k] < 0.2:
                lowZ.append([RA[k]*degtorad, dec[k]*degtorad, dist[k], disterror[k], redshifts[k]])
        elif redshifts[k] < 0.5:
                midZ.append([RA[k]*degtorad, dec[k]*degtorad, dist[k], disterror[k], redshifts[k]])
        elif redshifts[k] > 0.5:
                highZ.append([RA[k]*degtorad, dec[k]*degtorad, dist[k], disterror[k], redshifts[k]])
        k+=1

#RA = [0], dec = [1], dist = [2], disterr = [3], redshifts = [4]
#lowZ, midZ, highZ
extract = lambda x,y: [b[y] for b in x] ## extract(highZ, 1) for highZ declinations


##Test against H0
H0 = 68
Array = []

l = 0

while l < len(lowZ):
        H_calc = round((extract(lowZ, 4)[l]/extract(lowZ, 2)[l]) * fun.c) ##Finds H for given point
        Array.append([H_calc - H0, extract(lowZ, 0), extract(lowZ, 1)]) ## H0 deviation, RA, dec
        l+=1


boxx = np.array([335-360, 295-360, 335-360, 295-360, 335-360]) * 0.0174533
boxy = np.array([-30, 30, 30, -30, -30]) * 0.0174533

CelRA = []
Celdec = []
i=0
while i < len(boxx):
        x = SkyCoord(l=(boxx[i])*u.radian, b=boxy[i]*u.radian, frame='galactic')
        l = x.icrs.ra.radian
        if l > np.pi:
                l = l-2*np.pi
        CelRA.append(l)
        Celdec.append(x.icrs.dec.radian)
        i+=1

lx = np.array([129, 134, 165, 126, 127, 126, 129, 146]) * 0.0174533
bx = np.array([-18, -28, -19, -13, 14.3, 18.2, 8.6, -12.02, 18]) * 0.0174533
newRA =[]
newdec=[]
i=0
while i < len(lx):
        x = SkyCoord(l=lx[i]*u.radian, b=bx[i]*u.radian, frame='galactic')
        l = x.icrs.ra.radian
        if l > np.pi:
                l = l-2*np.pi
        newRA.append(l)
        newdec.append(x.icrs.dec.radian)
        i+=1
xbox = np.array([0, 0, 75, 75, 0])* 0.0174533
ybox = np.array([20, 90, 90, 20, 20])* 0.0174533

fig = plt.figure()#figsize=(10, 5))
ax = fig.add_subplot(111, projection="mollweide")
ax.grid(True)
ax.plot(extract(lowZ, 0), extract(lowZ, 1), '.')
#ax.plot(newRA, newdec, 'x', color='red')
ax.plot(xbox, ybox, '--', color='red')
for i, txt in enumerate(['A', 'B', 'C', 'D', 'E', 'F', 'H']):
        plt.annotate(txt, (newRA[i], newdec[i]))


#ax.scatter(extract(lowZ, 0), extract(lowZ, 1), s=10, c=extract(Array, 0), vmin=-30., vmax=30., cmap=plt.cm.seismic)
#ax.colorbar()
#for i, txt in enumerate(extract(Array, 0)):
#        plt.annotate(txt, (extract(lowZ, 0)[i], extract(lowZ, 1)[i]))

plt.xlabel('RA')
plt.ylabel('DEC')
plt.title('Spacial distribution of z<0.2 supernovae with annotated location of dipole')
plt.show()

RA = extract(lowZ, 0)
dec = extract(lowZ, 1)

'''
SgalRA = []
Sgaldec = []
i = 0
while i < len(RA):
        x = SkyCoord(ra=(RA[i])*u.radian, dec=dec[i]*u.radian, frame='icrs')
        l = x.galactic.l.radian
        if l > np.pi:
                l = l-2*np.pi
        SgalRA.append(l)
        Sgaldec.append(x.galactic.b.radian)
        i+=1


fig = plt.figure()#figsize=(10, 5))
ax = fig.add_subplot(111, projection="mollweide")
ax.grid(True)

ax.plot(SgalRA, Sgaldec, '.', color='black')
ax.plot(boxx, boxy, '-', color='red')
plt.xlabel('RA')
plt.ylabel('DEC')
plt.title('Galactic Coordinates')
plt.legend()
plt.show()
'''
