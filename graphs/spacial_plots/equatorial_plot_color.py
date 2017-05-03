import numpy as np
import matplotlib.pyplot as plt
##import functions
import sys
sys.path.insert(0, 'C:\Python34\masters\graphs')
import functions as fun

RA = fun.RA
dec = fun.dec
dist = fun.distMpc
disterror = fun.disterrMpc
redshifts = fun.redshifts

##Loop to convert the RA into -180<RA<180, so they're smoother
j=0
'''
while j < len(RA):
        if RA[j] > 180:
                RA[j] = RA[j] - 360
                #RA[j] = -RA[j]
        j+=1
'''
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

k = 0

while k < len(dec):
        if redshifts[k] < 0.1:
                lowZ.append([RA[k], dec[k], dist[k], disterror[k], redshifts[k]])
        elif redshifts[k] < 0.5:
                midZ.append([RA[k], dec[k], dist[k], disterror[k], redshifts[k]])
        elif redshifts[k] > 0.5:
                highZ.append([RA[k], dec[k], dist[k], disterror[k], redshifts[k]])
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


##Superclusters
sc = fun.sc
lowZsc = []
midZsc = []
highZsc = []

k = 0
while k < len(sc):
        if float(sc[k][3]) < 0.1:
                lowZsc.append([float(sc[k][1]), float(sc[k][2]), float(sc[k][3])])
        elif float(sc[k][3]) < 0.5:
                midZsc.append([float(sc[k][1]), float(sc[k][2]), float(sc[k][3])])
        else:
                highZsc.append([float(sc[k][1]), float(sc[k][2]), float(sc[k][3])])
        k+=1




plt.scatter(extract(lowZ, 0), extract(lowZ, 1), s=5, c=extract(Array, 0), vmin=-30., vmax=30., cmap=plt.cm.seismic)
#plt.scatter(extract(lowZsc, 0), extract(lowZsc, 1), color="yellow")
plt.grid(True)

#for i, txt in enumerate(extract(Array, 0)):
#        plt.annotate(txt, (extract(lowZ, 0)[i], extract(lowZ, 1)[i]))

plt.xlabel('RA')
plt.ylabel('DEC')
plt.title('Deviation from H0=68, 0 < z < 0.1')
plt.show()
