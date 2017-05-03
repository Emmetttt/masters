import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
##import functions
import sys
sys.path.insert(0, 'C:\Python34\masters\graphs')
import functions as fun


RA = fun.RA
dec = fun.dec
dist = fun.dist
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

##Loop to sort the data into deviation from mean
data = []
k = 0

while k < len(dec):
        data.append([RA[k], dec[k], dist[k], disterror[k], redshifts[k]])
        k+=1

#RA = [0], dec = [1], dist = [2], disterr = [3], redshifts = [4]
#lowZ, midZ, highZ
extract = lambda x,y: [b[y] for b in x] ## extract(highZ, 1) for highZ declinations


##Test against H0
H0 = 68
Array = [] ##H0 deviation, RA, dec

l = 0

while l < len(data):
        #H_calc = round((extract(data, 4)[l]/extract(data, 2)[l]) * fun.c) ##Finds H for given point
        Array.append([extract(data,2)[l] - fun.expectedD(extract(data, 4)[l]), extract(data, 0)[l], extract(data, 1)[l]]) ## H0 deviation, RA, dec
        l+=1 ##H_calculated - expected distance modulus at z





##x, y, size of marker, color gradients from deviation of H0, min gradient, max gradient, color scheme.
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.scatter(extract(data, 0), extract(data, 1), extract(data, 4), s=50, c=extract(Array, 0), vmin=-3., vmax=3., cmap=plt.cm.seismic)
ax.grid(True)

#for i, txt in enumerate(extract(Array, 0)):
#        plt.annotate(txt, (extract(lowZ, 0)[i], extract(lowZ, 1)[i]))

plt.xlabel('RA')
plt.ylabel('DEC')
plt.title('Deviation from H0=68')
plt.show()


