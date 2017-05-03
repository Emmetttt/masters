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

##Loop to split into low z < 0.1, mid z < 0.5, large z > 0.5
#RA = [0], dec = [1], dist = [2], disterr = [3], redshifts = [4]
lowZ = []
midZ = []
highZ = []

k = 0

while k < len(dec):
        if redshifts[k] < 0.05:
                lowZ.append([RA[k], dec[k], dist[k], disterror[k], redshifts[k]])
        elif redshifts[k] < 0.5:
                midZ.append([RA[k], dec[k], dist[k], disterror[k], redshifts[k]])
        elif redshifts[k] > 0.5:
                highZ.append([RA[k], dec[k], dist[k], disterror[k], redshifts[k]])
        k+=1


extract = lambda x,y: [b[y] for b in x] ## extract(highZ, 1) for highZ declinations


##Test against H0
H0 = 68
H_dev = []

l = 0

test = lowZ

while l < len(lowZ):
        H_calc = round((extract(test, 4)[l]/extract(test, 2)[l]) * fun.c) ##Finds H for given point
        H_dev.append(H_calc - H0) ## H0 deviation, RA, dec
        l+=1



##x, y, size of marker, color gradients from deviation of H0, min gradient, max gradient, color scheme.
plt.scatter(extract(test, 0), extract(test, 1), s=50, c=H_dev, vmin=-30., vmax=30., cmap=plt.cm.seismic)
plt.colorbar()
plt.grid(True)

#for i, txt in enumerate(extract(Array, 0)):
#        plt.annotate(txt, (extract(lowZ, 0)[i], extract(lowZ, 1)[i]))

plt.xlabel('RA')
plt.ylabel('DEC')
plt.title('Deviation from H0=68')
plt.show()
