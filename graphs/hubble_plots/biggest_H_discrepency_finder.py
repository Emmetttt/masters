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
loops = 10000
numberdetection = 7
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

##RA between -180 and 180
##Dec between -90 and 90
k = 0
maxH = 70
bigHRA = 0
bigHdec = 0
bignumber = 0
bigz = []
bigdist = []
bigdisterror = []

minH = 70
smallHRA = 0
smallHdec = 0
smallnumber = 0
smallz = []
smalldist = []
smalldisterror = []

while k < loops:
    pos = [(np.random.rand()*300)-150, (np.random.rand()*120)-60]
    SN = []
    degtorad = 0.0174533
    i = 0

    while i < len(redshifts):
        if redshifts[i] < maxredshift and redshifts[i] > minredshift:
            if dec[i] > (pos[1] - error) and dec[i] < (pos[1] + error) and RA[i] > (pos[0] - error) and RA[i] < (pos[0] + error):
                SN.append([redshifts[i], dist[i], disterror[i], float(RA[i])*degtorad, float(dec[i])*degtorad, name[i]])
        i+=1


    ##Plots
    SNred = np.array([x[0] for x in SN])
    SNdisterror = np.array([x[2] for x in SN])
    SNdist = np.array([x[1] for x in SN])
    SNRA = [x[3] for x in SN]
    SNdec = [x[4] for x in SN]

    if len(SNred) > numberdetection:

        H = fun.WLSF(SNred, SNdist, SNdisterror)
        if H > maxH:
            maxH = H
            bigHRA = pos[0]
            bigHdec = pos[1]
            bignumber = len(SNred)
            bigz = SNred
            bigdist = SNdist
            bigdisterror = SNdisterror

        if H < minH:
            minH = H
            smallHRA = pos[0]
            smallHdec = pos[1]
            smallnumber = len(SNred)
            smallz = SNred
            smalldist = SNdist
            smalldisterror = SNdisterror

    k+=1

print('Largest H: ', maxH, '\nBest RA: ', bigHRA, '\nBest dec: ', bigHdec, '\nNumber of SN: ', bignumber)
print('\nSmallest H: ', minH, '\nBest RA: ', smallHRA, '\nBest dec: ', smallHdec, '\nNumber of SN: ', smallnumber)

print('bigz = ', [x for x in bigz], '\nbigdist = ', [x for x in bigdist], '\nbigdisterror = ', [x for x in bigdisterror])
print('smallz = ', [x for x in smallz], '\nsmalldist = ', [x for x in smalldist], '\nsmalldisterror = ', [x for x in smalldisterror])
