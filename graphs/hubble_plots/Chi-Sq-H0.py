##import graphing
import numpy as np
from scipy.integrate import quad
##import functions
import sys
sys.path.insert(0, 'C:\Python34\masters\graphs')
import functions as fun

redshifts = fun.redshifts
distmod = fun.distmod
dist = fun.dist
disterror = fun.disterror

def residual(data, dataerr, redshift, h, m, de):
    def integrand(z, m, de): ##Integral in luminosity distance equation
        return 1/np.sqrt((m*(1+z)**3)+de)

    resid = 0
    for i in range(len(data)):    
        integral = quad(integrand, 0, redshift[i], args=(m,de))
        d = (fun.c*(1+redshift[i]))/h * integral[0]
        model = fun.dmod(d)+20
        resid = resid + (data[i]-model)**2/dataerr[i]**2 ##Sum of each residual

    return resid

minresidual = 10000000
bestH = 0
bestM = 0
bestDE = 0

i=0
while i < 20:
    h = 69.7 + i/100 ##H_0 known to be between 68 and 72 from previous tests
    k = 0
    while k < 50:
        m = (27.8 + k/100)/100 ##m known to be between 27 and 31 from previous tests
        de = 1 - m ## m + de = 1
        res = residual(dist, disterror, redshifts, h, m, de)
        if res  < minresidual:
            minresidual = res
            bestH = h
            bestM = m
            bestDE = de
        k+=1
    print(i)
    i+=1


##Error Calculation##
##Keeping either H0 or m constant, look for when chisq = 563.62570859524328
print('go')
res=0
i=0
while i < 10000:
    h = bestH + i/1000
    res = residual(dist, disterror, redshifts, h, bestM, bestDE)
    if res > minresidual + 6.17:
        break
    i+=1
print("H0 = ", bestH, ' ± ', i/1000)

res=0
i=0
while i < 4000:
    m = bestM + i/1000
    res = residual(dist, disterror, redshifts, bestH, m, bestDE)
    if res > minresidual + 6.17:
        break
    i+=1
print("Matter = ", bestM, ' ± ', i/1000)
print("Dark Energy = ", bestDE, ' ± ', i/1000)


'''
    
minerr = 10
dde = 0
i=0
while i < 1600:
    h = 69.58
    m = (25 + i/100)/100 
    de = 1 - m
    res = residual(dist, disterror, redshifts, h, m, de)/(len(dist)-1)
    print(np.abs(res - chisq))
    if np.abs(res - chisq) - 1 < minerr and np.abs(res - chisq) - 1 > 0:
        print(np.abs(res - chisq))
        dde = np.abs(m - bestM)
    i+=1
print('Best H0: ', bestH, ' ± ', dH, '\nBest Matter Density: ', bestM, ' ± ', dde, '\nBest Dark Energy Density: ', bestDE, ' ± ', dde, '\nChiSquared Statistic: ', minresidual/(len(dist)-1))
'''
