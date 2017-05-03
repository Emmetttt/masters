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

def HubbleCalculate(z, dist):
    ##Calculates the distance modulus to an object at z redshift. Assesses the integral
    ##between 0 and z
    ##To find the value in MPc, use fun.D(fun.HubbleIntegrate(0.2))
    def integrand(z, m, de):
        return 1/np.sqrt((m*(1+z)**3)+de)
    
    m=0.22
    de=0.7205
    H0 = 69.79
    a = 3*10**5/dist
    I = quad(integrand, 0, z, args=(m,de))
    H = I[0]*a*(1+z)
    return H



d = fun.D(fun.HubbleIntegrate(0.06))
#d = 260.90498780322844

x=HubbleCalculate(0.06, d)
print(x)
