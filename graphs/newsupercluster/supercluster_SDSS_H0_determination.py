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

def deltaH(SNered, SNedist):
    return (fun.HubbleIntegrate(float(SNered)) - (float(SNedist)))  / fun.HubbleIntegrate(float(SNered))*100

def residual(data, dataerr, redshift, h, m, de):
    def integrand(z, m, de): ##Integral in luminosity distance equation
        return 1/np.sqrt((m*(1+z)**3)+de)

    resid = 0
    for i in range(len(data)):    
        integral = quad(integrand, 0, redshift[i], args=(m,de))
        model = (fun.c*(1+redshift[i]))/h * integral[0]
        resid = resid + (data[i]-model)**2/dataerr[i]**2 ##Sum of each residual

    return resid

##SUPERNOVA##
SNeRA = fun.RA
SNedec = fun.dec
SNedist = fun.dist
SNedisterror = fun.disterror
SNedistMpc = fun.distMpc
SNered = fun.redshifts
degtorad = 0.0174533

scindist =[231.20647901755987, 313.3285724315589, 64.26877173170195, 66.06934480075964, 155.59656316050751, 93.32543007969925, 586.1381645140297, 625.1726927756849, 209.89398836235287, 151.3561248436207, 190.54607179632444, 289.7343587701333, 145.2111617587745, 94.62371613657949, 91.20108393559097, 80.53784411990677, 102.32929922807536, 61.65950018614835, 2147.8304741305287, 972.7472237769641, 2322.736796357105, 3221.068791283441, 2322.736796357105, 3664.3757464783334, 3944.573020752785, 2766.9416454115126, 3019.951720402019, 2376.8402866248834]
scindisterror = [10.046995285934162, 13.37581952749451, 2.643251481327342, 4.262538374242557, 7.355788581002857, 5.088043533756917, 24.14575343003212, 20.849781955063893, 10.319835538165398, 7.167281677831622, 8.899129726751417, 12.424952399684088, 6.893576514658381, 4.883104617139995, 5.241441605493734, 4.898044386093373, 3.7953805705134944, 3.2691340304879826, 113.42359680958145, 29.226256097455106, 249.87605985194773, 424.0241004040261, 333.168079802597, 444.99658761530435, 110.13233189630856, 98.32770594923642, 206.55330163127016, 130.53325356344868]
scinz = [0.0529, 0.0627, 0.0163, 0.0152, 0.0357, 0.022, 0.1441, 0.1299, 0.0546, 0.0321, 0.0421, 0.0651, 0.0325, 0.0233, 0.0208, 0.0151, 0.0242, 0.015, 0.3804, 0.1897, 0.45, 0.615, 0.465, 0.453, 0.671, 0.496, 0.583, 0.51]

minresidual = 10000000
bestH = 0

i=0
while i < 200:
    h = 60 + i/10 ##H_0 known to be between 68 and 72 from previous tests
    k = 0
    res = residual(scindist, scindisterror, scinz, h, 0.2795, 0.7205)
    if res  < minresidual:
        minresidual = res
        bestH = h
        i+=1
    i+=1

###((3*10^5)/a)*(x + (0.5*(1 + 0.7205 - (0.2795/2)))*x^2 + (1/6)*(1+0.7205 - (0.2795/2) - (3*((0.2795/2)-0.7205)^2+1))*x^3)
###gives H0 = 72.76


redshiftrange = np.linspace(min(scinz), max(scinz), 100)
H = [fun.D(fun.HubbleIntegrateVarH(x, bestH)) for x in redshiftrange]
Hglobal = [fun.D(fun.HubbleIntegrateVarH(x, 69.79)) for x in redshiftrange]

##Graph for redshift vs distance
plt.figure()
plt.errorbar(scinz, scindist, yerr=scindisterror, fmt='.', markersize=5, color="blue", label='SNe in superclusters')
plt.plot(redshiftrange, H, '-', color="blue", label='H0 = 71.6')
plt.plot(redshiftrange, Hglobal, '-', color="grey", label='H0 = 69.79')

plt.legend(loc="upper left")
plt.xlabel('Redshift')
plt.ylabel('Distance (Mpc)')
plt.title('Hubble Diagrams for SGC and its corresponding dipole')
plt.axis([0, max(scinz)*1.1, 0, max(scindist)*1.1])
plt.legend()
plt.show()
