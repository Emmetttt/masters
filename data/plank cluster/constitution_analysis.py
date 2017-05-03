from astropy.io import fits
import pickle ##to write to file
import sys
sys.path.insert(0, 'C:\Python34\masters\graphs')
import functions as fun

hdulist = fits.open('HFI_PCCS_SZ-union_R2.08.fits')
data = hdulist[1].data

##taken from hdulist[1].header
##[0] = Index                            [5] = Dec
##[1] = NameError                        [6] = SNR
##[2] = Galactic Longitude               [7] = ID
##[3] = Galactic Latitutde               [8] = Redshift
##[4] = RA                         

cluster = []

i = 0
while i < len(data):
    if str(data[i][9]) == 'nan':
        i+=1
    elif data[i][18] < 0:
        i+=1
        print("skip")
    else:
        cluster.append([data[i][1], data[i][18], data[i][4], data[i][5]])
        i+=1

k = 0
with open("cluster_data2.txt", "w") as text_file:
    while k < len(cluster):
        text_file.write("{}\n".format(cluster[k]))
        k+=1
