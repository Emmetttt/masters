##import matplotlib.pyplot as plt ##import graphing
import re ##Regex
import pickle ##to write to file

sc = open("cluster_south.txt", "r").read().split('\n')
scdata = []


for x in sc:
    scdata.append(x.split('\t'))

for x in scdata:
    i=0
    while i < len(x):
        if x[i] == '':
            del x[i]
        else:
            i+=1

scdata = scdata[:-1]

for x in scdata:
    RA = x[1]
    RA = RA.split(' ')
    RA[1] = float(RA[1]) + float(RA[2])/60
    RA[0] = float(RA[0]) + RA[1]/60
    x[1] = RA[0]*15
    

    dec = x[2]
    dec = dec.split(' ')
    dec[1] = float(dec[1]) + float(dec[2])/60
    dec[0] = float(dec[0]) + dec[1]/60
    x[2] = dec[0]

    x[3] = float(x[3])

i=0

with open("cluster_south_reduced.txt", "w") as text_file:
    while i < len(scdata)-1:
        text_file.write("{}\n".format(scdata[i]))
        i+=1
    if i == len(scdata)-1:
        text_file.write("{}".format(scdata[i]))

