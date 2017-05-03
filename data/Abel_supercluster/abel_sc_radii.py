##import matplotlib.pyplot as plt ##import graphing
import re ##Regex
import pickle ##to write to file
import sys
sys.path.insert(0, 'C:\Python34\masters\graphs')
import functions as fun

sc = open("abel_supercluster_radii.txt", "r").read().split('\n')
scdata = []

abelname = fun.abelname
abelz = fun.abelz

for x in sc:
    scdata.append(x.split('\t'))

k=0
while k < len(scdata):
    i=0
    while i < len(scdata[k]):
        if scdata[k][i] == '':
            del scdata[k][i]
        else:
            i+=1
    if len(scdata[k]) < 6:
        del scdata[k]
    else:
        k+=1

k=0
while k < len(scdata):
    z=None
    for y in scdata[k]:
        if str(y)[-1] == 'z':
            name = y[:-1]
            if name in abelname:
                index = abelname.index(name)
                z = abelz[index]
    scdata[k].append(z)
    if z == None:
        del scdata[k]
    else:  
        k+=1

i=0
while i < len(scdata):
    scdata[i] = scdata[i][-7:]
    i+=1

i=0
for x in scdata:
    RA = x[:2]
    RA[0] = float(RA[0]) + float(RA[1])/60

    dec = x[2:4]
    dec[0] = float(dec[0]) + float(dec[1])/60
    del x[:-3]
    x.append(RA[0]*15)
    x.append(dec[0])

    x[0] = float(x[0])
    x[2] = float(x[2])
    i=i+1

i=0
with open("abel_sc_reduced.txt", "w") as text_file:
    while i < len(scdata)-1:
        text_file.write("{}\n".format(scdata[i]))
        i+=1
    if i == len(scdata)-1:
        text_file.write("{}".format(scdata[i]))
