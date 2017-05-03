##import matplotlib.pyplot as plt ##import graphing
import re ##Regex
import pickle ##to write to file
import sys
sys.path.insert(0, 'C:\Python34\masters\graphs')
import functions as fun

sc = open("abel_supercluster.txt", "r").read().split('\n')
scdata = []

abelname = fun.abelname
abelz = fun.abelz

for x in sc:
    scdata.append(x.split('\t'))
for x in scdata:
    i=0
    while i < len(x):
        if x[i] == '':
            del x[i]
        else:
            i+=1

i=0

for x in scdata:
    for y in x:
        if y[-1] == 'z':
            name = y[:-1]
            if name in abelname:
                index = abelname.index(name)
                z = abelz[index]
    x.append(z)

while i < len(scdata):
    if len(scdata[i]) < 7:
        del scdata[i]
    else:
        scdata[i] = scdata[i][-5:]
        i+=1

i=0
for x in scdata:
    RA = x[:2]
    RA[0] = float(RA[0]) + float(RA[1])/60

    dec = x[2:4]
    dec[0] = float(dec[0]) + float(dec[1])/60
    del x[:-1]
    x.append(RA[0]*15)
    x.append(dec[0])

    x[0] = float(x[0])
    i=i+1

i=0
with open("abel_sc_reduced.txt", "w") as text_file:
    while i < len(scdata)-1:
        text_file.write("{}\n".format(scdata[i]))
        i+=1
    if i == len(scdata)-1:
        text_file.write("{}".format(scdata[i]))

