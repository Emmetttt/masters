##import matplotlib.pyplot as plt ##import graphing
import re ##Regex
import pickle ##to write to file

sc = open("SDSS_Cluster.txt", "r").read().split('\n')
scdata = []


for x in sc:
    scdata.append(x.split(' '))

for x in scdata:
    i=0
    while i < len(x):
        if x[i] == '':
            del x[i]
        else:
            i+=1

for x in scdata:
    x[1] = float(x[1])
    x[2] = float(x[2])
    x[3] = float(x[3])
    x[4] = float(x[4])

i=0

with open("SDSS_cluster_reduced.txt", "w") as text_file:
    while i < len(scdata)-1:
        text_file.write("{}\n".format(scdata[i]))
        i+=1
    if i == len(scdata)-1:
        text_file.write("{}".format(scdata[i]))
