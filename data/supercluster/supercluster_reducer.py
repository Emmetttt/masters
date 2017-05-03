##import matplotlib.pyplot as plt ##import graphing
import re ##Regex
import pickle ##to write to file

sc = open("superclusters.txt", "r").read().split('\n')
scdata = []

for x in sc:
    scdata.append(x.split(' '))

a = ['Supercluster', 'RA', 'Dec', 'z', 'Rmax (Mpc)', 'L_X', 'Multiplicity', 'f=50', 'ID']
scdata = [a] + scdata

i=0

with open("sc_reduced.txt", "w") as text_file:
    while i < len(scdata)-1:
        text_file.write("{}\n".format(scdata[i]))
        i+=1
    if i == len(scdata)-1:
        text_file.write("{}".format(scdata[i]))

