import re ##Regex
import pickle ##to write to file
import numpy as np

data = [eval(x) for x in open("C:\\Python34\\masters\\data\\union2\\Reduced_Data.txt", "r").read().split("\n")]
Union2names = open("C:\\Python34\\masters\\data\\union2\\Union2_names.txt", "r").read().split("\n")

'''
EssenceData = open("essence_data.txt", "r").read().split('\n')[:-1]
EssenceNames = []

i = 0
while i < len(EssenceData):
    x = EssenceData[i].split('\t')
    EssenceNames.append(x[2])
    i+=1
'''

i = 0
gooddata = []
gooddatanames = []
baddata = []

while i < len(data):
    if data[i][0] in Union2names:
        gooddata.append(data[i])
        gooddatanames.append(data[i][0])
        Union2names.remove(data[i][0])
    else:
        baddata.append(data[i])
    i+=1

