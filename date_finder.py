##import matplotlib.pyplot as plt ##import graphing
import re ##regex for separating the data

data = open("ListSne.txt", "r").read().split('\n') ##lines separated by \n
#data = data[5:] ##first five lines are not data, just comments

d = [['SN', 'Host Galaxy', 'Date', 'R/A Decl', 'Offset', 'Mag', 'Disc Ref', 'SN Pos', 'Posn Ref', 'Type', 'SN', 'Discoverers']]
i = 0

while i < len(data):
    d.append(re.split(r'\s{2,}', data[i])) ##regex splits every time if more than 2 spaces
    i+=1
