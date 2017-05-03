##import matplotlib.pyplot as plt ##import graphing
import re ##Regex
import pickle ##to write to file
import numpy as np

data = open("AllSNe.txt", "r").read().split('\\') ##lines separated by \n
#data = data[5:] ##first five lines are not data, just comments

MainDataSet = [['name', 'redshift', 'B-Band Mag', 'Strech', 'Color', 'dist mod', 'sample', 'Failed']]
i = 0
j = 0

while i < len(data): ##Tested, works well. Reduces lines from 2244 to 911, preserving all data
    if data[i] == '' or data[i] == 'nodata': ##Deletes missing cells
        del data[i]
    else:
        i+=1

while j < len(data): ##Tested, works well. Just appends all of data into MainDataSet
    MainDataSet.append(data[j].split(' & ')) ##data points are separated by \t
    if j > 0: #won't affect the 'name' header
        MainDataSet[j][0] = MainDataSet[j][0][2:] ##First datapoint comes as ' \n1992br', just trims the beginning
        j+=1
    else:
        j+=1

print("length: ", len(MainDataSet)) ##Length should be 911 + 1, true.

k = 0
delcount = 0
while k < len(MainDataSet): ##Tested, works well. Still has entries such as ['data ', '0.128', '6', 'f']
    if MainDataSet[k] == ['data ', ' ']:
        del MainDataSet[k]
        delcount+=1
    else:
        k+=1

print("length: ", len(MainDataSet), ", ", delcount, "empty cells deleted") ##Length should be 911 + 1, true.

k = 0
newdelcount = 0
while k < len(MainDataSet):
    if MainDataSet[k][0] == 'data ':
        del MainDataSet[k]
    else:
        k+=1

print("length: ", len(MainDataSet), ", ", newdelcount, "cells with 'datacount' and no data as name") ##Length should be 911 + 1, true.


#######Delete all that are not in the union2.1 Data set
Union2names = open("C:\\Python34\\masters\\data\\union2\\Union2_names.txt", "r").read().split("\n")

i = 0
while i < len(MainDataSet):
    if MainDataSet[i][0] in Union2names:
        i+=1
    else:
        del MainDataSet[i]


###Harvard DB positions

datesList = open("ListSne.txt", "r").read().split('\n') ##lines separated by \n
harvardDB = [['SN', 'Host Galaxy', 'Date', 'R/A Decl', 'Offset', 'Mag', 'Disc Ref', 'SN Pos', 'Posn Ref', 'Type', 'SN', 'Discoverers']]
i = 0

while i < len(datesList):
    harvardDB.append(re.split(r'\s{2,}', datesList[i])) ##regex splits every time if more than 2 spaces
    i+=1
    
####### Append main data set with name of galaxy, coords

k = 0 ##loops over all MainDataSet cells
l = 0 ##loops over harvards data
fail = [] ##counts failures


while k < len(MainDataSet):
    while l < len(harvardDB):
        if MainDataSet[k][0] == harvardDB[l][0] or MainDataSet[k][0].upper() == harvardDB[l][0] or MainDataSet[k][0].lower() == harvardDB[l][0]:
            MainDataSet[k].append(harvardDB[l][3]) ##appends coords
            MainDataSet[k].append(harvardDB[l][1]) ##appends galaxy name

            ##Loop to get SN type
            if len(harvardDB[l]) == 13:
                MainDataSet[k].append(harvardDB[l][10]) ##SN type
                #print(harvardDB[l][10], l)
            elif len(harvardDB[l]) == 12: ##Sometimes data is missing so length can be shorter
                if 'AJ' in harvardDB[l][9] or 'IAUC' in harvardDB[l][9]: ##Puts no data for SN with no type
                    MainDataSet[k].append('')
                else:
                    MainDataSet[k].append(harvardDB[l][9]) ##Appends the type of supernova to the main database
            elif len(harvardDB[l]) == 11:
                MainDataSet[k].append(harvardDB[l][8])
            elif len(harvardDB[l]) == 14:
                MainDataSet[k].append(harvardDB[l][11])
            l = len(harvardDB)

            
        elif l == len(harvardDB)-1:
            #print("No entry found for ", MainDataSet[k][0])
            fail.append(MainDataSet[k])
            l+=1
        else:
            l+=1
    k+=1
    l=0

#print("length: ", len(MainDataSet), ", ", fail, "cells with no position in the harvardDB") ##Length should be 911 + 1, true.

PosDataSet = [] ##All data with positions

m = 0
e = 0 ##How many passed
while m < len(MainDataSet): ##Appends new list with all complete data
    if len(MainDataSet[m]) > 10:            ##Only adds if number of data points if > 10. Might seem like it
        PosDataSet.append(MainDataSet[m])   ##deletes points w/ relevent data but it is okay
        e+=1
    #else:
        #print(MainDataSet[m])
    m+=1

m = 0
newcoords = []

def Errors(n): ##Makes the data python friendly
    n = n.replace("(", ",")
    n = n.replace(")", "")
    n = n.split(",")
    z = []
    for x in n:
        z.append(float(x))
    return z

while m < len(PosDataSet): ##Makes all the errors and positions readable
    ##ERRORS##
    PosDataSet[m][2] = Errors(PosDataSet[m][2])
    PosDataSet[m][3] = Errors(PosDataSet[m][3])
    PosDataSet[m][4] = Errors(PosDataSet[m][4])
    PosDataSet[m][5] = Errors(PosDataSet[m][5])
    
    ##POSITION##
    x = PosDataSet[m][9]
    x = [x[:7], x[8:]] ##separate RA and Dec
    ##
    RA = x[0] ##start manipulating RA
    y = RA[:2]
    z = RA[3:]
    RA = float(y) + (float(z)/60) ##decimalises the minute
    RA = round(RA,4) ##4 decimal points
    ##
    DEC = x[1]
    y = DEC[1:3]
    z = DEC[4:]
    D = float(y) + (float(z)/60)
    D = round(D,4)
    if DEC[0] == "-":
        D = 0 - D
    PosDataSet[m][9] = [RA, D]
    m+=1
    
i=0

#####Position data for Legacy Survey
LegacyData = open("legacy_names_pos.txt", "r").read().split('\n')
LegacyNewdata = []

for d in LegacyData:
    LegacyNewdata.append(d.split('\t'))

LegacyNewdata.append(['', 'SNLS 03D4ag', '', '22 14 45.808', '-17 44 22.99', ' ', ' ', ' ', ' ', '23.44', 'SNIa', '', '']) ##http://www.rochesterastronomy.org/sn2003/#2003gyhb
LegacyNewdata.append(['', 'SNLS 03D3ba', '', '14 16 33.448', '+52 20 32.19', ' ', ' ', ' ', ' ', '23.44', 'SNIa', '', ''])
LegacyNewdata.append(['', 'SNLS 03D3cd', '', '14 18 29.968', '+52 36 44.19', ' ', ' ', ' ', ' ', '23.44', 'SNIa', '', ''])
LegacyNewdata.append(['', 'SNLS 03D3cc', '', '14 19 25.198', '+52 32 25.79', ' ', ' ', ' ', ' ', '23.44', 'SNIa', '', ''])

LegName = [x[1] for x in LegacyNewdata]
LegName = [x[5:] for x in LegName]
LegType = [x[-3][2:] for x in LegacyNewdata]
#RA
LegRA = [x[3] for x in LegacyNewdata]
LegRA = [x.split(" ") for x in LegRA]
LegRA = [float(x[0]) + (float(x[1])/60) + (float(x[2])/3600) for x in LegRA]
#Dec
LegDec = [x[4] for x in LegacyNewdata]
LegDec = [x.split(" ") for x in LegDec]
LegDec = [float(x[0][0] + str(float(x[0][1:]) + (float(x[1])/60) + (float(x[2])/3600))) for x in LegDec]
#Pos
LegPos = [[x, y] for x, y in zip(LegRA, LegDec)]

u = 0
h = 0
success = 0
fail2 = []

while u < len(LegName):
    while h < len(fail):
        if LegName[u] == fail[h][0]:
            fail[h].append(LegPos[u])
            fail[h].append('')
            fail[h].append(LegType[u])
            #Fix errors#
            fail[h][2] = Errors(fail[h][2])
            fail[h][3] = Errors(fail[h][3])
            fail[h][4] = Errors(fail[h][4])
            fail[h][5] = Errors(fail[h][5])

            PosDataSet.append(fail[h])
            success += 1
            del fail[h]
            h = len(fail)
        elif h == len(fail)-1:
            fail2.append(fail[h])
        h+=1
    h=0
    u+=1
    

EssenceData = open("essence_data.txt", "r").read().split('\n')
NewEdata = []

i = 0
while i < len(EssenceData):
    x = EssenceData[i]
    x = x.split(' ')
    NewEdata.append(x)
    i+=1

NewEdata = NewEdata[1:]
NewEdata.append(['-', '02:04:56.09', '-03:49:03.67', 'p534', 'Ia', 'Ia-norm', '73.6', '82.4', '17.6', '0.0', '0.587', '0.590', '0.006'])
NewEdata.append(['-', '02:07:04.66', '-03:28:04.37', 'p528', 'Ia', 'Ia-norm', '73.6', '82.4', '17.6', '0.0', '0.587', '0.590', '0.006'])


EssName = [x[3] for x in NewEdata]
#RA
EssRA = [x[1] for x in NewEdata]
EssRA = [x.split(':') for x in EssRA]
EssRA = [float(x[0]) + (float(x[1])/60) + (float(x[2])/3600) for x in EssRA]
##dec
Essdec = [x[2] for x in NewEdata]
Essdec = [x.split(':') for x in Essdec]
Essdec = [float(x[0][0] + str(float(x[0][1:]) + (float(x[1])/60) + (float(x[2])/3600))) for x in Essdec]

EssPos = [[x, y] for x, y in zip(EssRA, Essdec)]

u = 0
h = 0
success = 0
fail2 = []

while u < len(EssName):
    while h < len(fail):
        if EssName[u] == fail[h][0]:
            fail[h].append(EssPos[u])
            fail[h].append('')
            fail[h].append('.')
            #Fix errors#
            fail[h][2] = Errors(fail[h][2])
            fail[h][3] = Errors(fail[h][3])
            fail[h][4] = Errors(fail[h][4])
            fail[h][5] = Errors(fail[h][5])

            PosDataSet.append(fail[h])
            success += 1
            del fail[h]
            h = len(fail)
        elif h == len(fail)-1:
            fail2.append(fail[h])
        h+=1
    h=0
    u+=1


##http://classic.sdss.org/supernova/snlist_confirmed.html
SDSSdata = open("SDSS_data.txt", "r").read().split('\n')
NewSDSSdata = []

i = 0
while i < len(SDSSdata):
    x = SDSSdata[i]
    x = x.split(' ')
    NewSDSSdata.append(x)
    i+=1

NewSDSSdata = NewSDSSdata[:-1]

SDSSName = [x[0][2:] for x in NewSDSSdata]
#RA
SDSSRA = [x[1] for x in NewSDSSdata]
SDSSRA = [x.split(':') for x in SDSSRA]
SDSSRA = [float(x[0]) + (float(x[1])/60) + (float(x[2])/3600) for x in SDSSRA]
##dec
SDSSdec = [x[2] for x in NewSDSSdata]
SDSSdec = [x.split(':') for x in SDSSdec]
SDSSdec = [float(x[0][0] + str(float(x[0][1:]) + (float(x[1])/60) + (float(x[2])/3600))) for x in SDSSdec]

SDSSPos = [[x, y] for x, y in zip(SDSSRA, SDSSdec)]

u = 0
h = 0
success = 0
fail2 = []

while u < len(SDSSName):
    while h < len(fail):
        if SDSSName[u] == fail[h][0]:
            fail[h].append(SDSSPos[u])
            fail[h].append('')
            fail[h].append('.')
            #Fix errors#
            fail[h][2] = Errors(fail[h][2])
            fail[h][3] = Errors(fail[h][3])
            fail[h][4] = Errors(fail[h][4])
            fail[h][5] = Errors(fail[h][5])

            PosDataSet.append(fail[h])
            success += 1
            del fail[h]
            h = len(fail)
        elif h == len(fail)-1:
            fail2.append(fail[h])
        h+=1
    h=0
    u+=1




ReissData = open("reiss_data.txt", "r").read().split('\n')
NewReissdata = []

i = 0
while i < len(ReissData):
    x = ReissData[i]
    x = x.split(' ')
    NewReissdata.append(x)
    i+=1


NewReissdata = NewReissdata[:-1]

ReissName = [x[0][3:] for x in NewReissdata]
#RA
ReissRA = [x[3] for x in NewReissdata]
ReissRA = [x.split(':') for x in ReissRA]
ReissRA = [float(x[0]) + (float(x[1])/60) + (float(x[2])/3600) for x in ReissRA]
##dec
Reissdec = [x[4] for x in NewReissdata]
Reissdec = [x.split(':') for x in Reissdec]
Reissdec = [float(x[0][0] + str(float(x[0][1:]) + (float(x[1])/60) + (float(x[2])/3600))) for x in Reissdec]

ReissPos = [[x, y] for x, y in zip(ReissRA, Reissdec)]

u = 0
h = 0
success = 0
fail2 = []

while u < len(ReissName):
    while h < len(fail):
        if ReissName[u] == fail[h][0]:
            fail[h].append(ReissPos[u])
            fail[h].append('')
            fail[h].append('.')
            #Fix errors#
            fail[h][2] = Errors(fail[h][2])
            fail[h][3] = Errors(fail[h][3])
            fail[h][4] = Errors(fail[h][4])
            fail[h][5] = Errors(fail[h][5])

            PosDataSet.append(fail[h])
            success += 1
            del fail[h]
            h = len(fail)
        elif h == len(fail)-1:
            fail2.append(fail[h])
        h+=1
    h=0
    u+=1

##https://arxiv.org/pdf/1010.5786.pdf
Hubbledata = open("hubble_data.txt", "r").read().split('\n')
NewHubbledata = []

i = 0
while i < len(Hubbledata):
    x = Hubbledata[i]
    x = x.split(' ')
    NewHubbledata.append(x)
    i+=1

NewHubbledata.append(['2003XX', '12:37:29.00', '+62:11:27.8'])

HubbleName = [x[0] for x in NewHubbledata]
#RA
HubbleRA = [x[1] for x in NewHubbledata]
HubbleRA = [x.split(':') for x in HubbleRA]
HubbleRA = [float(x[0]) + (float(x[1])/60) + (float(x[2])/3600) for x in HubbleRA]
##dec
Hubbledec = [x[2] for x in NewHubbledata]
Hubbledec = [x.split(':') for x in Hubbledec]
Hubbledec = [float(x[0][0] + str(float(x[0][1:]) + (float(x[1])/60) + (float(x[2])/3600))) for x in Hubbledec]

HubblePos = [[x, y] for x, y in zip(HubbleRA, Hubbledec)]

u = 0
h = 0
success = 0
fail2 = []

while u < len(HubbleName):
    while h < len(fail):
        if HubbleName[u] == fail[h][0]:
            fail[h].append(HubblePos[u])
            fail[h].append('')
            fail[h].append('.')
            #Fix errors#
            fail[h][2] = Errors(fail[h][2])
            fail[h][3] = Errors(fail[h][3])
            fail[h][4] = Errors(fail[h][4])
            fail[h][5] = Errors(fail[h][5])

            PosDataSet.append(fail[h])
            success += 1
            del fail[h]
            h = len(fail)
        elif h == len(fail)-1:
            fail2.append(fail[h])
        h+=1
    h=0
    u+=1




[print(x[0]) for x in fail]

PosDataSet = [['name', 'redshift', 'B-Band Mag', 'Strech', 'Color', 'dist mod', 'sample', 'Failed', '?', 'Position', 'Sn Type']] + PosDataSet

i=0
with open("Reduced_Data_update.txt", "w") as text_file:
    while i < len(PosDataSet)-1:
        text_file.write("{}\n".format(PosDataSet[i]))
        i+=1
    if i == len(PosDataSet)-1:
        text_file.write("{}".format(PosDataSet[i]))
