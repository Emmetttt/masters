EssenceData = open("essence_data.txt", "r").read().split('\n')
NewEdata = []

i = 0
while i < len(EssenceData):
    x = EssenceData[i]
    x = x.split(' ')
    NewEdata.append(x)
    i+=1

NewEdata = NewEdata[1:]

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
