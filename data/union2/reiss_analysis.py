ReissData = open("reiss_data.txt", "r").read().split('\n')
NewReissdata = []

i = 0
while i < len(ReissData):
    x = ReissData[i]
    x = x.split(' ')
    NewReissdata.append(x)
    i+=1


NewReissdata = NewReissdata[:-1]

ReissName = [x[1][2:] for x in NewReissdata]
#RA
ReissRA = [x[3] for x in NewReissdata]
ReissRA = [x.split(':') for x in ReissRA]
ReissRA = [float(x[0]) + (float(x[1])/60) + (float(x[2])/3600) for x in ReissRA]
##dec
Reissdec = [x[4] for x in NewReissdata]
Reissdec = [x.split(':') for x in Reissdec]
Reissdec = [float(x[0][0] + str(float(x[0][1:]) + (float(x[1])/60) + (float(x[2])/3600))) for x in Reissdec]

ReissPos = [[x, y] for x, y in zip(ReissRA, Reissdec)]

