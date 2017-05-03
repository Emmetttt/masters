##http://classic.sdss.org/supernova/snlist_confirmed.html
Hubbledata = open("hubble_data.txt", "r").read().split('\n')
NewHubbledata = []

i = 0
while i < len(Hubbledata):
    x = Hubbledata[i]
    x = x.split(' ')
    NewHubbledata.append(x)
    i+=1

NewHubbledata = NewHubbledata[:-1]

SDSSName = [x[0][2:] for x in NewHubbledata]
#RA
HubbleRA = [x[1] for x in NewHubbledata]
HubbleRA = [x.split(':') for x in HubbleRA]
HubbleRA = [float(x[0]) + (float(x[1])/60) + (float(x[2])/3600) for x in HubbleRA]
##dec
Hubbledec = [x[2] for x in NewHubbledata]
Hubbledec = [x.split(':') for x in Hubbledec]
Hubbledec = [float(x[0][0] + str(float(x[0][1:]) + (float(x[1])/60) + (float(x[2])/3600))) for x in Hubbledec]

HubblePos = [[x, y] for x, y in zip(HubbleRA, Hubbledec)]

