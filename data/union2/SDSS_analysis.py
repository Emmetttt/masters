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

