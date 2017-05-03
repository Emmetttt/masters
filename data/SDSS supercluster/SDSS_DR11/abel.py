##import matplotlib.pyplot as plt ##import graphing
import re ##Regex
import pickle ##to write to file

##Import CMASS_Norths
CMASS_North = open("C:\Python34\masters\data\SDSS supercluster\SDSS_DR11\CMASS_North\Clusters_info.txt", "r").read().split('\n')
CMASS_North = CMASS_North[2:-1]
CMASS_North_pos = open("C:\Python34\masters\data\SDSS supercluster\SDSS_DR11\CMASS_North\Clusters_skypos.txt", "r").read().split('\n')
CMASS_North_pos = CMASS_North_pos[1:-1]
cluster = [] ##redshift, RA, dec, name, radius

for x in CMASS_North_pos:
    x = x.split('\t')
    cluster.append([x[3], x[1], x[2], x[0]])

i = 0
while i < len(CMASS_North):
    z = CMASS_North[i].split('\t')
    cluster[i].append(z[6]) ##appends radius
    i+=1

##Import CMASS_South
CMASS_South = open("C:\Python34\masters\data\SDSS supercluster\SDSS_DR11\CMASS_South\Clusters_info.txt", "r").read().split('\n')
CMASS_South = CMASS_North[1:-1]
CMASS_South_pos = open("C:\Python34\masters\data\SDSS supercluster\SDSS_DR11\CMASS_South\Clusters_skypos.txt", "r").read().split('\n')
CMASS_South_pos = CMASS_North_pos[1:-1]

for x in CMASS_South_pos:
    x = x.split('\t')
    cluster.append([x[3], x[1], x[2], x[0]])

##KEEP i AS COUNTER FROM ABOVE
j = 0
while j < len(CMASS_South):
    z = CMASS_South[j].split('\t')
    cluster[i].append(z[6]) ##appends radius
    i+=1
    j+=1



##LOWZ North

##Import CMASS_South
LOWZ_North = open("C:\Python34\masters\data\SDSS supercluster\SDSS_DR11\LOWZ_North\Clusters_info.txt", "r").read().split('\n')
LOWZ_North = LOWZ_North[2:-1]
LOWZ_North_pos = open("C:\Python34\masters\data\SDSS supercluster\SDSS_DR11\LOWZ_North\Clusters_skypos.txt", "r").read().split('\n')
LOWZ_North_pos = LOWZ_North_pos[1:-1]

for x in LOWZ_North_pos:
    x = x.split('\t')
    cluster.append([x[3], x[1], x[2], x[0]])

##KEEP i AS COUNTER FROM ABOVE
j = 0
while j < len(LOWZ_North):
    z = LOWZ_North[j].split('\t')
    cluster[i].append(z[6]) ##appends radius
    i+=1
    j+=1


##LOWZ South

##Import CMASS_South
LOWZ_South = open("C:\Python34\masters\data\SDSS supercluster\SDSS_DR11\LOWZ_South\Clusters_info.txt", "r").read().split('\n')
LOWZ_South = LOWZ_South[2:-1]
LOWZ_South_pos = open("C:\Python34\masters\data\SDSS supercluster\SDSS_DR11\LOWZ_South\Clusters_skypos.txt", "r").read().split('\n')
LOWZ_South_pos = LOWZ_South_pos[1:-1]

for x in LOWZ_South_pos:
    x = x.split('\t')
    cluster.append([x[3], x[1], x[2], x[0]])

##KEEP i AS COUNTER FROM ABOVE
j = 0
while j < len(LOWZ_South):
    z = LOWZ_South[j].split('\t')
    cluster[i].append(z[6]) ##appends radius
    i+=1
    j+=1


############################################
############################################
############################################
################ISOLATED VOID###############
############################################
############################################
############################################

##Import CMASS_Norths
isovoid_CMASS_North = open("C:\Python34\masters\data\SDSS supercluster\SDSS_DR11\CMASS_North\IsolatedVoids_info.txt", "r").read().split('\n')
isovoid_CMASS_North = isovoid_CMASS_North[2:-1]
isovoid_CMASS_North_pos = open("C:\Python34\masters\data\SDSS supercluster\SDSS_DR11\CMASS_North\IsolatedVoids_skypos.txt", "r").read().split('\n')
isovoid_CMASS_North_pos = isovoid_CMASS_North_pos[1:-1]
void_iso = []

for x in isovoid_CMASS_North_pos:
    x = x.split('\t')
    void_iso.append([x[3], x[1], x[2], x[0]])

i = 0
while i < len(isovoid_CMASS_North):
    z = isovoid_CMASS_North[i].split('\t')
    void_iso[i].append(z[6]) ##appends radius
    i+=1
##Import CMASS_South
isovoid_CMASS_South = open("C:\Python34\masters\data\SDSS supercluster\SDSS_DR11\CMASS_South\IsolatedVoids_info.txt", "r").read().split('\n')
isovoid_CMASS_South = isovoid_CMASS_South[2:-1]
isovoid_CMASS_South_pos = open("C:\Python34\masters\data\SDSS supercluster\SDSS_DR11\CMASS_South\IsolatedVoids_skypos.txt", "r").read().split('\n')
isovoid_CMASS_South_pos = isovoid_CMASS_South_pos[1:-1]


for x in isovoid_CMASS_South_pos:
    x = x.split('\t')
    void_iso.append([x[3], x[1], x[2], x[0]])

##KEEP i AS COUNTER FROM ABOVE
j = 0
while j < len(isovoid_CMASS_South):
    z = isovoid_CMASS_South[j].split('\t')
    void_iso[i].append(z[6]) ##appends radius
    i+=1
    j+=1



##LOWZ North

##Import CMASS_South
isovoid_LOWZ_North = open("C:\Python34\masters\data\SDSS supercluster\SDSS_DR11\LOWZ_North\IsolatedVoids_info.txt", "r").read().split('\n')
isovoid_LOWZ_North = isovoid_LOWZ_North[2:-1]
isovoid_LOWZ_North_pos = open("C:\Python34\masters\data\SDSS supercluster\SDSS_DR11\LOWZ_North\IsolatedVoids_skypos.txt", "r").read().split('\n')
isovoid_LOWZ_North_pos = isovoid_LOWZ_North_pos[1:-1]

for x in isovoid_LOWZ_North_pos:
    x = x.split('\t')
    void_iso.append([x[3], x[1], x[2], x[0]])

##KEEP i AS COUNTER FROM ABOVE
j = 0
while j < len(isovoid_LOWZ_North):
    z = isovoid_LOWZ_North[j].split('\t')
    void_iso[i].append(z[6]) ##appends radius
    i+=1
    j+=1


##LOWZ South

##Import CMASS_South
isovoid_LOWZ_South = open("C:\Python34\masters\data\SDSS supercluster\SDSS_DR11\LOWZ_South\IsolatedVoids_info.txt", "r").read().split('\n')
isovoid_LOWZ_South = isovoid_LOWZ_South[2:-1]
isovoid_LOWZ_South_pos = open("C:\Python34\masters\data\SDSS supercluster\SDSS_DR11\LOWZ_South\IsolatedVoids_skypos.txt", "r").read().split('\n')
isovoid_LOWZ_South_pos = isovoid_LOWZ_South_pos[1:-1]

for x in isovoid_LOWZ_South_pos:
    x = x.split('\t')
    void_iso.append([x[3], x[1], x[2], x[0]])

##KEEP i AS COUNTER FROM ABOVE
j = 0
while j < len(isovoid_LOWZ_South):
    z = isovoid_LOWZ_South[j].split('\t')
    void_iso[i].append(z[6]) ##appends radius
    i+=1
    j+=1



############################################
############################################
############################################
################MINIMAL  VOID###############
############################################
############################################
############################################

##Import CMASS_Norths
minvoid_CMASS_North = open("C:\Python34\masters\data\SDSS supercluster\SDSS_DR11\CMASS_North\MinimalVoids_info.txt", "r").read().split('\n')
minvoid_CMASS_North = minvoid_CMASS_North[2:-1]
minvoid_CMASS_North_pos = open("C:\Python34\masters\data\SDSS supercluster\SDSS_DR11\CMASS_North\MinimalVoids_skypos.txt", "r").read().split('\n')
minvoid_CMASS_North_pos = minvoid_CMASS_North_pos[1:-1]
void_min = []

for x in minvoid_CMASS_North_pos:
    x = x.split('\t')
    void_min.append([x[3], x[1], x[2], x[0]])

i = 0
while i < len(minvoid_CMASS_North):
    z = minvoid_CMASS_North[i].split('\t')
    void_min[i].append(z[6]) ##appends radius
    i+=1
##Import CMASS_South
minvoid_CMASS_South = open("C:\Python34\masters\data\SDSS supercluster\SDSS_DR11\CMASS_South\MinimalVoids_info.txt", "r").read().split('\n')
minvoid_CMASS_South = minvoid_CMASS_South[2:-1]
minvoid_CMASS_South_pos = open("C:\Python34\masters\data\SDSS supercluster\SDSS_DR11\CMASS_South\MinimalVoids_skypos.txt", "r").read().split('\n')
minvoid_CMASS_South_pos = minvoid_CMASS_South_pos[1:-1]


for x in minvoid_CMASS_South_pos:
    x = x.split('\t')
    void_min.append([x[3], x[1], x[2], x[0]])

##KEEP i AS COUNTER FROM ABOVE
j = 0
while j < len(minvoid_CMASS_South):
    z = minvoid_CMASS_South[j].split('\t')
    void_min[i].append(z[6]) ##appends radius
    i+=1
    j+=1



##LOWZ North

##Import CMASS_South
minvoid_LOWZ_North = open("C:\Python34\masters\data\SDSS supercluster\SDSS_DR11\LOWZ_North\MinimalVoids_info.txt", "r").read().split('\n')
minvoid_LOWZ_North = minvoid_LOWZ_North[2:-1]
minvoid_LOWZ_North_pos = open("C:\Python34\masters\data\SDSS supercluster\SDSS_DR11\LOWZ_North\MinimalVoids_skypos.txt", "r").read().split('\n')
minvoid_LOWZ_North_pos = minvoid_LOWZ_North_pos[1:-1]

for x in minvoid_LOWZ_North_pos:
    x = x.split('\t')
    void_min.append([x[3], x[1], x[2], x[0]])

##KEEP i AS COUNTER FROM ABOVE
j = 0
while j < len(minvoid_LOWZ_North):
    z = minvoid_LOWZ_North[j].split('\t')
    void_min[i].append(z[6]) ##appends radius
    i+=1
    j+=1


##LOWZ South

##Import CMASS_South
minvoid_LOWZ_South = open("C:\Python34\masters\data\SDSS supercluster\SDSS_DR11\LOWZ_South\MinimalVoids_info.txt", "r").read().split('\n')
minvoid_LOWZ_South = minvoid_LOWZ_South[2:-1]
minvoid_LOWZ_South_pos = open("C:\Python34\masters\data\SDSS supercluster\SDSS_DR11\LOWZ_South\MinimalVoids_skypos.txt", "r").read().split('\n')
minvoid_LOWZ_South_pos = minvoid_LOWZ_South_pos[1:-1]

for x in minvoid_LOWZ_South_pos:
    x = x.split('\t')
    void_min.append([x[3], x[1], x[2], x[0]])

##KEEP i AS COUNTER FROM ABOVE
j = 0
while j < len(minvoid_LOWZ_South):
    z = minvoid_LOWZ_South[j].split('\t')
    void_min[i].append(z[6]) ##appends radius
    i+=1
    j+=1






i=0
with open("SDSS_sc_reduced.txt", "w") as text_file:
    while i < len(cluster)-1:
        text_file.write("{}\n".format(cluster[i]))
        i+=1
    if i == len(cluster)-1:
        text_file.write("{}".format(cluster[i]))

i=0
with open("SDSS_isovoid_reduced.txt", "w") as text_file:
    while i < len(void_iso)-1:
        text_file.write("{}\n".format(void_iso[i]))
        i+=1
    if i == len(void_iso)-1:
        text_file.write("{}".format(void_iso[i]))

i=0
with open("SDSS_minvoid_reduced.txt", "w") as text_file:
    while i < len(void_min)-1:
        text_file.write("{}\n".format(void_min[i]))
        i+=1
    if i == len(void_min)-1:
        text_file.write("{}".format(void_min[i]))
