#!/usr/bin/env python
"""
convert the mesh that Junle gave me

id lon lat

to the file needed by ICM package to calculate the GF

obsId obs stdev ODirE ODirN ODirU E N U CODE lon lat

"""
import ICM
import numpy as NP
import string

# parameters to map projection (this is true for Sendai EQ surface)
lon0 = 143.50
lat0 = 41.50
ProjPar = {'Type' : 'transversemercator'}
ProjPar['Origin'] = {'lon0' : lon0 , 'lat0': lat0}
MapProj = ICM.CoordTransformations.MapProjector(ProjPar)

# read the mesh file
id = []
lon = []
lat = []
DDirE = []
DDirN = []
DDirU = []
file = open('ITO_idLonLat.dat', 'r')
# first line is header
line = file.readline()
print line
for line in file:
   l = string.split(line)
   id.append( l[0] )
   lon.append( float(l[1]) )
   lat.append( float(l[2]) )
   DDirE.append( float(l[5]) )
   DDirN.append( float(l[6]) )
   DDirU.append( float(l[7]) )

   #line = '%s % -7.5f % -7.5f' %(id[-1], lon[-1] , lat[-1])
   #print line
file.close()

# Project the data points
lon = NP.array(lon)
lat = NP.array(lat)
E, N = MapProj.projectXY(lon, lat, inverse = False)


# compute the output file
file = open('ITO.EDKSsub', 'w')
# write heather line
label = 'obsId E N DDirE DDirN DDirU\n'
file.write(label)
# write the remaining data
for i in range(0,len(id)):
   line = '%s ' %(id[i])
   line += '% -7.5f % -7.5f ' %(E[i], N[i])
   line += '% -7.5f % -7.5f % -7.5f \n' %(DDirE[i], DDirN[i], DDirU[i])
   file.write(line)

file.close()


