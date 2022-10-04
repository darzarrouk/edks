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

# parameters to map projection (this is true for Maule EQ surface)
lon0 = 143.50
lat0 = 41.50
ProjPar = {'Type' : 'transversemercator'}
ProjPar['Origin'] = {'lon0' : lon0 , 'lat0': lat0}
MapProj = ICM.CoordTransformations.MapProjector(ProjPar)

# read the mesh file
id = []
lon = []
lat = []

file = open('Tohoku_PostSeismic_2011_10_25.dat', 'r')
# first line is header
line = file.readline()
print line
for line in file:
   l = string.split(line)
   id.append( l[0] )
   lon.append( float(l[1]) )
   lat.append( float(l[2]) )
   #line = '%s % -7.5f % -7.5f' %(id[-1], lon[-1] , lat[-1])
   #print line
file.close()

# Project the mesh points
lon = NP.array(lon)
lat = NP.array(lat)
E, N = MapProj.projectXY(lon, lat, inverse = False)


# compute the output file
file = open('JapanGPSObs.idEN', 'w')
# write heather line
label = 'obsId E N\n'
file.write(label)
# write the remaining data
for i in range(0,len(id)):
   # I am only using the data for the spatial window:
   # 34 <= lat <= 43, 138 <= lon <= 144
   #if lat[i] >= 34 and lat[i] <= 45 and lon[i] >= 138 and lon[i] <= 147: 
      line = '%s ' %(id[i])
      line += '% -7.5f % -7.5f\n' %(E[i], N[i])
      file.write(line)

file.close()


