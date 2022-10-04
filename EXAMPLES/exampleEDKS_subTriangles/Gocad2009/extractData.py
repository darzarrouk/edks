#!/usr/bin/env python
import numpy as NP
import ICM
import sys


GeomodProjName = 'TohokuOki'

TSfilename = sys.argv[1]
filenameOutTriangles = TSfilename + '.TriangleProp.txt'
filenameOutPoints = TSfilename + '.PointCoord.txt'
filenameFigureTriangles = TSfilename + '.figureTriangles.ps'
# instantiate the Map projector
lon0 = 143.50
lat0 = 41.50
ProjPar = {'Type' : 'transversemercator'}
ProjPar['Origin'] = {'lon0' : lon0 , 'lat0': lat0}
MapProj = ICM.CoordTransformations.MapProjector(ProjPar)


# read the Gocad TSurf (get Triangle connectivity and Triangle Coordinate System)
Zdown = True  # Z is positive upwards in this Gocad Project
Zscale = 1.0/1000.0
TSurf = ICM.Gocad2009.TSurfReader(TSfilename, Zdown, Zscale)
P = TSurf['P'] # Dict with Point instances
T = TSurf['T'] # Dict with Triangle Instances
CS = TSurf['CS'] # Coordinate System Instance
Tcalc = ICM.Geometry.TriangleCalculator(CS) 

# check that the normal of the triangle is pointing upwards
# up direction is [0,0,1]
upDir = NP.array([0,0,1])
for id in T.keys():
   n = Tcalc.unitNormal(T[id])
   if (NP.dot(upDir, n) < 0):
            T[id].reverseOrientation()

# instantiate the CoordinateTransformation (from E,N to Lon,Lat)


# calculate the Geometric properties for each triangle
a = []
s = []
d = []
Ec = []
Nc = []
Dc = []
Tid = []
idP1 = []
idP2 = []
idP3 = []
TriangleIDs = T.keys()
TriangleIDs.sort()

for id in TriangleIDs:
   Tid.append(id)
   a.append( Tcalc.Area(T[id]) )
   strike, strikeDir = Tcalc.strike(T[id])
   s.append( strike )
   dip, dipDir = Tcalc.dip(T[id])
   d.append( dip )
   Centroid = Tcalc.Centroid(T[id])
   Ec.append( Centroid[0] )
   Nc.append( Centroid[1] )
   Dc.append( -1.0 * Centroid[2] )
   P1, P2, P3 = T[id].getPoints()
   idP1.append(P1.id)
   idP2.append(P2.id)
   idP3.append(P3.id)

# compute points id set (without repeated values)
Pid = []
Pid.extend(idP1)
Pid.extend(idP2)
Pid.extend(idP3)
Pid = list(set(Pid))
Pid.sort()
# compute the coordinates of the points
Ep = []
Np = []
Dp = []
for i in range(0,len(Pid)):
   e, n, z = CS.getPointCoord(Pid[i])
   Ep.append(e)
   Np.append(n)
   Dp.append(-z)

   

# project cartesian positions back to Lon Lat
Ec = NP.array(Ec)
Nc = NP.array(Nc)
LONc, LATc = MapProj.projectXY(Ec, Nc, inverse = True)

Ep = NP.array(Ep)
Np = NP.array(Np)
LONp, LATp = MapProj.projectXY(Ep, Np, inverse = True)

# print the results to screen and file
f = open(filenameOutTriangles, 'w')
line = 'Tid lon lat E[km] N[km] dep[km] strike dip Area[km2] idP1 idP2 idP3\n'
print(line)
f.write(line)
for i in range(0, len(Ec)):
   line =  '%s % -7.02f % -7.02f ' %(Tid[i], LONc[i], LATc[i])
   line +=  '% -7.03f % -7.03f % -7.03f ' %(Ec[i], Nc[i], Dc[i])
   line += '% -7.03f % -7.03f % -10.03f ' %(s[i], d[i], a[i])
   line += '%s %s %s' %(idP1[i], idP2[i], idP3[i]) 
   print line
   line += '\n'
   f.write(line)

f.close()


f = open(filenameOutPoints, 'w')
line = 'Pid lon lat E[km] N[km] dep[km]\n'
print line
f.write(line)
for i in range(0, len(Pid)):
   line = '%s % -7.02f % -7.02f ' %(Pid[i], LONp[i], LATp[i])
   line += '% -7.03f % -7.03f % -7.03f' %(Ep[i], Np[i], Dp[i]) 
   print line
   line += '\n'
   f.write(line)

f.close()


# make some plots
import pylab as PL
fig = PL.figure(1)
s1 = PL.subplot(1,1,1)
s1.plot(Ec, Nc, 'r.')
s1.plot(Ep, Np, 'k.')
s1.axis('equal')
# plot the triangular mesh
for id in TriangleIDs:
   P1, P2, P3 = T[id].getPoints() 
   e1, n1, z1 = CS.getPointCoord(P1.id)
   e2, n2, z2 = CS.getPointCoord(P2.id)
   e3, n3, z3 = CS.getPointCoord(P3.id)
   s1.plot([e1, e2, e3, e1], [n1, n2, n3, n1], 'k', linewidth = 1)
fig.savefig(filenameFigureTriangles)


