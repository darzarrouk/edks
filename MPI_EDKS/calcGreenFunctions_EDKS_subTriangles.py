"""
    by Francisco Hernan Ortega Culaciati, Nov 09, 2011
    Seismological Laboratory
    California Institute of Technology
    
    Created      : Nov 09, 2011
    Last modified: Nov 09, 2011 

    Modification History:
    


tested to work with Python 2.7
"""

import numpy as NP
import cPickle
from ICM.Logger import LogFile
from ICM.Geometry import TriangleCalculator
from ICM.Simplices import Point
from ICM.Simplices import Line
from ICM.Simplices import Triangle
import string

def main(TriPropFile, TriPointsFile, ReceiverFile, useRecvDir = False, Amax = None):
   log = LogFile('log.calcEDKSsubGreenFunctions.txt')
   log.clearFile()
   log.addLine('loading data...')
   # read TriPropFile
   Tprop = {}
   T_IDs = [] # I save the order in the file so I can produce GF in same order
   file = open(TriPropFile, 'r')
   aux = file.readline() # first line is header
   fields = ['id', 'lon', 'lat', 'e', 'n', 'dep', 'strike', 'dip', 'area',\
   'idP1', 'idP2', 'idP3' ]
   for line in file:
      line = string.split(line)
      for i in range(1,9): # all the inner ones are float
         line[i] = float( line[i] )
      Tprop[line[0]] = dict(zip(fields, line))
      T_IDs.append(line[0])
   file.close()
   
   # read the TriPoints file
   Tpoints = {}
   file = open(TriPointsFile, 'r')
   aux = file.readline() # first line is header
   fields = ['id', 'lon', 'lat', 'e',  'n', 'dep']
   for line in file:   
      line = string.split(line) 
      for i in range(1:len(fields)):
         line[i] = float( line[i] )
      Tpoints[line[0]] = dict(zip(fields, line))
   file.close()

   # read the receivers.
   Recv = {}
   R_IDs = []
   file = open(ReceiverFile, 'r')
   aux = file.readline() # first line is header
   if useRecvDir:
      fields = ['id', 'e', 'n', 'ODirE', 'ODirN', 'ODirU']
   else:
      fields = ['id', 'e', 'n'] 
   for line in file:
      line = string.split(line)
      for i in range(1:len(fields)):
         line[i] = float( line[i] )
      Recv[line[0]] = dict(zip(fields, line))
      R_IDs.append( line[0] )
   file.close()
      
   # save log
   log.addLine('data successfully loaded from file ' + PickleDataFile)

   # now compute the coordinates of all the subtriangles
   # HERE

   # Compute TriDict and check for normal orientation
   TriDict = {}
   Tcalc = TriangleCalculator(TriCS)
   upDir = NP.array([0, 0, 1])

   for Tid in Tcon.keys():
      # instantiate the Points
      P1 = Point(Tcon[Tid][0])
      P2 = Point(Tcon[Tid][1])
      P3 = Point(Tcon[Tid][2])
      # instantiate the lines which form the sides of the triangle
      L12 = Line(P1, P2, 'L12')
      L23 = Line(P2, P3, 'L23')
      L31 = Line(P3, P1, 'L31')
      # instantiate the triangle
      TriDict[Tid] = Triangle(L12, L23, L31, Tid)

      # Check for consistent orientation of the triangles (normal must point up)
      # if the normal of the triangle is pointing down, i reverse its orientation.
      n = Tcalc.unitNormal(TriDict[Tid])
      if (NP.dot(upDir, n) < 0):
         TriDict[Tid].reverseOrientation()

   log.addLine(' Subdividing triangular Mesh ')
   # Now compute the location of the subfaults for all the Triangles
   # Subfaults locations are computed by subdividing triangles into similar ones, 
   # until stop criterion is met. The location is the centroid of the subTriangles.
   # the id of the subtriangles is the same as the one for the Master triangle.
   T_Ce = []
   T_Cn = []
   T_Cz = []
   T_A = []
   T_IDs = []
   for Tid in TriDict.keys():
      Tri = TriDict[Tid]
      Ce, Cn, Cz, A = getReducedTrianglesProp(Tri, TriCS, Amax = Amax)
      T_Ce.extend( T_Ce )
      T_Cn.extend( T_Cn )
      T_Cz.extend( T_Cz )
      T_A.extend( A )
      T_IDs.extend( [Tid] * len(A) ) # all IDs of subfaults are the same fot Tid

   msg = 'Triangle %s has %i subsources...'%(Tid, len(T_A))
   log.addLine(msg)
   
   # assemble the array with the east and north coordinates of the observation points
   obsID = ObsDict.keys()

   e = []
   n = []
   for key in obsID:
       Pcoord = ObsDict[key].getPoint(ObsCS)
       e.append(Pcoord[0])
       n.append(Pcoord[1])

   e = NP.array(e)
   n = NP.array(n)


   # Convert units
   Units2meters = method_par['EDKSunits']
   T_Ce = NP.array(T_Ce) * Units2meters
   T_Cn = NP.array(T_Cn) * Units2meters
   T_Cdep = -1.0 * NP.array(T_Cz) * Units2meters
   T_A = NP.array(T_A) * Units2meters * Units2meters

   





   edks = method_par['EDKSfilename']
   

   GF = EDKSsub_Tri_GF(TriDict, TriCS, ObsDict, ObsCS, Amax,\
                                              edks, Units2meters)      
   else:
      msg = ' Invalid method %s ...' %(method)
      log.addLine(msg)


   


###
def getReducedTrianglesProp(Tri, TriCS, Amax = None, depT2cLTratio = 4.0):
   """
   get the Triangle and its coordinate system, recursively subdivide it until 
   all the subdivided triangles have an area less or equal Amax. 
   Each triangle is subdivided into 4 triangles by adding the center points on
   the edge of the master triangle.
   
   This function returns 3 Numpy arrays with the coordinates of the centers of the
   subdivided triangles (east, north and depth) and its Area 

   if Amax is not given, it will be calculated based on the Saint-Venant principle using 
   as a reference the depth of the shallowest vertex of the triangle(depT). Thus, it will
   divide the triangle until the characteristic length of the triangle (cLT) follows:
   
                            depT >= depT2cLTratio * cLT

   the default value of depT2cLTratio is 4 ( to be conservative within the Saint-Venant
   principle (this assumes there are receivers at horizontal distance = 0.0, so the
   minimum distance from source to receiver is the depth ).

   if Amax is given, it will calculate the area of the triangle and subdivide until
   AreaTriangle <= Amax.

   """
   # calculate the Area of the triangle
   Tcalc = TriangleCalculator(TriCS)
   # get subdivideFlag for the triangle
   if Amax != None:
      Area = Tcalc.Area(Tri)
      subdivideFlag = Amax < Area # True if Amax < Area => it gets subdivided
   else: # Amax == None
      # get the Triangle's point coordinates
      Pids = Tri.getPointsID()
      Pcoords = [TriCS.getPointCoord(Pid) for Pid in Pids]
      # get the minimum depth
      Pdepths = [NP.abs(Pcoord[2]) for Pcoord in Pcoords] # Z is positive up.
      depT = NP.min(Pdepths)
      # geth the length os the longest side of the triangle
      Lside = []
      for i in [0,1,2]:
         Pi = Pcoords[i]
         j = NP.mod(i+1,3)
         Pj = Pcoords[j]
         Lij = NP.sqrt( NP.sum( (Pj - Pi)*(Pj - Pi) ) )
         Lside.append( Lij )
      cLT = NP.max( Lside )
      # calculate the subdivideFlag
      subdivideFlag =  depT2cLTratio * cLT > depT
      
   # to store the triangles center and area
   Ce = []
   Cn = []
   Cz = []
   A = []
   # decide if we split or not.
   if not(subdivideFlag): # we calculate the properties and return it
      Area = Tcalc.Area(Tri)
      Centroid = Tcalc.Centroid(Tri)
      Ce.append(Centroid[0])
      Cn.append(Centroid[1])
      Cz.append(Centroid[2])# Depth is positive downwards,U is positive upwards.
      A.append(Area)
      return [Ce, Cn, Cz, A]
   else: # we split into 4 child triangles and work recursively
      T1, T2, T3, T4, TCS = splitTriangle(Tri, TriCS)
      Ce1, Cn1 , Cz1, A1 = getReducedTrianglesProp(T1, TCS, Amax)
      Ce2, Cn2 , Cz2, A2 = getReducedTrianglesProp(T2, TCS, Amax)
      Ce3, Cn3 , Cz3, A3 = getReducedTrianglesProp(T3, TCS, Amax)
      Ce4, Cn4 , Cz4, A4 = getReducedTrianglesProp(T4, TCS, Amax)
      Ce.extend(Ce1)
      Ce.extend(Ce2)
      Ce.extend(Ce3)
      Ce.extend(Ce4)
      Cn.extend(Cn1)
      Cn.extend(Cn2)
      Cn.extend(Cn3)
      Cn.extend(Cn4)
      Cz.extend(Cz1)
      Cz.extend(Cz2)
      Cz.extend(Cz3)
      Cz.extend(Cz4)
      A.extend(A1)
      A.extend(A2)
      A.extend(A3)
      A.extend(A4)

   return [NP.array(Ce), NP.array(Cn), NP.array(Cz), NP.array(A)]


###   
def splitTriangle(Tri, TriCS):
   """
   takes the triangle Tri and split it into 4 triangles by adding the midpoints 
   of each side of the triangle.

   """
   # get the Point instances and coordinates of the master triangle
   P1, P2, P3 = Tri.getPoints()
   e1, n1, u1 = TriCS.getPointCoord(P1.id)
   e2, n2, u2 = TriCS.getPointCoord(P2.id)
   e3, n3, u3 = TriCS.getPointCoord(P3.id)

   # compute the midpoint coordinates for each side of the triangle
   e12, n12, u12 = calcMidPoint(e1, n1, u1, e2, n2, u2)
   e23, n23, u23 = calcMidPoint(e2, n2, u2, e3, n3, u3)
   e31, n31, u31 = calcMidPoint(e3, n3, u3, e1, n1, u1)

   # instantiate points and midpoints for sides of master triangle
   idP1 = 'P1'
   idP12 = 'P12'
   idP2 = 'P2'
   idP23 = 'P23'
   idP3 = 'P3'
   idP31 = 'P31'
   P1 = Point(idP1)
   P12 = Point(idP12)
   P2 = Point(idP2)
   P23 = Point(idP23)
   P3 = Point(idP3)
   P31 = Point(idP31)

   # instantiate new coordinate system
   TCS = PointCoordinates(CoordSystemName = 'LocalTriangleCoordinateSystem',\
                         xyzNames = ['E', 'N', 'U'])
   TCS.addPoint(e1, n1, u1, idP1)
   TCS.addPoint(e12, n12, u12, idP12)
   TCS.addPoint(e2, n2, u2, idP2)
   TCS.addPoint(e23, n23, u23, idP23)
   TCS.addPoint(e3, n3, u3, idP3)
   TCS.addPoint(e31, n31, u31, idP31)

   # instantiate the 4 child triangles T1, T2, T3, T4
   # T1
   L1 = Line(P1, P12, 'L1')
   L2 = Line(P12, P31, 'L2')
   L3 = Line(P31, P1, 'L3')
   T1 = Triangle(L1, L2, L3, 'T1')
   # T2
   L1 = Line(P12, P2, 'L1')
   L2 = Line(P2, P23, 'L2')
   L3 = Line(P23, P12, 'L3')
   T2 = Triangle(L1, L2, L3, 'T2')
   # T3
   L1 = Line(P31, P23, 'L1')
   L2 = Line(P23, P3, 'L2')
   L3 = Line(P3, P31, 'L3')
   T3 = Triangle(L1, L2, L3, 'T3')
   # T4
   L1 = Line(P12, P23, 'L1')
   L2 = Line(P23, P31, 'L2')
   L3 = Line(P31, P12, 'L3')
   T4 = Triangle(L1, L2, L3, 'T4')

   return [T1, T2, T3, T4, TCS]

###
def calcMidPoint(ei, ni, zi, ej, nj, zj):
   em = 0.5*(ei + ej)
   nm = 0.5*(ni + nj)
   zm = 0.5*(zi + zj)
   return [em, nm, zm]


