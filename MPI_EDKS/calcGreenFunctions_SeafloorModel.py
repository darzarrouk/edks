#!/usr/bin/env python2.6
"""
    by Junle Jiang, May 18, 2013
    Seismological Laboratory
    California Institute of Technology

    Created      : May 18, 2013
    Last modified:

    Modification History:



tested to work with Python 2.7
"""

import numpy as NP
from scipy import interpolate
import cPickle
from ICM.Logger import LogFile
from ICM.Geometry import TriangleCalculator
from ICM.Simplices import Point
from ICM.Simplices import Line
from ICM.Simplices import Triangle
from ICM.CoordinateSystem import PointCoordinates
from projectGFmatrices import projectGFmatrices

import string

def calcGreenFunctions_SeafloorModel(TriPropFile, TriPointsFile, ReceiverFile,\
                                         method_par, plotGeometry):
   """
   method_par has to contain the following info:
      useRecvDir : [False|True]
      Amax : [None|+float]
      EDKSunits : float (units to go from current length units to meters)
      EDKSfilename : name of the EDKS file with the kernels.
      prefix : a name prefix for the output files.
   """


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
   P_IDs = []
   file = open(TriPointsFile, 'r')
   aux = file.readline() # first line is header
   fields = ['id', 'lon', 'lat', 'e',  'n', 'dep', 'No', 'TriIDs'] # uncertain No of triangles
   for line in file:
      line = string.split(line)
      for i in range(1,len(fields)-1):
         line[i] = float( line[i] )
      line[len(fields)-1] = line[len(fields)-1:]
      if len(line) > len(fields):
         line[len(fields):] = []
      Tpoints[line[0]] = dict(zip(fields, line))
      P_IDs.append(line[0])
   file.close()

   # read the receivers.
   Recv = {}
   R_IDs = []
   file = open(ReceiverFile, 'r')
   aux = file.readline() # first line is header
   fields = ['id', 'e', 'n']
   for line in file:
      line = string.split(line)
      for i in range(1, len(fields)):
         line[i] = float( line[i] )
      Recv[line[0]] = dict(zip(fields, line))
      R_IDs.append( line[0] )
   file.close()

   # save log
   log.addLine('data successfully loaded ... ')

   # assemble the array with the east and north coordinates of the observation points
   eR = []
   nR = []
   for Rid in R_IDs:
       eR.append(Recv[Rid]['e'])
       nR.append(Recv[Rid]['n'])


   msg = '%i receivers to compute displacements...' %(len(eR))
   log.addLine(msg)

   # now compute the coordinates of all the subtriangles
   log.addLine(' Subdividing triangular Mesh ')

   if plotGeometry:
      MasterTriangles4Plot = []


   eST = []
   nST = []
   dST = []
   aST = []
   wST = []
   idST = []
   strikeST = []
   dipST = []
   slipST = []

   for Pid in P_IDs:
      No_nghbr = Tpoints[Pid]['No']
      # print No_nghbr
      T_nghbr = Tpoints[Pid]['TriIDs']
      # calculate for each neighboring triangle
      for Tid in T_nghbr[0:int(No_nghbr)]:
         # initialize coordinate system for triangles
         TriCS = PointCoordinates(CoordSystemName = 'TriangleCS',\
                            xyzNames = ['E', 'N', 'U'])
         idP1 = Tprop[Tid]['idP1']
         idP2 = Tprop[Tid]['idP2']
         idP3 = Tprop[Tid]['idP3']
         TriCS.addPoint(Tpoints[idP1]['e'], Tpoints[idP1]['n'], -1.0*Tpoints[idP1]['dep'],\
                      idP1)
         TriCS.addPoint(Tpoints[idP2]['e'], Tpoints[idP2]['n'], -1.0*Tpoints[idP2]['dep'],\
                      idP2)
         TriCS.addPoint(Tpoints[idP3]['e'], Tpoints[idP3]['n'], -1.0*Tpoints[idP3]['dep'],\
                      idP3)
         # instantiate Triangle calculator
         Tcalc = TriangleCalculator(TriCS)

         # instantiate the Points
         P1 = Point(idP1)
         P2 = Point(idP2)
         P3 = Point(idP3)
         # instantiate the lines which form the sides of the triangle
         L12 = Line(P1, P2, 'L12')
         L23 = Line(P2, P3, 'L23')
         L31 = Line(P3, P1, 'L31')
         # instantiate the triangle
         Tri = Triangle(L12, L23, L31, Tid)
         # Check for consistent orientation of the triangles (normal must point up)
         # if the normal of the triangle is pointing down, i reverse its orientation.
         upDir = NP.array([0, 0, 1])
         n = Tcalc.unitNormal(Tri)
         if (NP.dot(upDir, n) < 0):
            Tri.reverseOrientation()
         # compute the subfaults for the triangle
         # get the maximum Area for the triangle
         if Amax == None:
            print 'Have to specify Amax!'
         else:
            TriAmax = Amax
         Ce, Cn, Cz, A  = getReducedTrianglesProp(Tri, TriCS, TriAmax)
         eST.extend(Ce)
         nST.extend(Cn)
         dST.extend(-1.0 * Cz)
         aST.extend(A)
         # record the id for point instead of triangle
         idST.extend( [Pid] * len(A) )

         # get angles
         # phi, phiDir = Tcalc.strike(Tri)
         # strikeST.extend(phi * NP.ones(len(A)))
         # dip, dipDir = Tcalc.dip(Tri)
         # dipST.extend(dip * NP.ones(len(A)))

         ## assign weighted slip
         # slipST.extend(1.0 * NP.ones(len(A)))
         Cw = getWeight4SubTriangles(Tri, TriCS, Pid, Ce, Cn, Cz)
         slipST.extend(1.0 * NP.array(Cw))

         if plotGeometry:
            MasterTriangles4Plot.append([Tri, TriCS])

         msg = 'Triangle %s has %i subsources...'%(Tid, len(A))
         log.addLine(msg)


   # Convert units
   Units2meters = method_par['EDKSunits']
   eST = NP.array(eST) * Units2meters
   nST = NP.array(nST) * Units2meters
   dST = NP.array(dST) * Units2meters
   aST = NP.array(aST) * Units2meters * Units2meters
   eR = NP.array(eR) * Units2meters
   nR = NP.array(nR) * Units2meters
   # NP.savetxt('test-eST.txt',eST)
   # NP.savetxt('test-eR.txt',eR)
   # return

   if plotGeometry:
      import pylab as PL
      import mpl_toolkits.mplot3d.axes3d as plot3
      fig = PL.figure(1)
      s = plot3.Axes3D(fig)
      s.plot3D(eST/Units2meters, nST/Units2meters, -1.0*dST/Units2meters, '.r',\
		markersize = 2)
      for MasterTriangle in MasterTriangles4Plot:
         Mtri, MtriCS = MasterTriangle
         plotTriangle(Mtri, MtriCS, s, color = 'k')

      s.plot(1.0 * eR/Units2meters, 1.0 * nR/Units2meters, 'og')
      PL.show()

   #PL.plot(-nST, eST, '.k', markersize = 2)
   #PL.plot(-nR, eR, 'or')
   #PL.show()

   # the rest just convert to arrays.
   # strikeST = NP.array(strikeST)
   # dipST = NP.array(dipST)
   slipST = NP.array(slipST)

   # run layered_disloc_sub on all the subtriangles
   edks = method_par['EDKSfilename']
   prefix = method_par['prefix']

   # interp triangle deformation to observation points
   log.addLine('calculating Seafloor Green functions...')
   # log.addline('Shape: ' + str(idST.shape) + '; Rev:' + str(eR.shape))
   Gu = interpTri2Rec(idST, eST, nST, slipST, aST, eR, nR)

   # save the GF matrices.

   # write the matrices
   log.addLine('writing output matrix...')
   # Gu
   file = open('Gu.txt','w')
   Nrows, Ncols = Gu.shape
   print 'Nrows = ' + str(Nrows) + ', Ncols = ' + str(Ncols)
   for row in range(0,Nrows):
      for col in range(0,Ncols):
         value = '%s ' %(Gu[row][col])
         file.write(value)
      file.write('\n')
   file.close()

   return

###
def interpTri2Rec(IDs, xs, ys, slip, A, xr, yr):
   """
   interpolate from tent-like patches to seafloor receiver grid
   """
   log0 = LogFile('log.interpTri2Rec.txt')
   log0.clearFile()

   nrec = len(xr)
   np = len(set(IDs))
   ntsp = len(xs) # total number of point sources (sub sources)

   # compute the mapping between the string IDs and non decreasing positive integer number
   setOfAlreadyStoredIDs = set() # this is just for testing the right order of the IDs.
   sortedListOfFiniteFaultIDs = []
   ident = []
   i = 1
   IDprev = IDs[0]
   ident.append(i)
   sortedListOfFiniteFaultIDs.append(IDprev)
   setOfAlreadyStoredIDs.add(IDprev)
   NumSubSources = {}
   NumSubSources[IDprev] = 1
   u = NP.zeros((nrec, np))

   log0.addLine(" The Shape of Gu.txt: " + str(u.shape))

   for k in range(1, ntsp):
      if IDs[k] == IDprev:
         ident.append(i)
         NumSubSources[IDprev] += 1

      else: # the current ID is a new one.
         if IDs[k] in setOfAlreadyStoredIDs: # this is an error
            raise ValueError('The source IDs are not in the right order, see the help')
         else:
            IDprev = IDs[k]
            i += 1
            ident.append(i)
            sortedListOfFiniteFaultIDs.append(IDs[k])
            setOfAlreadyStoredIDs.add(IDs[k])
            NumSubSources[IDs[k]] = 1

   from matplotlib.nxutils import points_inside_poly

   k = 0
   for ind, j in enumerate(sortedListOfFiniteFaultIDs):
      log0.addLine(' Working on patch ' + str(ind) + ': ' + str(j))
      nSubSrc = NumSubSources[j]
      Xs = NP.zeros(nSubSrc)
      Ys = NP.zeros(nSubSrc)
      us = NP.zeros(nSubSrc)
      ur = NP.zeros(len(xr))
      for m in range(0, nSubSrc):
         Xs[m] = xs[k+m]
         Ys[m] = ys[k+m]
         us[m] = slip[k+m]

      # print Us
      margin = 1000
      xmin = Xs.min() - margin
      xmax = Xs.max() + margin
      ymin = Ys.min() - margin
      ymax = Ys.max() + margin

      # Bounding regions for the source region
      log0.addLine(' Before mask')
      mask = (xr > xmin) & (xr < xmax) & (yr > ymin) & (yr < ymax)
      log0.addLine(' After mask ' + str(sum(mask)))

      vert = convex_hull(zip(Xs, Ys))
      log0.addLine(' Finding vertices for convex hull:' + str(vert))

      # Bounding convex polygons for the receiver region
      gmask = points_inside_poly(NP.vstack((xr[mask], yr[mask])).T, vert)
      log0.addLine(' Finding points inside polygon:' + str(gmask))

      # if ind==2:
      #   NP.savetxt('test-us.txt', us)
      #   NP.savetxt('test-ur.txt', NP.vstack((xr[mask], yr[mask])).T)
      # add in corners before interpolation (scipy 0.9)
      # corner = NP.zeros((4,3))
      # corner[0,:] = [xmin, ymin, 0]
      # corner[1,:] = [xmin, ymax, 0]
      # corner[2,:] = [xmax, ymin, 0]
      # corner[3,:] = [xmax, ymax, 0]
      # us = NP.vstack((us, corner))

      log0.addLine(' Before interpolate, number of sub-sources: ' + str(NumSubSources[j]))
      ur0 = ur[mask][gmask]
      xr0 = xr[mask][gmask]
      yr0 = yr[mask][gmask]

      # ur0[gmask] = interpolate.griddata((us[:,0], us[:,1]), us[:,2], (xr[mask][gmask], yr[mask][gmask]),
      #                                 method='cubic', fill_value=0)
      ur0 = interpolate.griddata((Xs, Ys), us, (xr0, yr0), method='cubic', fill_value=0)

      ur[mask][gmask] = ur0
      log0.addLine(' After interpolate: ' + str(ur0))
      u[:, ind] = ur
      k = k + NumSubSources[j]

   # NP.savetxt('Output_sortedIDs.txt',sortedListOfFiniteFaultIDs, fmt='%s')
   # NP.savetxt('Output_NumSubSources.txt',NumSubSources.items(), fmt='%s')

   return u

###
def convex_hull(points):
    """Computes the convex hull of a set of 2D points.

    Input: an iterable sequence of (x, y) pairs representing the points.
    Output: a list of vertices of the convex hull in counter-clockwise order,
      starting from the vertex with the lexicographically smallest coordinates.
    Implements Andrew's monotone chain algorithm. O(n log n) complexity.
    """

    # Sort the points lexicographically (tuples are compared lexicographically).
    # Remove duplicates to detect the case we have just one unique point.
    points = sorted(set(points))

    # Boring case: no points or a single point, possibly repeated multiple times.
    if len(points) <= 1:
        return points

    # 2D cross product of OA and OB vectors, i.e. z-component of their 3D cross product.
    # Returns a positive value, if OAB makes a counter-clockwise turn,
    # negative for clockwise turn, and zero if the points are collinear.
    def cross(o, a, b):
        return (a[0] - o[0]) * (b[1] - o[1]) - (a[1] - o[1]) * (b[0] - o[0])

    # Build lower hull
    lower = []
    for p in points:
        while len(lower) >= 2 and cross(lower[-2], lower[-1], p) <= 0:
            lower.pop()
        lower.append(p)

    # Build upper hull
    upper = []
    for p in reversed(points):
        while len(upper) >= 2 and cross(upper[-2], upper[-1], p) <= 0:
            upper.pop()
        upper.append(p)

    # Concatenation of the lower and upper hulls gives the convex hull.
    # Last point of each list is omitted because it is repeated at the beginning of the other list.
    return lower[:-1] + upper[:-1]


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
   Area = Tcalc.Area(Tri)
   # get subdivideFlag for the triangle
   if Amax != None:
      subdivideFlag = Amax < Area # True if Amax < Area => it gets subdivided
   else: # Amax == None
      # get the Triangle's point coordinates
      Pids = Tri.getPointsID()
      Pcoords = [TriCS.getPointCoord(Pid) for Pid in Pids]
      # get the minimum depth
      Pdepths = [NP.abs(Pcoord[2]) for Pcoord in Pcoords] # Z is positive up.
      depT = NP.min(Pdepths)
      # the caracteristic length of the Maximum triangle from Saint-venant's principle
      cLT = 1.0 * depT / depT2cLTratio
      DiC = NP.sqrt(3.0) * cLT / 3.0 # diameter of the inscribed circle in equi Tri
      Amax = NP.pi * DiC * DiC / 4.0  #  The Maximum area of the triangle is the area of
                                      # the inscibed circle in an equilateral triangle
                                      # with side equal to cLT
      # calculate the subdivideFlag
      subdivideFlag = Amax < Area

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
      return [NP.array(Ce), NP.array(Cn), NP.array(Cz), NP.array(A)]
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
def getWeight4SubTriangles(Tri, TriCS, Pid, Ce, Cn, Cz):
   """
    calculate the weight for each subtriangles given its location in the tent element
   """
   # get the Point instances and coordinates of the master triangle
   PtList = Tri.getPointsID()
   print Pid
   print PtList
   e1, n1, u1 = TriCS.getPointCoord(Pid)
   e2, n2, u2 = TriCS.getPointCoord(list(set(PtList)-set([Pid]))[0])
   e3, n3, u3 = TriCS.getPointCoord(list(set(PtList)-set([Pid]))[1])

   Cw = []
   for i in range(0, len(Ce)):
       ep = Ce[i]
       np = Cn[i]
       zp = Cz[i]
       L23 = ((e3-e2)**2 + (n3-n2)**2)**0.5
       distp = ((np-n2)*(e3-e2) - (ep-e2)*(n3-n2))/L23
       dist1 = ((n1-n2)*(e3-e2) - (e1-e2)*(n3-n2))/L23
       Cw.append(distp/dist1)
   return Cw


###
def calcMidPoint(ei, ni, zi, ej, nj, zj):
   em = 0.5*(ei + ej)
   nm = 0.5*(ni + nj)
   zm = 0.5*(zi + zj)
   return [em, nm, zm]


###
def plotTriangle(T, TCS, s, color = 'k'):
   # plot the triangle
   P1, P2, P3 = T.getPoints()
   x1, y1, u1 = TCS.getPointCoord(P1.id)
   x2, y2, u2 = TCS.getPointCoord(P2.id)
   x3, y3, u3 = TCS.getPointCoord(P3.id)

   x = [x1, x2, x3, x1]
   y = [y1, y2, y3, y1]
   u = [u1, u2, u3, u1]
   s.plot(x, y, u, color)



if __name__ == '__main__':

   from EDKSsubParams import *
   parNames = ['useRecvDir', 'Amax', 'EDKSunits', 'EDKSfilename', 'prefix']
   parValues = [ useRecvDir ,  Amax ,  EDKSunits ,  EDKSfilename ,  prefix ]
   method_par = dict(zip(parNames, parValues))

   calcGreenFunctions_SeafloorModel(TriPropFile, TriPointsFile, ReceiverFile,\
                                         method_par, plotGeometry)


