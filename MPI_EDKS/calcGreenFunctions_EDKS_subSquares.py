#!/usr/bin/env python
"""
    by Francisco Hernan Ortega Culaciati, Nov 16, 2011
    Seismological Laboratory
    California Institute of Technology
    
    Created      : Nov 16, 2011
    Last modified: Nov 16, 2011 

    Modification History:
    


tested to work with Python 2.7
"""

import numpy as NP
import cPickle
from ICM.Logger import LogFile
from ICM.Geometry import RectangleCalculator
from ICM.Rectangles import rectanglePatch
from layered_disloc_sub import layered_disloc_sub
from projectGFmatrices import projectGFmatrices

import string

def calcGreenFunctions_EDKS_subSquares(SquaresPropFile, ReceiverFile,\
                                         method_par, plotGeometry, w_ascii=True, w_bin=True):
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
   # read SquaresPropFile
   Rprop = {}
   Rec_IDs = [] # I save the order in the file so I can produce GF in same order
   file = open(SquaresPropFile, 'r')
   aux = file.readline() # first line is header
   fields = ['e', 'n', 'dep', 'strike', 'dip', 'area', 'id']
   for line in file:
      line = string.split(line)
      for i in range(0,6): # all the left ones are float
         line[i] = float( line[i] )
      Rprop[line[-1]] = dict(zip(fields, line))
      Rec_IDs.append(line[-1])
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
   if useRecvDir:
      ODirE = []
      ODirN = []
      ODirU = []
   for Rid in R_IDs: 
       eR.append(Recv[Rid]['e'])
       nR.append(Recv[Rid]['n'])
       if useRecvDir:
          ODirE.append(Recv[Rid]['ODirE'])
          ODirN.append(Recv[Rid]['ODirN'])
          ODirU.append(Recv[Rid]['ODirU'])

   msg = '%i receivers to compute displacements...' %(len(eR))
   log.addLine(msg)

   # now compute the coordinates of all the sub rectangles
   log.addLine(' Subdividing 3Dsquare Mesh ')
   # instantiate rectangle calculator
   Rcalc = RectangleCalculator()

   if plotGeometry:
      subRectangles4Plot = []
      Rectangles4Plot = []
 
   eSR = []
   nSR = []
   dSR = []
   aSR = []
   idSR = []
   strikeSR = []
   dipSR = []
   slipSR = []
   for Rec_id in Rec_IDs:
      # instantiate the rectangle
      rEc = Rprop[Rec_id]['e']
      rNc = Rprop[Rec_id]['n']
      rZc = -1.0 * Rprop[Rec_id]['dep']
      L = NP.sqrt( Rprop[Rec_id]['area'] )
      W = L
      strike = Rprop[Rec_id]['strike']
      dip = Rprop[Rec_id]['dip']

      Rec = rectanglePatch( rEc, rNc, rZc, L, W, strike, dip )
      

      # compute the subfaults for the rectangle
      subRectangles = Rcalc.getReducedRectangles(Rec, eR, nR, Amax = Amax)
      eSR.extend([subR.Xc for subR in subRectangles])
      nSR.extend([subR.Yc for subR in subRectangles])
      dSR.extend([-1.0*subR.Zc for subR in subRectangles])
      aSR.extend([subR.L * subR.W for subR in subRectangles])
      idSR.extend( [Rec_id] * len(subRectangles) )
      # get angles
      strikeSR.extend([subR.strike for subR in subRectangles])
      dipSR.extend([subR.dip for subR in subRectangles])
      # assign unit slip
      slipSR.extend(1.0 * NP.ones(len(subRectangles)))

      msg = 'Square %s has %i subsources...'%(Rec_id, len(subRectangles))
      log.addLine(msg)
   
      if plotGeometry:
         subRectangles4Plot.extend(subRectangles)
         Rectangles4Plot.append(Rec)

   # Convert units
   Units2meters = method_par['EDKSunits']
   eSR = NP.array(eSR) * Units2meters
   nSR = NP.array(nSR) * Units2meters
   dSR = NP.array(dSR) * Units2meters
   aSR = NP.array(aSR) * Units2meters * Units2meters
   eR = NP.array(eR) * Units2meters
   nR = NP.array(nR) * Units2meters

   if plotGeometry:
      import pylab as PL
      import mpl_toolkits.mplot3d.axes3d as plot3
      fig = PL.figure(1)
      s = plot3.Axes3D(fig)
      for subRectangle in subRectangles4Plot:
         Rcalc.plotRectangle(s, subRectangle, Rcalc, color = 'r') 
      for MasterRectangle in Rectangles4Plot:
         Rcalc.plotRectangle(s, MasterRectangle, Rcalc, color = 'k')

      s.plot(1.0 * eR/Units2meters, 1.0 * nR/Units2meters, 'og')
      PL.show() 
   
   # the rest just convert to arrays.
   strikeSR = NP.array(strikeSR)
   dipSR = NP.array(dipSR)
   slipSR = NP.array(slipSR)
   
   # run layered_disloc_sub on all the subrectangles
   edks = method_par['EDKSfilename']
   prefix = method_par['prefix']

   # strike slip GFs
   log.addLine('calculating Strike-Slip Green functions...')
   prefixSS = prefix + '_SS_'
   rakeSR = 0.0 * NP.ones(len(eSR))
   GeSS, GnSS, GuSS = layered_disloc_sub(idSR, eSR, nSR, dSR, strikeSR, dipSR,\
                      rakeSR, slipSR, aSR, eR, nR, edks, prefixSS)  

   # dip slip GFs
   log.addLine('calculating upDip-Slip Green functions...')
   prefixDS = prefix + '_DS_'
   rakeSR = 90.0 * NP.ones(len(eSR))   
   GeDS, GnDS, GuDS = layered_disloc_sub(idSR, eSR, nSR, dSR, strikeSR, dipSR,\
                      rakeSR, slipSR, aSR, eR, nR, edks, prefixDS)  


   # if useRecvDir is True I need to project the GF matrices
   if useRecvDir:
      G_SS, G_DS = projectGFmatrices(GeSS, GnSS, GuSS, GeDS, GnDS, GuDS,\
 				     ODirE, ODirN, ODirU)  
 

   if w_bin: # write GF matrixes into binary files
      GeDS = GeDS.astype(NP.float32)
      GeDS.tofile('%s_GeDS.dat'%prefix)

      GeSS = GeSS.astype(NP.float32)
      GeSS.tofile('%s_GeSS.dat'%prefix)

      GnDS = GnDS.astype(NP.float32)
      GnDS.tofile('%s_GnDS.dat'%prefix)

      GnSS = GnSS.astype(NP.float32)
      GnSS.tofile('%s_GnSS.dat'%prefix)
      
      GuDS = GuDS.astype(NP.float32)
      GuDS.tofile('%s_GuDS.dat'%prefix)
      
      GuSS = GuSS.astype(NP.float32)
      GuSS.tofile('%s_GuSS.dat'%prefix)   

   if w_ascii:  # write GF matrices in ASCII format
      # GeDS
      file = open('GeDS.txt','w')
      Nrows, Ncols = GeDS.shape
      print 'Nrows = ' + str(Nrows) + ', Ncols = ' + str(Ncols)
      for row in range(0,Nrows):
         for col in range(0,Ncols):
            value = '%s ' %(GeDS[row][col])
            file.write(value)
         file.write('\n')
      file.close()

      #GeSS
      file = open('GeSS.txt','w')
      Nrows, Ncols = GeSS.shape
      print 'Nrows = ' + str(Nrows) + ', Ncols = ' + str(Ncols)
      for row in range(0,Nrows):
         for col in range(0,Ncols):
            value = '%s ' %(GeSS[row][col])
            file.write(value)
         file.write('\n')
      file.close()

      #####
      # GnDS
      file = open('GnDS.txt','w')
      Nrows, Ncols = GnDS.shape
      print 'Nrows = ' + str(Nrows) + ', Ncols = ' + str(Ncols)
      for row in range(0,Nrows):
         for col in range(0,Ncols):
            value = '%s ' %(GnDS[row][col])
            file.write(value)
         file.write('\n')
      file.close()

      #GnSS
      file = open('GnSS.txt','w')
      Nrows, Ncols = GnSS.shape
      print 'Nrows = ' + str(Nrows) + ', Ncols = ' + str(Ncols)
      for row in range(0,Nrows):
         for col in range(0,Ncols):
            value = '%s ' %(GnSS[row][col])
            file.write(value)
         file.write('\n')
      file.close()

      #####
      # GuDS
      file = open('GuDS.txt','w')
      Nrows, Ncols = GuDS.shape
      print 'Nrows = ' + str(Nrows) + ', Ncols = ' + str(Ncols)
      for row in range(0,Nrows):
         for col in range(0,Ncols):
            value = '%s ' %(GuDS[row][col])
            file.write(value)
         file.write('\n')
      file.close()

      #GuSS
      file = open('GuSS.txt','w')
      Nrows, Ncols = GuSS.shape
      print 'Nrows = ' + str(Nrows) + ', Ncols = ' + str(Ncols)
      for row in range(0,Nrows):
         for col in range(0,Ncols):
            value = '%s ' %(GuSS[row][col])
            file.write(value)
         file.write('\n')
      file.close()


   # if direction is used.
   if useRecvDir:
      
      if w_bin:    # write GF matrixes into binary files
         G_SS = G_SS.astype(NP.float32)
         G_SS.tofile('%s_G_SS.dat'%prefix)
         
         G_DS = G_DS.astype(NP.float32)
         G_DS.tofile('%s_G_DS.dat'%prefix)

      if w_ascii:  # write GF matrices in ASCII format
         #G_SS
         file = open('G_SS_proj.txt','w')
         Nrows, Ncols = G_SS.shape
         print 'Nrows = ' + str(Nrows) + ', Ncols = ' + str(Ncols)
         for row in range(0,Nrows):
            for col in range(0,Ncols):
               value = '%s ' %(G_SS[row][col])
               file.write(value)
            file.write('\n')
         file.close()

         #G_DS
         file = open('G_DS_proj.txt','w')
         Nrows, Ncols = G_DS.shape
         print 'Nrows = ' + str(Nrows) + ', Ncols = ' + str(Ncols)
         for row in range(0,Nrows):
            for col in range(0,Ncols):
               value = '%s ' %(G_DS[row][col])
               file.write(value)
            file.write('\n')
         file.close()
            
      return GeSS, GeDS, GnSS, GnDS, GuSS, GuDS, G_SS, G_DS      

   else:

      return GeSS, GeDS, GnSS, GnDS, GuSS, GuDS


if __name__ == '__main__':

   import os,sys
   sys.path.append(os.getcwd())
   from EDKSsubParams import *
   parNames = ['useRecvDir', 'Amax', 'EDKSunits', 'EDKSfilename', 'prefix']
   parValues = [ useRecvDir ,  Amax ,  EDKSunits ,  EDKSfilename ,  prefix ]
   method_par = dict(zip(parNames, parValues))

   calcGreenFunctions_EDKS_subSquares(SquaresPropFile, ReceiverFile, method_par, plotGeometry)
