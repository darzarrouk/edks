#!/usr/bin/env python
"""
test that compares the GFs of EDKS and OKADA for an homogeneous half space in a predefined grid and source.
"""

import numpy as NP
import string
from ICM.Physics.Okada85.okada85 import calcOkada85
import os
from layered_disloc_sub import layered_disloc_sub
import pylab as PL

###
def main():
   print(" THIS MAY TAKE A WHILE TO RUN !!! ")
   # units to meters
   km = 1000.0

   # create the equispaced grid in KMs for the receivers
   Xmin = -100
   Xmax = 100
   Nx = 50
   x = NP.linspace(Xmin, Xmax, Nx)
   Ymin = -100
   Ymax = 100
   Ny = 50
   y = NP.linspace(Ymin, Ymax, Ny)
   # create the mesh and reshape it to vectors.
   XX, YY = NP.meshgrid(x,y)
   xR = NP.reshape(XX, (Ny*Nx,))
   yR = NP.reshape(YY, (Ny*Nx,))

   
   # define the source coordinates and properties
   xSmin = -10.
   xSmax = 10.
   nSx = 15
   xxS, yyS = NP.meshgrid(NP.linspace(xSmin, xSmax, nSx), NP.linspace(xSmin, xSmax, nSx)) 
   xS = NP.reshape(xxS, (nSx*nSx,))
   yS = NP.reshape(yyS, (nSx*nSx,))
   zSmin = 10.
   zSmax = 40.
   zS = NP.linspace(zSmin, zSmax, yS.size) # depth
   aS = 0.9 * zS * zS # source area
   strike = NP.linspace(2.,360., yS.size) 
   dip =  NP.linspace(5.0, 90, yS.size)
   rake = NP.linspace(-179., 179.99, yS.size)
   slip = 1.0 * NP.ones(yS.shape)
   # define the properties of the half space (I am taking the first layer from the Tohoku
   # Oki velocity model).
   vp = 4.4 # km/s
   vs = 2.51 # km/s
   rho = 3.3 # kg/m3
   Mu = rho * vs * vs * km * km # kg/m/s
   Lambda = rho * vp * vp * km * km - 2 * Mu # kg/m/s
   Nu = Lambda / ( 2. * ( Lambda + Mu ) ) # Poisson Ratio
   
   # Calculate the Okada displacements, (Okada code requires coordinates in KMs)
   # create container for different sources:
   OKADA_disp = []
   for i in range(0, len(xS)):
      # in order to model a point source L=W = 1km and slip = slip[i] * aS(i)
      divLW = 100.0
      L = 1./divLW
      W = 1./divLW
      Ux, Uy, Uz = calcOkada85(slip[i]*aS[i]*divLW*divLW, strike[i], dip[i], \
           rake[i], L, W, xS[i], yS[i], zS[i], xR, yR, Nu)
      OKADA_disp.append([Ux, Uy, Uz])
   
   # Now calculate the EDKS displacements (EDKS requires all the units in meters)
   # first calculate the EDKS kernels
   # velocity model file
   edks_model="""1 %.1f
%.2f   %.2f   %.2f   0.0
""" %(km, rho, vp, vs) # last 0.0 is layer thickness (indicating half space).
   file = open('test_1a.model','w')
   file.write(edks_model)
   file.close()
   # edks_config file
   Rmax = 1.1 * NP.sqrt(2) * (Xmax - Xmin) # 1.5 is a security factor
   dR = 0.5
   depths = NP.exp(NP.linspace( NP.log(zSmin-2), NP.log(zSmax+2), 200 ) )
   deps = ""
   for dep in depths:
      deps += '%.3f '%(dep)
   edks_config = """BIN_tab4E        ${EDKS_HOME}/bin/tab4E
BIN_build_edks   ${EDKS_HOME}/bin/build_edks
ModPrefix        test_1a
dR               %.3f
Rmax             %.3f
depths           %s
""" %(dR, Rmax, deps)
   file = open('test_1a.edks_config','w')
   file.write(edks_config)
   file.close()
   # call MPI_EDKS to calculate the GFs.
   #cmd = "mpirun -n 5 MPI_EDKS.py test_1a.edks_config"
   #status = os.system(cmd)
   #if status != 0:
   #   raise ValueError('Unable to compute the EDKS kernels !!!')
   #else:
   #   os.system('rm -rf *.asc')
   # now compute the displacements for the given sources at the receivers.
   # must have units in meters
   print("calculating the displacements using layered_disloc")
   idS = [str(bla) for bla in range(0, len(xS))]
   ux, uy, uz = layered_disloc_sub(idS, xS*km, yS*km, zS*km, strike, dip, rake, slip,\
                 aS*km*km, xR*km, yR*km, 'test_1a.edks', 'test_1a') 
   
   # reassemble the matrices for comparison
   EDKS_disp = [[ux.transpose()[i][:], uy.transpose()[i][:], uz.transpose()[i][:]]\
                for i in range(0, len(xS)) ]

   # now compute the norm of the difference for all the sources
   MeanUdiff_X = []
   MeanUdiff_Y = []
   MeanUdiff_Z = [] 
   SdevUdiff_X = []
   SdevUdiff_Y = []
   SdevUdiff_Z = []
   SlopeDevUedksUokada_X = [] # delta = NP.abs(45 - NP.arctan(m2) * 180 / NP.pi)
   SlopeDevUedksUokada_Y = [] # delta = NP.abs(45 - NP.arctan(m2) * 180 / NP.pi)
   SlopeDevUedksUokada_Z = [] # delta = NP.abs(45 - NP.arctan(m2) * 180 / NP.pi)
   MaxDiff_X = []
   MaxDiff_Y = []
   MaxDiff_Z = []

   for i in range(0, len(xS)):
      UXokada, UYokada, UZokada = OKADA_disp[i]
      UXedks, UYedks, UZedks = EDKS_disp[i]
      MeanUdiff_X.append( NP.mean(UXokada - UXedks) )
      MeanUdiff_Y.append( NP.mean(UYokada - UYedks) )
      MeanUdiff_Z.append( NP.mean(UZokada - UZedks) )
      SdevUdiff_X.append( NP.std(UXokada - UXedks) )
      SdevUdiff_Y.append( NP.std(UYokada - UYedks) )
      SdevUdiff_Z.append( NP.std(UZokada - UZedks) )
      mX, nX = rectLineLSQfit(UXedks, UXokada)
      mY, nY = rectLineLSQfit(UYedks, UYokada)
      mZ, nZ = rectLineLSQfit(UZedks, UZokada)
      SlopeDevUedksUokada_X.append( NP.abs(45 - NP.arctan(mX) * 180 / NP.pi) )
      SlopeDevUedksUokada_Y.append( NP.abs(45 - NP.arctan(mY) * 180 / NP.pi) )
      SlopeDevUedksUokada_Z.append( NP.abs(45 - NP.arctan(mZ) * 180 / NP.pi) )
      MaxDiff_X.append( NP.max( NP.abs( UXokada - UXedks ) ) )
      MaxDiff_Y.append( NP.max( NP.abs( UYokada - UYedks ) ) )
      MaxDiff_Z.append( NP.max( NP.abs( UZokada - UZedks ) ) )

   # plot these quantities
   # first concatenate all the displacements and plot
 
   PL.figure()
   s = PL.subplot(4,3,1)
   s.plot(MeanUdiff_X, 'ob')
   s.set_title('mean(diff(Ux))')
   s = PL.subplot(4,3,2)
   s.plot(MeanUdiff_Y, 'ob')
   s.set_title('mean(diff(Uy))')
   s = PL.subplot(4,3,3)
   s.plot(MeanUdiff_Z, 'ob')
   s.set_title('mean(diff(Uz))')
   s = PL.subplot(4,3,4)
   s.plot(SdevUdiff_X, 'ob')
   s.set_title('Sdev(diff(Ux))')
   s = PL.subplot(4,3,5)
   s.plot(SdevUdiff_Y, 'ob')
   s.set_title('Sdev(diff(Uy))')
   s = PL.subplot(4,3,6)
   s.plot(SdevUdiff_Z, 'ob')
   s.set_title('Sdev(diff(Uz))')
   s = PL.subplot(4,3,7)
   s.plot(SlopeDevUedksUokada_X, 'ob')
   s.set_title('Slope_EDKS/OKADA(Ux)')
   s = PL.subplot(4,3,8)
   s.plot(SlopeDevUedksUokada_Y, 'ob')
   s.set_title('Slope_EDKS/OKADA(Uy)')
   s = PL.subplot(4,3,9)
   s.plot(SlopeDevUedksUokada_Z, 'ob')
   s.set_title('Slope_EDKS/OKADA(Uz)')
   s = PL.subplot(4,3,10)
   s.plot(MaxDiff_X, 'ob')
   s.set_title('Max(abs(diff(Ux)))')
   s = PL.subplot(4,3,11)
   s.plot(MaxDiff_Y, 'ob')
   s.set_title('Max(abs(diff(Uy)))')
   s = PL.subplot(4,3,12)
   s.plot(MaxDiff_Z, 'ob')
   s.set_title('Max(abs(diff(Uz)))')


   PL.show()

   
def plotBox1(s, X, Y1, markersize, labelX, labelY1, title,\
             useZoom = False, Fsigma = 1):
   import numpy as NP
   test1 = NP.isnan(X)
   test2 = NP.isnan(Y1)
   test = test1 + test2
   test = 1 - test
   I = NP.where(test)[0]
   X = X[I]
   Y1 = Y1[I]
   Y2 = Y2[I]

   valMin = NP.min([NP.min(X), NP.min(Y1)])
   valMax = NP.max([NP.max(X), NP.max(Y1)])
   # do a LS fit to each X,Y pair and put m, n in label
   m1, n1 = rectLineLSQfit(X, Y1)
   delta1 = NP.abs(45 - NP.arctan(m1) * 180 / NP.pi)
   L1x = NP.array([valMin, valMax])
   L1y = L1x * m1 + n1
   legendY1 = labelY1+': delta= %.2f (%.1f%%), n= %.1E'%(delta1,(m1-1)*100,n1)
   s.plot(X, Y1, 'xr', markersize = markersize, label = legendY1)
   s.plot([valMin, valMax], [valMin, valMax],'-k')
   s.legend(loc = 'upper left')
   s.set_xlabel(labelX)
   s.set_ylabel(labelY1)
   s.set_title(title)
   s.plot(L1x, L1y, ':r')
   s.axis('equal')

   if useZoom:
      meanX = NP.mean(X)
      sigmaX = NP.std(X)
      aux = list(Y1)
      meanY = NP.mean(aux)
      sigmaY = NP.std(aux)
      sigma = NP.max([sigmaX, sigmaY])
      Xmin, Xmax, Ymin, Ymax = s.axis()
      XminZoom = meanX - Fsigma * sigma
      XmaxZoom = meanX + Fsigma * sigma
      YminZoom = meanY - Fsigma * sigma
      YmaxZoom = meanY + Fsigma * sigma
      s.set_xlim((XminZoom, XmaxZoom))
      s.set_ylim((YminZoom, YmaxZoom))


def plotBox3(s, X, Y1, Y2, markersize, labelX, labelY1, labelY2, title,\
             useZoom = False, Fsigma = 1):
   import numpy as NP
   test1 = NP.isnan(X)
   test2 = NP.isnan(Y1)
   test3 = NP.isnan(Y2)
   test = test1 + test2 + test3
   test = 1 - test
   I = NP.where(test)[0]
   X = X[I]
   Y1 = Y1[I]
   Y2 = Y2[I]

   valMin = NP.min([NP.min(X), NP.min(Y1), NP.min(Y2)])
   valMax = NP.max([NP.max(X), NP.max(Y1), NP.max(Y2)])
   # do a LS fit to each X,Y pair and put m, n in label
   m1, n1 = rectLineLSQfit(X, Y1)
   delta1 = NP.abs(45 - NP.arctan(m1) * 180 / NP.pi)
   L1x = NP.array([valMin, valMax])
   L1y = L1x * m1 + n1
   m2, n2 = rectLineLSQfit(X, Y2)
   delta2 = NP.abs(45 - NP.arctan(m2) * 180 / NP.pi)
   L2x = NP.array([valMin, valMax])
   L2y = L2x * m2 + n2
   legendY1 = labelY1+': delta= %.2f (%.1f%%), n= %.1E'%(delta1,(m1-1)*100,n1)
   legendY2 = labelY2+': delta= %.2f (%.1f%%), n= %.1E'%(delta2,(m2-1)*100,n2)
   s.plot(X, Y1, 'xr', markersize = markersize, label = legendY1)
   s.plot(X, Y2, '+b', markersize = markersize, label = legendY2)
   s.plot([valMin, valMax], [valMin, valMax],'-k')
   s.legend(loc = 'upper left')
   s.set_xlabel(labelX)
   s.set_ylabel(labelY1 + ',' + labelY2)
   s.set_title(title)
   s.plot(L1x, L1y, ':r')
   s.plot(L2x, L2y, ':b')
   s.axis('equal')

   if useZoom:
      meanX = NP.mean(X)
      sigmaX = NP.std(X)
      aux = list(Y1)
      aux.extend(Y2)
      meanY = NP.mean(aux)
      sigmaY = NP.std(aux)
      sigma = NP.max([sigmaX, sigmaY])
      Xmin, Xmax, Ymin, Ymax = s.axis()
      XminZoom = meanX - Fsigma * sigma
      XmaxZoom = meanX + Fsigma * sigma
      YminZoom = meanY - Fsigma * sigma
      YmaxZoom = meanY + Fsigma * sigma
      s.set_xlim((XminZoom, XmaxZoom))
      s.set_ylim((YminZoom, YmaxZoom))




###
def rectLineLSQfit(X, Y):
   import numpy as NP
   test1 = NP.isnan(X)
   test2 = NP.isnan(Y)
   test = test1 + test2
   test = 1 - test
   I = NP.where(test)[0]
   X = X[I]
   Y = Y[I]
   L = len(X)
   G = NP.ones((L,2))
   G.transpose()[0][:] = X
   aux = NP.dot(G.transpose(), G)
   aux = NP.linalg.inv(aux)
   aux = NP.dot(aux, G.transpose())
   return NP.dot(aux, Y)
   

   

if __name__ == '__main__':
   main()
