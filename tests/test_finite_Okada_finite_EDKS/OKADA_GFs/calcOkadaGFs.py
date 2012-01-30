#!/usr/bin/env python
"""
Creates Okada GFs for a single finite square patch, and also saves the files needed to compute the GFs using the 3D square wrapper for EDKS
"""

import numpy as NP
import string
from ICM.Physics.Okada85.okada85 import calcOkada85
import pylab as PL

###
def main():

   # create the equispaced grid in KMs for the receivers
   Xmin = -300
   Xmax = 300
   Nx = 100
   x = NP.linspace(Xmin, Xmax, Nx)
   Ymin = -300
   Ymax = 300
   Ny = 100
   y = NP.linspace(Ymin, Ymax, Ny)
   # create the mesh and reshape it to vectors.
   XX, YY = NP.meshgrid(x,y)
   xR = NP.reshape(XX, (Ny*Nx,))
   yR = NP.reshape(YY, (Ny*Nx,))

   
   # define the source coordinates and properties
   xSmin = -20.
   xSmax = 20.
   nSx = 10
   xxS, yyS = NP.meshgrid(NP.linspace(xSmin, xSmax, nSx), NP.linspace(xSmin, xSmax, nSx)) 
   xS = NP.reshape(xxS, (nSx*nSx,))
   yS = NP.reshape(yyS, (nSx*nSx,))
   zSmin = 10.
   zSmax = 40.
   zS = NP.linspace(zSmin, zSmax, yS.size) # depth
   L = zS * 0.75 # square patch
   W = zS * 0.75
   strike = NP.linspace(2.,360., yS.size) 
   dip =  NP.linspace(5.0, 45, yS.size)
   rake = NP.linspace(-90., 90, yS.size)
   slip = 1.0 * NP.ones(yS.shape)
   # define the properties of the half space (I am taking the first layer from the Tohoku
   # Oki velocity model).
   #vp = 4.4 # km/s
   #vs = 2.51 # km/s
   #rho = 3.3 # kg/m3
   #Mu = rho * vs * vs * km * km # kg/m/s
   #Lambda = rho * vp * vp * km * km - 2 * Mu # kg/m/s
   #Nu = Lambda / ( 2. * ( Lambda + Mu ) ) # Poisson Ratio
   Nu = 0.25 # In EDKS I am using a poissonian solid
   # Calculate the Okada displacements, (Okada code requires coordinates in KMs)
   # create container for different sources:
   GeSS = NP.zeros((len(xS), len(xR)))
   GnSS = NP.zeros((len(xS), len(xR)))
   GuSS = NP.zeros((len(xS), len(xR)))
   GeDS = NP.zeros((len(xS), len(xR)))
   GnDS = NP.zeros((len(xS), len(xR)))
   GuDS = NP.zeros((len(xS), len(xR)))

   
   for i in range(0, len(xS)):
      print 100.*i/len(xS) 
      rake = 0 * NP.ones(yS.shape)
      Ux, Uy, Uz = calcOkada85(slip[i], strike[i], dip[i], \
           rake[i], L[i], W[i], xS[i], yS[i], zS[i], xR, yR, Nu)
      GeSS[i][:] = Ux
      GnSS[i][:] = Uy
      GuSS[i][:] = Uz

      rake = 90 * NP.ones(yS.shape)
      Ux, Uy, Uz = calcOkada85(slip[i], strike[i], dip[i], \
           rake[i], L[i], W[i], xS[i], yS[i], zS[i], xR, yR, Nu)
      GeDS[i][:] = Ux
      GnDS[i][:] = Uy
      GuDS[i][:] = Uz

   
   GeSS = GeSS.transpose()
   GeDS = GeDS.transpose()
   GnSS = GnSS.transpose()
   GnDS = GnDS.transpose()
   GuSS = GuSS.transpose()
   GuDS = GuDS.transpose()

   # save GFs
   NP.savetxt('O_GeSS.txt', GeSS)
   NP.savetxt('O_GnSS.txt', GnSS)
   NP.savetxt('O_GuSS.txt', GuSS)
   NP.savetxt('O_GeDS.txt', GeDS)
   NP.savetxt('O_GnDS.txt', GnDS)
   NP.savetxt('O_GuDS.txt', GuDS)

   # write input files for EDKS 3D squares GF calculation
   # 3d SquareGridFile
   # file format 
   #E[km] N[km] Dep[km] strike dip Area ID
   file = open('3D_squares.END','w')
   header = '#E[km] N[km] Dep[km] strike dip Area ID \n'
   file.write(header)
   for i in range(0, len(xS)):
      line  = '%.3f %.3f %.3f ' % (xS[i], yS[i], zS[i])
      line += '%.3f %.3f %.3f %04i\n' % (strike[i], dip[i], L[i]*W[i], i) 
      file.write(line)
   file.close()

   # write the receiver file
   header = 'obsId E N\n'
   file = open('receivers.idEN','w')
   file.write(header)
   for i in range(0, len(xR)):
      line = '%05i %.3f %.3f \n' %(i, xR[i], yR[i])
      file.write(line)
   file.close()


if __name__ == '__main__':
   main()
