import numpy as NP
from layered_disloc import layered_disloc
import pylab as PL

def test1():
   # units
   km = 1000.0 # in meters
   # definitions
   edks = 'tohoku_3.edks'
   prefix = 'Test_tohoku_3_'
   BIN_EDKS = '${HOME}/dev/edks/bin'
   # define grid
   xmin = -150 * km
   xmax = 150 * km
   Nx = 3
   x = NP.linspace(xmin,xmax,Nx)
   ymin = -150 * km
   ymax = 150 * km
   Ny = 3
   y = NP.linspace(ymin,ymax,Ny)

   XX, YY = NP.meshgrid(x,y)

   # reshape to form vectors.
   x = NP.reshape(XX, (Ny*Nx,))
   y = NP.reshape(YY, (Ny*Nx,))

   # source coordinates
   xS = NP.array([0.0, 0.0]) * km
   yS = NP.array([0.0, 0.0]) * km
   zS = NP.array([30, 30]) * km
   aS = NP.array([300 * 100, 300 * 100]) * km * km
   strike = NP.array([0.0, 0.0])
   dip = NP.array([90.0, 90.0])
   SSlip = NP.array([0.0, 0.0]) # in meters
   DSlip = NP.array([1.0, 1.0]) # in meters

   # calculate the GF's
   rake = 0 * NP.ones(SSlip.shape)
   slip = 1 * NP.ones(SSlip.shape)
   GFeSS, GFnSS, GFzSS = layered_disloc(xS, yS, zS, strike, dip, rake, slip, aS, x, y,\
                                     edks, prefix, BIN_EDKS)  

   rake = 90 * NP.ones(SSlip.shape)
   slip = 1 * NP.ones(SSlip.shape)
   GFeDS, GFnDS, GFzDS = layered_disloc(xS, yS, zS, strike, dip, rake, slip, aS, x, y,\
                                     edks, prefix, BIN_EDKS)

   # compute forward displacement calculation.

   dE = GFeSS * SSlip[0] + GFeDS * DSlip[0]
   dN = GFnSS * SSlip[0] + GFnDS * DSlip[0]
   dZ = GFzSS * SSlip[0] + GFzDS * DSlip[0]

   print dE

   dEE = NP.reshape(dE, (Ny, Nx))
   dNN = NP.reshape(dN, (Ny, Nx))
   dZZ = NP.reshape(dZ, (Ny, Nx))

   # plot the vector fields.
   import pylab as PL
   #PL.figure()
   norm = NP.sqrt(dE * dE + dN * dN)
   #PL.quiver(x, y, dE/norm, dN/norm)
   #PL.axis('equal')
   #PL.figure()
   #PL.quiver(x, y, dZ * 0.0, dZ/NP.abs(dZ))
   #PL.axis('equal')

   PL.figure()
   PL.subplot(1,3,1)
   PL.pcolor(XX, YY, dEE)
   PL.colorbar()
   PL.axis('equal')
   PL.subplot(1,3,2)
   PL.pcolor(XX, YY, dNN)
   PL.colorbar()
   PL.axis('equal')
   PL.subplot(1,3,3)
   PL.pcolor(XX, YY, dZZ)
   PL.axis('equal')
   PL.colorbar()


   PL.show()

if __name__ == '__main__':
   test1()
