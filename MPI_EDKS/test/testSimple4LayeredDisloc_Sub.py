import numpy as NP
from layered_disloc_sub import layered_disloc_sub
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
   Nx = 300
   x = NP.linspace(xmin,xmax,Nx)
   ymin = -150 * km
   ymax = 150 * km
   Ny = 300
   y = NP.linspace(ymin,ymax,Ny)

   XX, YY = NP.meshgrid(x,y)

   # reshape to form vectors.
   x = NP.reshape(XX, (Ny*Nx,))
   y = NP.reshape(YY, (Ny*Nx,))

   # source coordinates
   NPatches = 256
   NumSourcesPerPatch = 4096*2
   Nsources = NumSourcesPerPatch * NPatches
   xS = 0.0 * NP.ones(Nsources) * km
   yS = 0.0 * NP.ones(Nsources) * km
   zS = 30.0 * NP.ones(Nsources) * km
   aS = 300*100.0 * NP.ones(Nsources) * km * km
   strike =  0.0 * NP.ones(Nsources)
   dip = 90.0 * NP.ones(Nsources)  
   SSlip = 0.0 * NP.ones(Nsources)  # in meters
   DSlip = 1.0 * NP.ones(Nsources) * km # in meters
   IDs = []
   for i in range(0,NPatches):
      IDs.extend([str(i)]*NumSourcesPerPatch )
 
   # calculate the GF's
   rake = 0 * NP.ones(SSlip.shape)
   slip = 1 * NP.ones(SSlip.shape)
   GFeSS, GFnSS, GFzSS = layered_disloc_sub(IDs, xS, yS, zS, strike, dip, rake, slip, aS,\
                                        x, y, edks, prefix, BIN_EDKS)  

   rake = 90 * NP.ones(SSlip.shape)
   slip = 1 * NP.ones(SSlip.shape)
   GFeDS, GFnDS, GFzDS = layered_disloc_sub(IDs, xS, yS, zS, strike, dip, rake, slip, aS,\
                                        x, y, edks, prefix, BIN_EDKS)

   # compute forward displacement calculation.

   dE = GFeSS * SSlip[0] + GFeDS * DSlip[0]
   dN = GFnSS * SSlip[0] + GFnDS * DSlip[0]
   dZ = GFzSS * SSlip[0] + GFzDS * DSlip[0]

   print dE[0]
   print dN[0]
   print dZ[0]
   print "  ----- "
   print dE
   return None
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
