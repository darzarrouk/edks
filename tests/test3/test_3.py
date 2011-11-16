#!/usr/bin/env python
"""
test that compares the GFs of EDKS and OKADA for an homogeneous half space in a predefined grid and source.
"""

import numpy as NP
import string
from ICM.Physics.Okada85.okada85 import calcOkada85
import os
from layered_disloc_sub import layered_disloc_sub
from layered_disloc import layered_disloc
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
   xS = NP.array([0.0])
   yS = NP.array([0.0])
   zS = NP.array([48.0])
   Ls = NP.array([2.])
   Ws = NP.array([2. ])
   aS = NP.array([Ls * Ws]) # source area
   strike = NP.array([0.0]) 
   dip =  NP.array([90.0])
   rake = NP.array([90.0])
   slip = NP.array([1.0])
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
   OUx, OUy, OUz = calcOkada85(slip[0], strike[0], dip[0], \
                rake[0], Ls, Ws, xS[0], yS[0], zS[0], xR, yR, Nu)
   
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
   zSmin = 5.0
   zSmax = 50.0
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
   cmd = "mpirun -n 5 MPI_EDKS.py test_1a.edks_config"
   #status = os.system(cmd)
   #if status != 0:
   #   raise ValueError('Unable to compute the EDKS kernels !!!')
   #else:
   #   os.system('rm -rf *.asc')



   # now compute the displacements for the given sources at the receivers.
   # must have units in meters
   print("calculating the displacements using layered_disloc_sub")
   # here I am going to divide the rectangle into several subfaults
   # this is only valid for strike = 0, dip = 90
   NyS = 30
   NzS = 30
   dy = 1.0 * Ls / NyS
   dz = 1.0 * Ws / NzS
   yS = NP.linspace(-Ls/2.0 + dy/2.0 , Ls/2.0 - dy/2.0, NyS-1) 
   zS = NP.linspace( zS[0] - Ws/2.0 + dz/2.0 , zS[0] + Ws/2.0 - dz/2.0, NzS-1)
   yyS, zzS = NP.meshgrid(yS, zS)
   yS = NP.reshape(yyS, yyS.size)
   zS = NP.reshape(zzS, zzS.size)
   xS = NP.zeros(zS.shape)
   aS = (1.0 * aS[0]/ len(zS)) * NP.ones(zS.shape)
   strike = strike[0] * NP.ones(zS.shape)
   dip = dip[0] * NP.ones(zS.shape)
   rake = rake[0] * NP.ones(zS.shape)
   slip = slip[0] * NP.ones(zS.shape)

   IDs = ["MAIN" for i in range(0, len(xS))]
   EUx, EUy, EUz = layered_disloc_sub(IDs, xS*km, yS*km, zS*km, strike, dip, rake, slip,\
                 aS*km*km, xR*km, yR*km, 'test_1a.edks', 'test_1a') 
   EUx2, EUy2, EUz2 = layered_disloc(xS*km, yS*km, zS*km, strike, dip, rake, slip,\
                 aS*km*km, xR*km, yR*km, 'test_1a.edks', 'test_1a')
   NEUx = NP.sum(EUx2, 1)
   NEUy = NP.sum(EUy2, 1)
   NEUz = NP.sum(EUz2, 1)


   PL.figure()
   s = PL.subplot(1,3,1)
   s.plot( NEUx, EUx, '.r')
   valMin = NP.min([ NP.min(NEUx), NP.min(EUx)])
   valMax = NP.max([ NP.max(NEUx), NP.max(EUx)])
   s.plot([valMin, valMax], [valMin, valMax],'-k')
   PL.axis('equal')

   s = PL.subplot(1,3,2)
   s.plot( NEUy, EUy, '.r')
   valMin = NP.min([NP.min(NEUy), NP.min(EUy)])
   valMax = NP.max([NP.max(NEUy), NP.max(EUy)])
   s.plot([valMin, valMax], [valMin, valMax],'-k')
   PL.axis('equal')

   s = PL.subplot(1,3,3)
   s.plot( NEUz, EUz, '.r')
   valMin = NP.min([NP.min(NEUz), NP.min(EUz)])
   valMax = NP.max([NP.max(NEUz), NP.max(EUz)])
   s.plot([valMin, valMax], [valMin, valMax],'-k')
   PL.axis('equal')


   PL.figure()
   s = PL.subplot(1,3,1)
   s.plot( OUx, EUx, '.r')
   valMin = NP.min([ NP.min(OUx), NP.min(EUx)])
   valMax = NP.max([ NP.max(OUx), NP.max(EUx)])
   s.plot([valMin, valMax], [valMin, valMax],'-k')
   PL.axis('equal')

   s = PL.subplot(1,3,2)
   s.plot( OUy, EUy, '.r')
   valMin = NP.min([NP.min(OUy), NP.min(EUy)])
   valMax = NP.max([NP.max(OUy), NP.max(EUy)])
   s.plot([valMin, valMax], [valMin, valMax],'-k')
   PL.axis('equal')

   s = PL.subplot(1,3,3)
   s.plot( OUz, EUz, '.r')
   valMin = NP.min([NP.min(OUz), NP.min(EUz)])
   valMax = NP.max([NP.max(OUz), NP.max(EUz)])
   s.plot([valMin, valMax], [valMin, valMax],'-k')
   PL.axis('equal')

   PL.show()  


 

if __name__ == '__main__':
   main()
