#!/usr/bin/env python
"""
test that compares the GFs of EDKS and OKADA for an homogeneous half space in a predefined grid and source.
"""

import numpy as NP
import string
from ICM.Physics.Okada85.okada85 import calcOkada85
import os

###
def main():
   # units to meters
   km = 1000.0

   # create the equispaced grid in KMs for the receivers
   Xmin = -200
   Xmax = 200
   Nx = 200
   x = NP.linspace(Xmin, Xmax, Nx)
   Ymin = -200
   Ymax = 200
   Ny = 200
   y = NP.linspace(Ymin, Ymax, Ny)
   # create the mesh and reshape it to vectors.
   XX, YY = NP.meshgrid(x,y)
   xR = NP.reshape(XX, (Ny*Nx,))
   yR = NP.reshape(YY, (Ny*Nx,))

   
   # define the source coordinates and properties
   xSmin = -50.
   xSmax = 50.
   nSx = 5
   xxS, yyS = NP.meshgrid(NP.linspace(xSmin, xSmax, nSx), NP.linspace(xSmin, xSmax, nSx)) 
   xS = NP.reshape(xxS, (nSx*nSx,))
   yS = NP.reshape(yyS, (nSx*nSx,))
   zSmin = 10.
   zSmax = 40.
   zS = NP.linspace(zSmin, zSmax, yS.size) # depth
   aS = 0.9 * zS * zS # source area
   strike = NP.linspace(10.,360., yS.size) 
   dip =  NP.linspace(12.0, 90.0, yS.size)
   rake = NP.linspace(-150., 180., yS.size)
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
      L = NP.sqrt(aS[i])
      W = L
      Ux, Uy, Uz = calcOkada85(slip[i], strike[i], dip[i], \
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
   Rmax = 1.5 * NP.sqrt(2) * (xSmax - xSmin) # 1.5 is a security factor
   dR = 0.5
   depths = NP.exp(NP.linspace( NP.log(zSmin-1), NP.log(zSmax+1), 2 * len(zS) ) )
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
   cmd = "mpirun -n 2 MPI_EDKS.py test_1a.edks_config"
   status = os.system(cmd)
   if status != 0:
      raise ValueError('Unable to compute the EDKS kernels !!!')
   else:
      os.system('rm -rf *.asc')
   # now compute the displacements for the given sources at the receivers.
   


if __name__ == '__main__':
   main()
