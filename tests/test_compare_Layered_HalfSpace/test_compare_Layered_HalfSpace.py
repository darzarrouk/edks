#!/usr/bin/env python
"""
uses the tohoku_nied_2.edks kernels and another one generated with a layered half space
(poisson_half_space.edks) to compare the predictions of both models for a thrust focal mechanism 

"""
import numpy as NP
import pylab as PL
import matplotlib.pyplot as plt

from layered_disloc_sub import layered_disloc_sub

km = 1000.0
# focal mechanism.
strike = 0.0
dip = 17.0
rake = 90.0
slip = 1.0
area = 30.0

# receivers
Rmin = 4.5*4
Rmax = 150.0  
dR = 1.5/10.0
eR = NP.arange(Rmin, Rmax, dR) * km
nR = NP.zeros(eR.shape) * km
Nrec = len(eR)

# depths of sources
Dmin = 5.0
Dmax = 75.0
dD = 0.1
depS = NP.arange(Dmin, Dmax, dD) * km
Ndep = len(depS)

# setup the sources
strike = strike * NP.ones(Ndep)
dip = dip * NP.ones(Ndep)
rake = rake * NP.ones(Ndep)
slip = slip * NP.ones(Ndep)
area = area * NP.ones(Ndep) * km * km
eS = NP.zeros(Ndep) * km
nS = NP.zeros(Ndep) * km
idS = [str(i) for i in range(0, Ndep)]

# calculate the displacements for tohoku_nied_2 layered space
edks = 'tohoku_nied_2.edks'
prefix = 'Disp_tohoku_nied_2'
BIN_EDKS = '${HOME}/dev/edks/bin'

Te, Tn, Tu = layered_disloc_sub(idS, eS, nS, depS, strike, dip, rake, slip, area,\
                                                  eR, nR, edks, prefix, BIN_EDKS)

Te = Te.transpose()
Tn = Tn.transpose()
Tu = Tu.transpose()

# calculate the homogeneus layered space (poisson solid)
edks = 'poisson_half_space.edks'
prefix = 'Disp_poisson_half_space'
BIN_EDKS = '${HOME}/dev/edks/bin'

He, Hn, Hu = layered_disloc_sub(idS, eS, nS, depS, strike, dip, rake, slip, area,\
                                                  eR, nR, edks, prefix, BIN_EDKS)

He = He.transpose()
Hn = Hn.transpose()
Hu = Hu.transpose()


print Ndep, Nrec
print Te.shape, Tn.shape, Tu.shape, He.shape, Hn.shape, Hu.shape

# save arrays
NP.savetxt('Te.dat', Te)
NP.savetxt('Tn.dat', Tn)
NP.savetxt('Tu.dat', Tu)

NP.savetxt('He.dat', He)
NP.savetxt('Hn.dat', Hn)
NP.savetxt('Hu.dat', Hu)

NP.savetxt('Dist.dat', eR/km)
NP.savetxt('Dep.dat', depS/km)

