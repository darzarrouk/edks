#!/usr/bin/env python
import numpy as NP
import string
import os


# first concatenate the GF matrices
os.system('paste -d\  GeSS.txt GeDS.txt > Ge.txt')
os.system('paste -d\  GnSS.txt GnDS.txt > Gn.txt')
os.system('paste -d\  GuSS.txt GuDS.txt > Gu.txt')

os.system('cat Ge.txt Gn.txt Gu.txt > G.quickInv')

# now get the datafile
# I need to save 3 files, 
# d_GPS.dat = one long list with dE;dN;dU
# Wv_GPS.dat = one with the weights 1/sigE^2; 1/sigN^2; 1/sigU^2
# data_GPS.dat = stnm,lat,lon,ux,uy,uz

# file with the name and coordinates of the sites used
# format: obsId obs stdev ODirE ODirN ODirU E N U CODE lon lat
VectorObsFilename = 'VecObs/JapanGPS/JapanGPSObs.idEN'
# actual sites with observations and uncertainties.
# format : CODE lon lat dE dN dU SdE SdN SdU
dataInFile = 'VecObs/JapanGPS/Tohoku_PostSeismic_2011_10_25.dat'

# I only need the site CODE from the VecObs file, in the same order, so 
# I can build up the 2 other files
VecObs = []
file = open(VectorObsFilename, 'r')
# skip header
line = file.readline()
for line in file:
   line = string.split(line)
   CODE = line[0]
   VecObs.append(CODE)
file.close()

# read now all the data from the dataInFile
dataIn = {}
file = open(dataInFile, 'r')
# skip header
line = file.readline()
for line in file:
   line = string.split(line)
   CODE = line[0]
   info = line[1:]
   dataIn[CODE] = info
file.close()


# d_GPS.dat:  = one long list with dE;dN;dU
file = open('d_GPS.quickInv', 'w')
# for dE
for CODE in VecObs:
   dE = dataIn[CODE][2]
   file.write(dE + '\n')
# for dN
for CODE in VecObs:
   dN = dataIn[CODE][3]
   file.write(dN + '\n')
# for dU
for CODE in VecObs:
   dU = dataIn[CODE][4]
   file.write(dU + '\n')
file.close()

# # Wv_GPS.dat = one with the weights 1/sigE^2; 1/sigN^2; 1/sigU^2
file = open('Wv_GPS.quickInv', 'w')
# for E
for CODE in VecObs:
   SdE = float(dataIn[CODE][5])
   file.write(str(1.0/(SdE**2.0)) + '\n')
# for N
for CODE in VecObs:
   SdN = float(dataIn[CODE][6])
   file.write(str(1.0/(SdN**2.0)) + '\n')
# for U
for CODE in VecObs:
   SdU = float(dataIn[CODE][7])
   file.write(str(1.0/(SdU**2.0)) + '\n')
file.close()

# data_GPS.dat = stnm,lat,lon,ux,uy,uz
file = open('data_GPS.quickInv', 'w')
for CODE in VecObs:
   info = dataIn[CODE]
   lon = info[0] 
   lat = info[1] 
   dE = info[2] 
   dN = info[3] 
   dU = info[4]
   line = '%s %s %s %s %s %s\n' %(CODE, lat, lon, dE, dN, dU)
   file.write(line)
file.close()




