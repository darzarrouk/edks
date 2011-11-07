#!/usr/bin/env python
import numpy as NP
from layered_disloc_sub import layered_disloc_sub
import string 

def read_srcFile(filename):
   file = open(filename,'r')
   # initialize containers
   e = []
   n = []
   dep = []
   strike = []
   dip = []
   area = []
   IDs = []
   for line in file:
      line = string.split(line)
      e.append( float( line[0] ) )
      n.append( float( line[1] ) )
      dep.append( float( line[2] ) )
      strike.append( float( line[3] ) )
      dip.append( float( line[4] ) )
      area.append( float( line[5] ) )
      IDs.append( line[6] )
   file.close()
   # convert to numpy arrays
   e = NP.array(e)
   n = NP.array(n)
   dep = NP.array(dep)
   strike = NP.array(strike)
   dip = NP.array(dip)
   area = NP.array(area)
   return [e, n, dep, strike, dip , area, IDs]


def read_recvFile(filename):
   file = open(filename,'r')
   # initialize containers
   e = []
   n = []
   for line in file:
      line = string.split(line)
      e.append( float( line[0] ) )
      n.append( float( line[1] ) )

   file.close()
   # convert to numpy arrays
   e = NP.array(e)
   n = NP.array(n)
   return [e, n]


def calc_layered_GF(recvFilename, srcFilename, edks, OutPrefix, units2m = 1000):

   # location of binaries
   BIN_EDKS = '${EDKS_HOME}/bin'
   # load source properties
   eS, nS, depS, sS, dS , aS , idS= read_srcFile(srcFilename)

   # load receivers
   eR, nR = read_recvFile(recvFilename)

   # change units to meters
   eS = eS * units2m
   nS = nS * units2m
   depS = depS * units2m
   aS = aS * units2m * units2m
   eR = eR * units2m
   nR = nR * units2m

   # calculate the GF's
   # left lateral strike slip unit dislocation
   rake = 0 * NP.ones(eS.shape)
   slip = 1 * NP.ones(eS.shape)
   GFeSS, GFnSS, GFzSS = layered_disloc_sub(eS, nS, depS, sS, dS, rake, slip,\
                          aS, eR, nR, edks, OutPrefix, BIN_EDKS)  

   # thrust (updip) dip slip unit dislocation
   rake = 90 * NP.ones(eS.shape)
   slip = 1 * NP.ones(eS.shape)
   GFeDS, GFnDS, GFzDS = layered_disloc_sub(eS, nS, depS, sS, dS, rake, slip,\
                          aS, eR, nR, edks, OutPrefix, BIN_EDKS)         


   return [GFeSS, GFnSS, GFzSS, GFeDS, GFnDS, GFzDS]


if __name__ == '__main__':
   import sys

   usage = """
   usage: 

      %s recvFilename srcFilename edks OutPrefix units2m

   """ %(sys.argv[0])

   if len(sys.argv) != 6:
      print usage
      exit()

   scriptName = sys.argv[0]
   recvFilename = sys.argv[1]
   srcFilename = sys.argv[2]
   edks = sys.argv[3]
   OutPrefix = sys.argv[4]
   units2m = float(sys.argv[5])

   GFeSS, GFnSS, GFzSS, GFeDS, GFnDS, GFzDS = calc_layered_GF(recvFilename,\
                                      srcFilename, edks, OutPrefix, units2m)
  
   # save the matrices into text files
   NP.savetxt(OutPrefix+'GFeSS.edksGF', GFeSS)
   NP.savetxt(OutPrefix+'GFnSS.edksGF', GFnSS)
   NP.savetxt(OutPrefix+'GFzSS.edksGF', GFzSS)

   NP.savetxt(OutPrefix+'GFeDS.edksGF', GFeDS)
   NP.savetxt(OutPrefix+'GFnDS.edksGF', GFnDS)
   NP.savetxt(OutPrefix+'GFzDS.edksGF', GFzDS)

    

