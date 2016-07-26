#!/usr/bin/env python
"""
    by Francisco Hernan Ortega Culaciati
    Seismological Laboratory
    California Institute of Technology
    ortega@gps.caltech.edu
    
    Created      : Aug 26, 2011
    Last modified: Jul 25, 2016

    Modification History:
        -  python3 compatibility (T. Ragon)

"""
from subprocess import Popen, PIPE
import time
import numpy as NP

###
def runEDKS(params):

   output_filename = '%(ModPrefix)s_%(dep)f_%(dR)f_%(Rmax)f.asc' % params
   params['output_filename'] = output_filename
   
   # create shell command for EDKS
   # $BIN/tab4E ${MODEL} ${DEPTH} ${DELTAR} ${RMAX} ${output_filename}
   cmd = '%(BIN_tab4E)s %(ModPrefix)s.model %(dep)f ' % params
   cmd += '%(dR)f %(Rmax)f %(output_filename)s > /dev/null' % params
   #cmd += '%(dR)f %(Rmax)f %(output_filename)s >' % params
   import os
   print(cmd)
   os.system(cmd)
   
   # get the number of lines in the output_file
   params['nl'] = len(open(params['output_filename'],'r').readlines())

   #logdata = edks(params)
   logdata = '%(dep)f %(dR)f %(Rmax)f %(output_filename)s %(nl)i' % params
   return logdata
   
###
def readEDKSconfigFile(filename):
   # for the moment just return the values for a simple example
   # but later I will need to read them from a configuration file.
   data = {}
   file = open(filename,'r')
   lines = file.readlines()
   for line in lines:
      line = line.split()
      if len(line) > 0 and line[0] != '#':
         data[line[0]] = line[1:]
   # now remove from list and transform to float if necessary 
   data['BIN_tab4E'] = data['BIN_tab4E'][0]
   data['BIN_build_edks'] = data['BIN_build_edks'][0]
   data['ModPrefix'] = data['ModPrefix'][0]
   data['dR'] = float(data['dR'][0])
   data['Rmax'] = float(data['Rmax'][0])
   data['depths'] = [float(dep) for dep in data['depths']]
   
   return data


def runEDKS2(model_name):
    filename = model_name+'.edks_config'
    data = readEDKSconfigFile(filename)
    
    params = {}
    params['BIN_tab4E'] = data['BIN_tab4E']
    params['BIN_build_edks'] = data['BIN_build_edks']
    params['ModPrefix'] = data['ModPrefix']
    params['dR'] = data['dR']
    params['Rmax'] = data['Rmax']
    
    for i in data['depths']:
        params['depths'] = i
        runEDKS(params)
    
    return
