#!/usr/bin/env python
"""
    by Francisco Hernan Ortega Culaciati
    Seismological Laboratory
    California Institute of Technology
    ortega@gps.caltech.edu
    
    Created      : Aug 26, 2011
    Last modified: Jul 25, 2016 

    Modification History:
    
       - Python3 compatibility (T. Ragon)

"""

from LogFile import LogFile
from os import system
from mpi4py import MPI
import edks
from numpy import argsort

# Define function for worker process in MPI
def EDKSworker(comm, TagIn, TagOut):
   """
   This function executes EDKS given the parameters in data.

   rank of boss process is assumed to be 0

   comm: the MPI communicator to use.
   TagIn: MPI_TAG (int) for receiving data from Boss process.
   TagOut: MPI_TAG (int) for sending data to Boss process.

   It expects a message (comm.send) from Boss process that contains a 
   dictionary with the following keys:
      - 'ActiveProc' : Boolean variable that determines if the Worker is going
		       to perform EDKS calculations (True) or if the Boss is 
		       terminating the Worker process (False)

      - 'EDKSparam' : data structure to pass to EDKS.

   """
   rank = comm.Get_rank()
   
   # Notify boss that worker is ready to begin work
   comm.send( {'rank': rank }, dest = 0, tag = TagOut )
   # receive data structure from Boss process
   data = comm.recv(source = 0, tag = TagIn)
   while data['ActiveProc']:
      # compute EDKS
      result = edks.runEDKS(data['EDKSparam']) 
      # send the results and rank to Boss
      ToSend = {'rank' : rank, 'result' : result}
      comm.send( ToSend, dest = 0, tag = TagOut )
      # wait for new data and instructions from Boss
      data = comm.recv(source = 0, tag = TagIn) 
      
   # If i Reached this point the Boss has terminated this process
   return data['ActiveProc'] # this should be always False


# define function for Boss process in MPI
# Boss is in charge of distributing a balanced load among processors.
def EDKSboss(comm, EDKSparam, TagIn, TagOut):
   """
   This functions distributes a balanced workload among the available
   processors. 
   """
   rank = comm.Get_rank()
   size = comm.Get_size()

   # the partitioning is done by depths.
   depths = EDKSparam.pop('depths')    
   depths.sort(reverse = True) # So depths.pop() goes from small to high depth
   # check for errors
   if rank != 0:
      raise ValueError("EDKSboss can only be run in root process.")
      return None

   Logger = EDKSparam['Logger']
   Logger.addLine('I am Boss with rank = %i, I have %i Workers'%(rank, size - 1))

   # create dictionary to keep track of active processes
   ranks = range(1,size)
   flags = [ True for elem in ranks ]
   ActiveProcs = dict(zip(ranks,flags))
   
   # to store results
   results = []
   
   # iterate over depths and send to processes reporting results
   print(ActiveProcs.values())
   #while (ActiveProcs.values()).count(True) > 0 :
   while True in ActiveProcs.values():
      # listen for a process reporting results
      # for the moment i dont need to do anything with results.
      dataRECV = comm.recv( source = MPI.ANY_SOURCE, tag = TagIn )
      pRank = dataRECV['rank']
      if 'result' in dataRECV:
         results.append( dataRECV['result'] )         
         Logger.addLine('received results from Worker %i'%(pRank))
      else:
         Logger.addLine("Worker %i Initialized"%(pRank))
         
      # send next depth to the process that just reported
      if len(depths) > 0:
         dep = depths.pop()
         EDKSparam['dep'] = dep
         ToSend = {'ActiveProc': True, 'EDKSparam': EDKSparam}
         msg = "Boss sending depth %.3f to worker %i ...\n"%(dep, pRank)
         msg += "    %i depths waiting for workers."%(len(depths))
         Logger.addLine(msg)
         
      else:
         ToSend = {'ActiveProc': False} #shuts down the worker
         ActiveProcs[pRank] = False
         Logger.addLine("Worker %i being Terminated"%(pRank))

      comm.send( ToSend, dest = pRank, tag = TagOut ) 
      
   # return results
   return results



###
def MPI_EDKS(comm, EDKSparam, TagBoss2Worker = 1979, TagWorker2Boss = 28 ):

   rank = comm.Get_rank()
   size = comm.Get_size()

   if rank == 0 : # I am boss process
      Logger = LogFile('EDKS_status.log') # to keep track of what it is doing 
      Logger.clearFile()
      if 'Logger' not in EDKSparam:
         EDKSparam['Logger'] = Logger
      else:
         raise ValueError("EDKS parameter name 'Logger' is a reserved name")
         return None
      # call the Boss process
      ndep = len(EDKSparam['depths'])
      ModPref = EDKSparam['ModPrefix']
      # run EDKS calculations and save the EDKS log data
      Logger.addLine('Starting EDKS kernels calculation')
      logData = EDKSboss( comm, EDKSparam, TagWorker2Boss, TagBoss2Worker )
      Logger.addLine('EDKS kernels calculation done... Now saving EDKS log data.')
      logFilename = ModPref + '.log'
      logfile = open(logFilename, 'w')
      logfile.write(ModPref + '.model\n')
      logfile.write("%d\n"%ndep)
      depths = [float(line.split()[0]) for line in logData]
      Idep = argsort(depths)
      for i in Idep:
         logfile.write('%f '%depths[i])
      logfile.write('\n')
      for i in Idep:
         logfile.write(logData[i] + '\n')
      logfile.close()
      Logger.addLine('running build_edks...')
      edksFilename = ModPref + '.edks'
      system('rm -rf %s hdr.%s'%(edksFilename, edksFilename))
      system(EDKSparam['BIN_build_edks'] + ' ' + logFilename + ' ' + edksFilename)
      Logger.addLine("That's all folks !!!...")
   else : # I am a worker

      EDKSworker( comm, TagBoss2Worker, TagWorker2Boss )


###
def Main(configFilename):

   # Initialize MPI communicator
   comm = MPI.COMM_WORLD
   rank = comm.Get_rank()
   EDKSconfig = None
   if rank == 0: # only Boss job needs the EDKSconfig data.
      # read configuration file   
      EDKSconfig = edks.readEDKSconfigFile(configFilename)
   
   if comm.Get_size() > 1:
      MPI_EDKS(comm, EDKSconfig)
   else:
      msg = "number of processors must be > 1"
      raise ValueError(msg)  

if __name__ == '__main__':
   import sys
   if len(sys.argv) == 2:
      configFilename = sys.argv[1]
      Main(configFilename)

   else:
      print("  usage: \n")
      print("        " + sys.argv[0] + "  configFile.edks_config")
      exit()




