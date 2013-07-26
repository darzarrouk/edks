import os
import struct # for reading and writing binary files.
import numpy as NP
def layered_disloc_sub(IDs, xs, ys, zs, strike, dip, rake, slip, A, xr, yr, edks,\
                   prefix, BIN_EDKS = '${EDKS_BIN}'):
   """
   [ux, uy, uz] = layered_disloc(xs, ys, zs, strike, dip, rake, slip,
                                 A, xr, yr, edks)

   
    --- INPUT ---
    --- SOURCE INFO
    --- 1D NUMPY arrays, length = number of fault patches
    IDs      list of strings,IDs of  point sources (see below for a detailed explanation)
    xs       m, east coord to center of fault patch
    ys       m, north coord to center of fault patch
    zs       m,depth coord to center of fault patch (+ down) 
    strike   deg, clockwise from north 
    dip      deg, 90 is vertical 
    rake     deg, 0 left lateral strike slip, 90 up-dip slip 
    slip     m, slip in the rake direction
    A        m2, surface area of fault patch 
    --- RECEIVER INFO
    1D arrays, length = number of receivers
    xr       m, east coordinate of receivers 
    yr       m, north coordinate of receivers 
    --- ELASTIC STRUCTURE INFO
    edks     string, full name of edks file, e.g., halfspace.edks
    --- FILE NAMING 
    prefix   string, prefix for the IO files generated by sum_layered
    BIN_EDKS string, folder (full path) where EDKS binaries are located. 
    --- OUTPUT ---
    --- 2D arrays, (receivers, fault patches)
    ux     m, east displacement
    uy     m, west displacement
    uz     m, up displacement (+ up)

   Explanation of IDs:
      the source ID (IDs) is used to be able to represent a finite source by a set
      of "well defined" point sources (ex: point sources modeling a triangular or 
      rectangular finite source). 
      - If you want to use this code only to calculate independent point sources, just
      give a different ID to all the sources.
      - If you want to use this code to approximate several finite dislocations, you need
      to define and assign a different ID to each finite source. The sources with 
      equal IDs will be added to compute the surface displacements of the finite 
      dislocation. Then the code will return only the displacements corresponding to 
      the one of the finite dislocation, in the same order as the specified IDs.
      IMPORTANT: The equal IDs must be contiguous in order to ensure that the order
      in which the output is computed is the same. 
      Ex: -  a good list of source IDs is 
             [id1, id1,..., id1, id2, id2,..., id2, idj,..., idj, ..., idN,..,idN] 
          - a BAD list of source (you should not do this) would be:
             [id1,id2, id3, ... , idN, id1, id2, .. idN, id1, id3, id8] 

   NOTE ON THE NUMBER OF CORES USED BY sum_layered_sub:
       in order to set the number of cores used (by default openMP uses all the 
       available cores) you must set the environment variable OMP_NUM_THREADS
       with the number of cores (threads) that you want. 

   """
   # some values I need to define for compatibility with "sum_layered""
   # Since sum_layered was originally thought for modelling rectangular patches.
   nrec = len(xr) # number of receivers
   np = len(set(IDs)) # total number of finite faults 
   ntsp = len(xs) # total number of point sources (sub sources)
   # compute the mapping between the string IDs and non decreasing positive integer number
   setOfAlreadyStoredIDs = set() # this is just for testing the right order of the IDs.
   sortedListOfFiniteFaultIDs = []
   ident = []
   i = 1
   IDprev = IDs[0]
   ident.append(i)
   sortedListOfFiniteFaultIDs.append(IDprev)
   setOfAlreadyStoredIDs.add(IDprev)
   NumSubSources = {}
   NumSubSources[IDprev] = 1
   for k in range(1, ntsp):
      if IDs[k] == IDprev:
         ident.append(i)
         NumSubSources[IDprev] += 1
          
      else: # the current ID is a new one.
         if IDs[k] in setOfAlreadyStoredIDs: # this is an error
            raise ValueError('The source IDs are not in the right order, see the help')
         else:
            IDprev = IDs[k]
            i += 1
            ident.append(i)
            sortedListOfFiniteFaultIDs.append(IDs[k])
            setOfAlreadyStoredIDs.add(IDs[k])
            NumSubSources[IDs[k]] = 1
   nspp = NP.max(NumSubSources.values()) # maximum number of subsources

   BIN_FILE_FMT_real4 = 'f' # python float = C/C++ float = Fortran 'real*4' 
   BIN_FILE_FMT_int4 = 'i' # python int, fortran 'integer*4''
   NBYTES_FILE_FMT = 4  # a Fortran (real*4) uses 4 bytes.

   # Define filenames:
   file_rec = prefix + '.rec'
   file_pat = prefix + '.pat'
   file_dux = prefix + '_ux.dis'
   file_duy = prefix + '_uy.dis'
   file_duz = prefix + '_uz.dis'
   cmd = 'rm -f %s %s %s %s %s' %(file_rec, file_pat, file_dux, file_duy, file_duz)
   os.system(cmd) # delete the files if exist
   
   # write receiver location file (observation points)
   temp = [xr, yr]
   file = open(file_rec, 'wb') 
    
   for k in range(0, nrec):
      for i in range(0, len(temp)):
         file.write( struct.pack( BIN_FILE_FMT_real4, temp[i][k] ) )       
   file.close() 
 
   # write point sources information
   temp = [xs, ys, zs, strike, dip, rake, A, slip, ident]
   file = open(file_pat, 'wb');
   for k in range(0, ntsp):
      for i in range(0, len(temp)-1):
         file.write( struct.pack( BIN_FILE_FMT_real4, temp[i][k] ) )
      file.write( struct.pack( BIN_FILE_FMT_int4, temp[-1][k] ) )
      
   file.close()
 
   # call sum_layered
   print "before sum_layered"
   cmd = '%s/sum_layered_sub %s %s %i %i %i %i'%(BIN_EDKS, \
                           edks, prefix, nrec, np, ntsp, nspp)
   print cmd
   os.system(cmd)
   print "after sum_layered"

   # read sum_layered output Greens function
   # ux
   file = open(file_dux, 'rb')
   ux = NP.zeros((nrec, np))
   for j in range(0, np):
      for i in range(0, nrec):
         byteVal = file.read(NBYTES_FILE_FMT)
         if byteVal != '':
            ux[i][j] = struct.unpack('f', byteVal)[0]
         else:
            raise ValueError(' Premature EOF in %s, something nasty happened'%(file_dux))
   file.close()

   # uy
   file = open(file_duy, 'rb')
   uy = NP.zeros((nrec, np))
   for j in range(0, np):
      for i in range(0, nrec):
         byteVal = file.read(NBYTES_FILE_FMT)
         if byteVal != '':
            uy[i][j] = struct.unpack('f', byteVal)[0]
         else:
            raise ValueError('Premature EOF in %s, something nasty happened'%(file_duy))
   file.close()

   # uz
   file = open(file_duz, 'rb')
   uz = NP.zeros((nrec, np))
   for j in range(0, np):
      for i in range(0, nrec):
         byteVal = file.read(NBYTES_FILE_FMT)
         if byteVal != '':
            uz[i][j] = struct.unpack('f', byteVal)[0]
         else:
            raise ValueError('Premature EOF in %s, something nasty happened'%(file_duz))
   file.close()

   # remove IO files.
   cmd = 'rm -f %s %s %s %s %s' %(file_rec, file_pat, file_dux, file_duy, file_duz)
   os.system(cmd)  

   # return the computed displacements for each sources
   return [ux, uy, uz]


