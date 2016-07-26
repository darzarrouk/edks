#!/usr/bin/env python
import numpy as NP
###
def readLayerThicknessFromVelModel(modelFilename):
   file = open( modelFilename, 'r' )
   lines = file.readlines()
   file.close()
   # first line contains number of layers 
   line = lines[0].split()
   Nlayers = int(line[0])
   UnitsFactor = float(line[1])
   # now read layer thickness
   LTh = [float( ((lines[i]).split())[3] ) for i in range(1, 1 + Nlayers)]
   return LTh

###
def createEDKSconfigFile(modelFilename):
   # get the model prefix from modelFilename
   I = modelFilename.find('.model')
   ModelPrefix = modelFilename[:I]
   # get the layer Thickness
   LayerThick = NP.array(readLayerThicknessFromVelModel(modelFilename))
   DepDown = NP.cumsum(LayerThick)
   DepUp = DepDown - LayerThick
   
   MaxDepth = eval(raw_input("Max depth for sources (km)? "))
   
   # go through each layer and compute the depths of the layer:
   SourceDepths = []
   for i in range(0, len(LayerThick) -1): # we dont want the halfspace yet.
      dH = LayerThick[i]
      msg = 'Layer %i with thickness %.3f, please indicate Number of sources (>1):'\
             %(i+1, dH)
      Ndepths = NP.array(eval(raw_input(msg)))
      #Ndepths = NP.max( [NP.ceil( 1.0*dH/dZ ) + 1. , 2.] )
      dZeff = (1.0*dH) / Ndepths  
      print( '   --> depth spacing for layer is %.3f [km] ...'%(dZeff))
      depths = NP.linspace(DepUp[i] + dZeff/2.0,\
                           DepDown[i] - dZeff/2.0, Ndepths - 1)
      SourceDepths.extend(depths)

   if MaxDepth > NP.max(DepDown):
      dH = MaxDepth - NP.max(DepDown)
      msg = 'Layer %i with thickness %.3f, (Half space)'%(len(LayerThick), dH)\
             + ' please indicate Number of sources (>1):'
      Ndepths = NP.array(eval(raw_input(msg)))
      #Ndepths = NP.max( [NP.ceil( 1.0*dH/dZ ) + 1. , 2.] )
      dZeff = (1.0*dH) / Ndepths
      depths = NP.linspace(NP.max(DepDown) + dZeff/2.0,\
                           MaxDepth - dZeff/2.0, Ndepths - 1)
      SourceDepths.extend(depths)
   SourceDepths.extend([MaxDepth])
   # now only select depths shallower than MaxDepth
   SourceDepths = NP.array(SourceDepths)
   test = SourceDepths <= MaxDepth
   SourceDepths = SourceDepths[test]
   print( SourceDepths , len(SourceDepths) )

   # define the rest of the variables
   Rmax = eval(raw_input("Maximum horizontal source - receiver distance (km)? "))  
   dR = eval(raw_input("Step for Horizontal  source - receiver distance (km)? ")) 
   BIN_DIR = '${EDKS_HOME}/bin'
   ans = raw_input("Folder with EDKS binaries [" + BIN_DIR + "] ?")
   if len(ans) > 0:
      BIN_DIR = ans
   BIN_tab4E = BIN_DIR + '/tab4E'
   BIN_build_edks = BIN_DIR + '/build_edks'
   
   # write the configuration file
   configFilename = ModelPrefix + '.edks_config' 
   cfile = open(configFilename, 'w')
   text = """# This is a table that defines the input required for EDKS
# any line whose first character is '#' will be ignored.
# The fields required are:
# BIN_tab4E : full path of the tab4E executable.
# BIN_build_edks: full path of the build_edks executable.
# ModPrefix : Name of the velocity model file prefix 
#             if velocity model file is named "tohoku_3.model"
#             then ModPrefix is "tohoku_3"
# dR : step for calculating the horizontal distance of each station (km).
# Rmax : Maximum horizontal distance for station (km).
# depths : list of source depths (km)
# 
# Here follows an example of the file """
   cfile.write(text + '\n')
   cfile.write('BIN_tab4E        ' + BIN_tab4E + '\n')
   cfile.write('BIN_build_edks   ' + BIN_build_edks + '\n')
   cfile.write('ModPrefix        ' + ModelPrefix + '\n')
   cfile.write('dR               ' + str(dR) + '\n')
   cfile.write('Rmax             ' + str(Rmax) + '\n')
   cfile.write('depths          ')
   for dep in SourceDepths:
      cfile.write(' % .3f'%(dep) )

   cfile.close()
   print ('\n \n **********************************************\n \n')
   print ('EDKS configuration file "' + configFilename + '" has been created...')
   print ('Number of depths : ' + str(len(SourceDepths)))


if __name__ == '__main__':
   import sys
   if len(sys.argv) == 2:
      vel_model_Filename = sys.argv[1]
      createEDKSconfigFile(vel_model_Filename)

   else:
      print( "  usage: \n")
      print( "        " + sys.argv[0] + "  edksVelocityModel.model")
      exit()


###
def createEDKSconfigFile2(modelFilename, MaxDepth, Nd, Rmax, dR):
   # get the model prefix from modelFilename
   I = modelFilename.find('.model')
   ModelPrefix = modelFilename[:I]
   # get the layer Thickness
   LayerThick = NP.array(readLayerThicknessFromVelModel(modelFilename))
   DepDown = NP.cumsum(LayerThick)
   DepUp = DepDown - LayerThick
   
     
   # go through each layer and compute the depths of the layer:
   SourceDepths = []
   for i in range(0, len(LayerThick) -1): # we dont want the halfspace yet.
      dH = LayerThick[i]
      #msg = 'Layer %i with thickness %.3f, please indicate Number of sources (>1):'\
       #       %(i+1, dH)
      Ndepths = NP.array(Nd[i])
      #Ndepths = NP.max( [NP.ceil( 1.0*dH/dZ ) + 1. , 2.] )
      dZeff = (1.0*dH) / Ndepths  
      #print '   --> depth spacing for layer is %.3f [km] ...'%(dZeff)
      depths = NP.linspace(DepUp[i] + dZeff/2.0,\
                           DepDown[i] - dZeff/2.0, Ndepths - 1)
      SourceDepths.extend(depths)

   if MaxDepth > NP.max(DepDown):
      dH = MaxDepth - NP.max(DepDown)
      #msg = 'Layer %i with thickness %.3f, (Half space)'%(len(LayerThick), dH)\
      #       + ' please indicate Number of sources (>1):'
      Ndepths = NP.array(Nd[len(Nd)-1])
      #Ndepths = NP.max( [NP.ceil( 1.0*dH/dZ ) + 1. , 2.] )
      dZeff = (1.0*dH) / Ndepths
      depths = NP.linspace(NP.max(DepDown) + dZeff/2.0,\
                           MaxDepth - dZeff/2.0, Ndepths - 1)
      SourceDepths.extend(depths)
   SourceDepths.extend([MaxDepth])
   # now only select depths shallower than MaxDepth
   SourceDepths = NP.array(SourceDepths)
   test = SourceDepths <= MaxDepth
   SourceDepths = SourceDepths[test]
   print( SourceDepths , len(SourceDepths) )

   # define the rest of the variables
   BIN_DIR = '${EDKS_HOME}/bin'
#   ans = raw_input("Folder with EDKS binaries [" + BIN_DIR + "] ?")
#   if len(ans) > 0:
#      BIN_DIR = ans
   BIN_tab4E = BIN_DIR + '/tab4E'
   BIN_build_edks = BIN_DIR + '/build_edks'
   
   # write the configuration file
   configFilename = ModelPrefix + '.edks_config' 
   cfile = open(configFilename, 'w')
   text = """# This is a table that defines the input required for EDKS
# any line whose first character is '#' will be ignored.
# The fields required are:
# BIN_tab4E : full path of the tab4E executable.
# BIN_build_edks: full path of the build_edks executable.
# ModPrefix : Name of the velocity model file prefix 
#             if velocity model file is named "tohoku_3.model"
#             then ModPrefix is "tohoku_3"
# dR : step for calculating the horizontal distance of each station (km).
# Rmax : Maximum horizontal distance for station (km).
# depths : list of source depths (km)
# 
# Here follows an example of the file """
   cfile.write(text + '\n')
   cfile.write('BIN_tab4E        ' + BIN_tab4E + '\n')
   cfile.write('BIN_build_edks   ' + BIN_build_edks + '\n')
   cfile.write('ModPrefix        ' + ModelPrefix + '\n')
   cfile.write('dR               ' + str(dR) + '\n')
   cfile.write('Rmax             ' + str(Rmax) + '\n')
   cfile.write('depths          ')
   for dep in SourceDepths:
      cfile.write(' % .3f'%(dep) )

   cfile.close()
   print ('\n \n **********************************************\n \n')
   print ('EDKS configuration file "' + configFilename + '" has been created...')
   print ('Number of depths : ' + str(len(SourceDepths)))


    
    
    
    









