# File with the triangles properties
SquaresPropFile = './Gocad2009/Mainshock_TohokuOki_V3.ts.3DsquareGrid.END'
# File with id, E[km], N[km] coordinates of the receivers.
ReceiverFile = './VecObs/JapanGPS/JapanGPSObs.idEN'
# read receiver direction (# not yet implemented)
useRecvDir = False # leave this to False for the moment (will be available in the future)
# Maximum Area to subdivide triangles. If None, uses Saint-Venant's principle.
Amax = None # None computes Amax automatically. float, uses this value for all triangles
EDKSunits = 1000.0 # to convert from kilometers to meters
EDKSfilename = 'tohoku_3.edks'
prefix = 'Mainshock_3DS_Tohoku_'
plotGeometry = True # set to False if you are running in a remote Workstation
