# File with the triangles properties
SquaresPropFile = './Gocad2009/Mainshock_TohokuOki_V3.ts.3DsquareGrid_ShiftedUp2.8KM.END'
# File with id, E[km], N[km] coordinates of the receivers.
ReceiverFile = './VecObs/ITO/ITO.EDKSsub'
# read receiver direction (# not yet implemented)
useRecvDir = True # if True it will require the observation direction in Receiver File. 
# Maximum Area to subdivide triangles. If None, uses Saint-Venant's principle.
Amax = None # None computes Amax automatically. float, uses this value for all triangles
EDKSunits = 1000.0 # to convert from kilometers to meters
EDKSfilename = 'tohoku_nied_2_shifted2.8km.edks'
prefix = 'Mainshock_3DS_Tohoku_SATO_'
plotGeometry = True # set to False if you are running in a remote Workstation
