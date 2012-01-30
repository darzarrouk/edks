# File with the triangles properties
SquaresPropFile = './3D_squares.END'
# File with id, E[km], N[km] coordinates of the receivers.
ReceiverFile = './receivers.idEN'
# read receiver direction (# not yet implemented)
useRecvDir = False # leave this to False for the moment (will be available in the future)
# Maximum Area to subdivide triangles. If None, uses Saint-Venant's principle.
Amax = None # None computes Amax automatically. float, uses this value for all triangles
EDKSunits = 1000.0 # to convert from kilometers to meters
EDKSfilename = 'poisson_half_space.edks'
prefix = 'EDKS_HalfSpace_Poisson_'
plotGeometry = True # set to False if you are running in a remote Workstation
