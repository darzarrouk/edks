import numpy as NP
def projectGFmatrices(GeSS, GnSS, GuSS, GeDS, GnDS, GuDS, DDirE, DDirN, DDirU):
   """
   take the 6 GF matrices and compute G_SS and G_DS in which the predicted 
   displacements are projected into the direction of the receivers.

   G?SS and G?DS have dimension {Nobs}x{Npatches}
   DDir? have dimension {Nobs}

   The rows of the six input Green Function matrices and the elements of DDir? must 
   be corresponding (ie row i is related to the same observation as the i-th
   element of DDir?.

   returns the projected matrices.

   return [G_SS, G_DS]

   """

   Nobs, Npatch = GeSS.shape

   # Initialize projected GF matrices
   G_SS = NP.zeros((Nobs, Npatch))
   G_DS = NP.zeros((Nobs, Npatch))

   # do the projection

   for i in range(0, Nobs):
      G_SS[i][:] = GeSS[i][:] * DDirE[i]\
                 + GnSS[i][:] * DDirN[i]\
                 + GuSS[i][:] * DDirU[i]
 
      G_DS[i][:] = GeDS[i][:] * DDirE[i]\
                 + GnDS[i][:] * DDirN[i]\
                 + GuDS[i][:] * DDirU[i]

   return [G_SS, G_DS]
