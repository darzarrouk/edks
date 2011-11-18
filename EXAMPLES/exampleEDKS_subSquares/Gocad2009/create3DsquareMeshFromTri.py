#!/usr/bin/env python
import numpy as NP
import string
import pylab as PL

# read the data with the properties of the triangles
filename = 'Mainshock_TohokuOki_V3.ts.TriangleProp.txt'
file = open(filename, 'r')
# skip header
file.readline()
Tid = []
e = []
n = []
d = []
phi = []
dip = []
# now read all the elements.
for line in file:
   line = string.split(line)
   Tid.append(line[0])
   e.append(float(line[3]))
   n.append(float(line[4]))
   d.append(float(line[5]))
   phi.append(float(line[6]))
   dip.append(float(line[7])) 
file.close()
e = NP.array(e)
n = NP.array(n)
d = NP.array(d)
phi = NP.array(phi)
dip = NP.array(dip)

# median strike
Mphi = NP.median(phi)
print "Median strike angle = ", Mphi
PL.hist(phi)
PL.title('Strike Histogram')

# Compute the center of the reference frame to perform the coordinate system rotation
Me = NP.median(e)
Mn = NP.median(n)

# translate the coordinates
te = e - Me
tn = n - Mn
td = d - 0.0

# make the rotation
Mphi_rad = Mphi * NP.pi / 180.0
tp = te * NP.cos(Mphi_rad) - tn * NP.sin(Mphi_rad)
ts = te * NP.sin(Mphi_rad) + tn * NP.cos(Mphi_rad)

PL.figure()
PL.plot(te, tn,'.r',tp, ts, '.k')
#PL.axis('equal')

# fit a polinomial function to the set of data.
deg = 4 # degree of the polinomial function
# The polinomial is td = P(tp)
# create the design matrix
G = NP.zeros((len(tp), deg + 1))
# constant
for i in range(0,len(tp)):
   G[i][0] = 1
# powers of polinomials
for i in range(0,len(tp)):
   for j in range(1,deg+1):
      G[i][j] = (tp[i])**(j)

# get the polinomial coefficients
c = NP.dot(G.transpose(), td)
print c.shape, (NP.linalg.inv(NP.dot(G.transpose(), G))).shape, c.shape
c = NP.dot(NP.linalg.inv(NP.dot(G.transpose(), G)), c)

tp_p = NP.linspace(tp.min(), tp.max(), 100000)
G_p = NP.zeros((len(tp_p), deg + 1))
# constant
for i in range(0,len(tp_p)):
   G_p[i][0] = 1
# powers of polinomials
for i in range(0,len(tp_p)):
   for j in range(1,deg+1):
      G_p[i][j] = (tp_p[i])**(j)

td_p = NP.dot(G_p, c)


# compute the total length of the polinomial "Width of the surface".
dL = NP.sqrt((tp_p[1:] - tp_p[:-1])**2 + (td_p[1:] - td_p[:-1])**2)
Lp = [0.0]
Lp.extend(NP.cumsum(dL))

# Total along strike length of the fault
Ls = NP.max(ts) - NP.min(ts)

# Define length of side of square
side = 29.0  # in [km]

# compute the location of the edges in the profile
#Pedges = []
#Dedges = []
#Pedges.append(tp_p[0])
#Dedges.append(td_p[0])
#lastI = 0
#for i in range(0,len(Lp)):
#   if Lp[i] - Lp[lastI] > side:
#      lastI = i
#      Pedges.append(tp_p[i])
#      Dedges.append(td_p[i])
#Pedges = NP.array(Pedges)
#Dedges = NP.array(Dedges)
# compute the location of the edges in the profile
Pedges = []
Dedges = []
Pedges.append(tp_p[0])
Dedges.append(td_p[0])
for i in range(0,len(tp_p)):
   dist = NP.sqrt((Pedges[-1] - tp_p[i])**2 + (Dedges[-1] - td_p[i])**2)
   if dist > side: # this assumes that tp_p spacing is very dense (1/10000 in this case)
      Pedges.append(tp_p[i])
      Dedges.append(td_p[i])
Pedges = NP.array(Pedges)
Dedges = NP.array(Dedges)



Pcenters = 0.5 * (Pedges[1:] + Pedges[:-1])
Dcenters = 0.5 * (Dedges[1:] + Dedges[:-1])
Dip = (180.0/NP.pi) * NP.arctan((Dedges[1:] - Dedges[:-1])/(Pedges[1:] - Pedges[:-1]))
NsquaresPerProf = len(Pcenters)
print "NsquaresPerProf = ", NsquaresPerProf

# plot
PL.figure()
PL.plot(tp, -td, '.k', tp_p, -td_p, '-r', Pcenters, -Dcenters, 'or',Pedges, -Dedges,'^r')
PL.grid('on')
PL.axis('equal')

# extrude the centers up to max and min ts
Sedges = []
Sedges.extend(NP.arange(-side/2.0, NP.min(ts) - 4*side, -side)) # added 3side to extend itto the north from the previous one
Sedges.sort()
Sedges.extend(NP.arange(side/2.0, NP.max(ts), side))
Sedges = NP.array(Sedges)
Scenters = 0.5 * (Sedges[1:] + Sedges[:-1])
NumProfiles = Scenters.size
print "NumProfiles = ", NumProfiles
print 'Number of "squares" = ', NumProfiles*NsquaresPerProf
# compute coordinates and properties for all grid points
S = []
P = []
D = []
Gstrike = []
Gdip = []
Garea = []
for iS in range(0, len(Scenters)):
   for iP in range(0, len(Pcenters)):
      S.append(Scenters[iS])
      P.append(Pcenters[iP])
      D.append(Dcenters[iP])
      Gstrike.append(Mphi)
      Gdip.append(Dip[iP])
      Garea.append(side*side)
Gs = NP.array(S)
Gp = NP.array(P)
Gd = NP.array(D)
Gstrike = NP.array(Gstrike)
Gdip = NP.array(Gdip)
Garea = NP.array(Garea)

# rotate and translate back the grid.
Ge = Gp * NP.cos(Mphi_rad) + Gs * NP.sin(Mphi_rad)
Gn = - Gp * NP.sin(Mphi_rad) + Gs * NP.cos(Mphi_rad)
Ge = Ge + Me
Gn = Gn + Mn
import mpl_toolkits.mplot3d.axes3d as p3
fig = PL.figure()
plot3 = p3.Axes3D(fig)
plot3.scatter3D(e, n, -d, color = 'red',)
plot3.scatter3D(Ge, Gn, -Gd, color = 'blue', )
#plot3.axis('equal')
#PL.axis('equal')


# save the results.
filename = filename + '.3DsquareGrid'
file = open(filename, 'w')
header = 'E[km] N[km] Dep[km] strike dip Area ID\n'
file.write(header)
for i in range(0, len(Ge)):
   line='% 12.3f % 12.3f % 7.2f % 6.2f % 6.2f % 6.1f %05i\n'%(Ge[i],Gn[i],Gd[i],Gstrike[i],Gdip[i],Garea[i],i)
   file.write(line)
file.close() 




PL.show()





