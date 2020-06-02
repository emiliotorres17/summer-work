#!/usr/bin/env python3
import sys
import numpy as np
import matplotlib.pyplot as plt

data    = np.loadtxt('vel-copy.txt')
print(data.shape)
u       = data[0:32, 0:32]
v       = data[32:64, 0:32]
print(data.shape)
print(u.shape)
print(v.shape)
dx      = 1.0/32
omega   = np.gradient(v, dx, edge_order=2)[1] -\
            np.gradient(u, dx, edge_order=2)[0]
#---------------------------------------------------------------------#
# Domain variables                                                    #
#---------------------------------------------------------------------#
dim     = omega.shape
x       = np.linspace(0.0, 1.0, dim[0])
y       = np.linspace(0.0, 1.0, dim[1])
[X,Y]   = np.meshgrid(x,y)
#---------------------------------------------------------------------#
# Plot settings                                                       #
#---------------------------------------------------------------------#
lower   = 0.0
upper   = 1.0
dp      = (upper-lower)/10.0
con     = [0.0, -0.5, 0.5, -1.0, 1.0, -2.0, 2.0, -3.0, 3.0, -4.0, -5.0]
con     = np.array([-5.0, -4.0, -3.0, -2.0, -1.0, -0.5, 0.0, 0.5,\
                    1.0, 2.0, 3.0])
#con = np.arange(0.0, 1.0, 1.0/500.0)
#---------------------------------------------------------------------#
# Generating contour plots                                            #
#---------------------------------------------------------------------#
cnt     = plt.contour(X,Y, omega, con,\
                        cmap='jet', extend='both')
for c in cnt.collections:
    c.set_edgecolor('face')
plt.colorbar()
#plt.savefig(media_path + 'vorticity-' + str(dim[0]-1) + '.png')
plt.show()


