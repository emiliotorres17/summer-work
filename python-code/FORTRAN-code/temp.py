#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt

x           = np.linspace(0.0, 1.0, 33)
[X1, X2]    = np.meshgrid(x,x)                             
v           = np.loadtxt('../data/v-velocity.dat')[1:34,:]
u           = np.loadtxt('../data/u-velocity.dat')[:,1:34]
omega       = np.loadtxt('../data/vorticity.dat')
umag        = np.sqrt(v**2.0+u**2.0)                                
cnt         = plt.quiver(X1, X2, np.transpose(u), np.transpose(v))#np.arange(-15.0, 0.0, 15.0/500), cmap='jet', extend='both')
#`for c in cnt.collections:
#`    c.set_edgecolor('face')
#`plt.colorbar()
plt.show()
