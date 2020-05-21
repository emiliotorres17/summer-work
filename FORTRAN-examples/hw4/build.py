#!/usr/bin/env python3
from subprocess import call
import numpy as np
import matplotlib.pyplot as plt
call(['clear'])
u           = np.loadtxt('data.dat')
x           = np.loadtxt('x.dat')
y           = np.loadtxt('y.dat')
[X1, X2]    = np.meshgrid(x,y)
c1          = 1.0
upper       = c1
lower       = 0.0
dp          = (upper-lower)/500.0
cnt         = plt.contourf(X1, X2, np.transpose(u), np.arange(lower, upper, dp),\
                    cmap='jet', extend='both')
for c in cnt.collections:
    c.set_edgecolor('face')
plt.colorbar()
plt.savefig('u-velocity.pdf')
plt.clf()
