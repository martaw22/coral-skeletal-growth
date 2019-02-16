# -*- coding: utf-8 -*-
"""
Created on Wed Dec 12 16:27:17 2018

@author: Marta
"""

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import pickle

fig = plt.figure()
ax = fig.gca(projection='3d')
ax.set_aspect('equal')

u,v = np.mgrid[0:2*np.pi:20j, 0:np.pi/2:10j]
x = np.cos(u)*np.sin(v)
y = np.sin(u)*np.sin(v)
z = np.cos(v)
ax.plot_wireframe(x,y,z, color='r')



fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

# Make data
#useful website for understanding spherical coordinates: https://mathinsight.org/spherical_coordinates
#so radius determines the size of the spheres, and the angles determine how much of a sphere they are (ie half, full, other shapes, etc)
r=10                               ##radius
u = np.linspace(0, 2*np.pi, 100)   ##theta in spherical coordinates - (how far around the x and y axes in a full circle we go)
v = np.linspace(0, np.pi/2, 100)     ##phi in spherical coordinates - determines how far along z axis sphere is (so forms a full half circle from fully positive to fully negative z when at full pi units)
x = r * np.outer(np.cos(u), np.sin(v))
y = r * np.outer(np.sin(u), np.sin(v))
z = r * np.outer(np.ones(np.size(u)), np.cos(v))

# Plot the surface
ax.plot_surface(x, y, z, color='b')

plt.show()


with open('/Users/Marta/Documents/Python/Plots/modelnucleisphere_om10.png', 'wb') as pickle_file:
    ax = pickle.load((pickle_file))
plt.show()

