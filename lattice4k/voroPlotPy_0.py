#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 12 16:55:33 2023

@author: hossein
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial import Voronoi, voronoi_plot_2d

rng = np.random.default_rng()
# points = rng.random((10,2))
# points = np.random.random((10,2))
points = np.loadtxt("test.dat", delimiter=' ')

LxLy = np.loadtxt("LxLy.dat", dtype=float, delimiter=' ')
L = LxLy[0]

total_points = points.copy()

# copy_direction_list = [[-1,-1], [-1, 0], [-1,+1], [0, +1], [+1, +1], [+1, 0], [+1, -1], [0,-1] ]
copy_direction_list = []

for copy_counter in range(0,len(copy_direction_list)):
    copy_dir = copy_direction_list[copy_counter]
    
    copy_points = points.copy()
    
    copy_points[:,0] += L*copy_dir[0]
    copy_points[:,1] += L*copy_dir[1]
    
    total_points = np.concatenate((total_points, copy_points), axis=0)

total_points_reversed = 0.0*total_points
total_points_reversed[:,0]=total_points[:,1].copy()
total_points_reversed[:,1]=total_points[:,0].copy()

vor = Voronoi(total_points_reversed)
# fig = voronoi_plot_2d(vor)
# fig = voronoi_plot_2d(vor, show_vertices=False, line_colors='orange',
#                       line_width=2, line_alpha=0.6, point_size=2)

fig = plt.figure()

# Create Voronoi plot
voronoi_plot_2d(vor, show_vertices=False, line_colors='orange',
                line_width=0.5, line_alpha=0.6, point_size=0.5)

# ax = fig.gca()
# plt.xlim(-L,2*L)
# plt.ylim(-L,2*L)
# ax.invert_yaxis()

# ax = fig.add_subplot(1, 1, 1)
# Set x and y limits
# L = 10  # Example limit value
# plt.xlim(-L, 2*L)
# plt.ylim(-L, 2*L)
plt.xlim(0, L)
plt.ylim(0, L)
plt.gca().invert_yaxis()

plt.plot([0,0],[0,L],linestyle='dashed', linewidth=0.6, color='black')
plt.plot([0,L],[0,0],linestyle='dashed', linewidth=0.6, color='black')
plt.plot([0,L],[L,L],linestyle='dashed', linewidth=0.6, color='black')
plt.plot([L,L],[0,L],linestyle='dashed', linewidth=0.6, color='black')
# # plt.show()
for i in range(len(points)):
    x = points[i,0]
    y = points[i,1]
    plt.text(y, x, i, fontsize=1)


# # ax.invert_yaxis()
plt.xlabel('Y')
plt.ylabel('X')
plt.grid(which='both')
# plt.axis("equal")
plt.gca().set_aspect("equal")
# # ax.set_xlim((-L,2*L))
# # ax.set_xlim(-L,2*L)
plt.savefig("voro_test.png", dpi=500)
# plt.show()
