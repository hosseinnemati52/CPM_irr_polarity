#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 15 16:18:52 2023

@author: hossein
"""
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
import numpy as np
from matplotlib.collections import PatchCollection
from scipy.spatial import ConvexHull
import math

# vertices_file = "vertices.dat"
# com_file = "test.dat"

# with open(com_file, 'r') as file:
#     text = file.read()
#     print(text)
# L=200
# L = float(input("L: "))
# NumCells = int(input("NumCells: "))
com_file = np.loadtxt("test.dat", dtype=float, delimiter=' ')

LxLy = np.loadtxt("LxLy.dat", dtype=float, delimiter=' ')
Lx = LxLy[0]
Ly = LxLy[1]
NumCells = np.shape(com_file)[0]

# vertices_file = np.loadtxt("vertices.dat", dtype=float, delimiter=' ')

def graham_scan(points):
    # Find the point with the lowest y-coordinate
    lowest = min(points, key=lambda p: (p[1], p[0]))

    # Sort the remaining points by angle
    def angle(p):
        dx = p[0] - lowest[0]
        dy = p[1] - lowest[1]
        # return math.atan2(dy, dx)
        if np.sqrt(dx*dx+dy*dy) ==0:
            ang = -1
        else:
            cos_val = dx/np.sqrt(dx*dx+dy*dy)
            ang = math.acos(cos_val)
            
        return ang
    
    sorted_points = sorted(points, key=angle)

    # Build the convex hull
    hull = [lowest, sorted_points[0]]
    for p in sorted_points[1:]:
        while len(hull) > 1 and orientation(hull[-2], hull[-1], p) <= 0:
            hull.pop()
        hull.append(p)

    return hull

def orientation(p, q, r):
    dx1 = q[0] - p[0]
    dy1 = q[1] - p[1]
    dx2 = r[0] - q[0]
    dy2 = r[1] - q[1]
    return dx1 * dy2 - dx2 * dy1

# def orientation(p, q, r):
#     val = ((q[1] - p[1]) * (r[0] - q[0]) -
#            (q[0] - p[0]) * (r[1] - q[1]))
#     if val == 0:
#         return 0  # collinear
#     elif val > 0:
#         return 1  # clock wise
#     else:
#         return 2  # counterclock wise

with open('vertices.dat', 'r') as vertices_file:
    
    vertices_dict=dict()
    
    # loop over the lines in the file
    for line in vertices_file:
        # do something with the line
        if ":" in line:
            index_str = line.split(":")[0]
            index = int(index_str)
            vertices_dict[index]=[]
            vertex_c=0
            continue
        
        if ("*" not in line):
            # vertices_dict[index].append((float(line.split("\t")[1]), float(line.split("\t")[0])))
            vertices_dict[index].append((float(line.split("\t")[0]), float(line.split("\t")[1])))
        else:
            continue


# # create a convex polygon
# vertices = [(0, 0), (1, 1), (2, 0), (1, -1), (1,0)]
# poly = Polygon(vertices, facecolor='blue', edgecolor='black')

# # create a plot and add the polygon to it
# fig, ax = plt.subplots()
# ax.add_patch(poly)

# # set the limits of the plot
# ax.set_xlim(-1, 3)
# ax.set_ylim(-2, 2)

# # show the plot
# plt.show()

polygons = []
values = []
for i in range(np.shape(com_file)[0]):
    polygons.append(np.array(vertices_dict[i]))
    values.append(i+1)
# define some example polygons and values
# polygons = [
#     [(0, 0), (1, 0), (1, 1), (0, 1)],  # square
#     [(1, 1), (2, 1), (2, 2), (1, 2)],  # square
#     [(1, 0), (1, 2), (0, 1)]]
# values = [1.7, 1.3, 0.2]

# create patches for each polygon and assign colors based on values
# patches = []
# colors = []
# for i, poly in enumerate(polygons):
#     patch = Polygon(poly, True)
#     patches.append(patch)
#     colors.append(values[i])

# # create a PatchCollection and set the colors
# pc = PatchCollection(patches, cmap='rainbow')
# pc.set_array(np.array(colors))

# # plot the polygons
# fig, ax = plt.subplots()
# ax.add_collection(pc)
# ax.autoscale()
# plt.show()
cmap = plt.get_cmap('rainbow')
vmin = 0
vmax = NumCells
# vmax = 16
norm = plt.Normalize(vmin, vmax)


# initConfig = np.loadtxt('initConfig.csv', delimiter=',')

copy_direction_list = [[-1,-1], [-1, 0], [-1,+1], [0, +1], [+1, +1], [+1, 0], [+1, -1], [0,-1] ]

# exit()
# fig, ax = plt.subplots()
for i in range(NumCells):
    color = cmap(norm(values[i]))
    # color = cmap(norm(initConfig[i]))
    # color = cmap(norm(1))
    main_points = polygons[i]
    points = main_points.copy()
    hull = ConvexHull(points)    
    convex_hull = graham_scan(points)
    convex_hull_plot = convex_hull.copy()
    for point_c in range(len(convex_hull_plot)):
        temp = convex_hull_plot[point_c][0].copy()
        convex_hull_plot[point_c][0] = convex_hull_plot[point_c][1].copy()
        convex_hull_plot[point_c][1] = temp.copy()
    plt.fill(*zip(*convex_hull_plot), color=color, alpha=1, edgecolor='black', linewidth=0.2)
    
    

    for copy_counter in range(0,8):
        
        copy_dir = copy_direction_list[copy_counter]
        
        points[:,0] = main_points[:,0].copy() + Lx*copy_dir[1]
        points[:,1] = main_points[:,1].copy() + Ly*copy_dir[0]

        margin_length = 5*np.sqrt(Lx*Ly/(NumCells))
        condition = (np.abs(points[0,0]-0.5*Lx)<0.5*Lx+margin_length) and \
                    (np.abs(points[0,1]-0.5*Ly)<0.5*Ly+margin_length)
        
        if condition: 
            hull = ConvexHull(points)
            convex_hull = graham_scan(points)
            convex_hull_plot = convex_hull.copy()
            for point_c in range(len(convex_hull_plot)):
                temp = convex_hull_plot[point_c][0].copy()
                convex_hull_plot[point_c][0] = convex_hull_plot[point_c][1].copy()
                convex_hull_plot[point_c][1] = temp.copy()
            plt.fill(*zip(*convex_hull_plot), color=color, alpha=1, edgecolor='black', linewidth=0.2)
        

    # plot the convex hull
    
    # ax.plot(points[:, 0], points[:, 1], 'o')
    # for simplex in hull.simplices:
    #     # ax.plot(points[simplex, 0], points[simplex, 1], 'k-', linewidth=0.5)
    #     # ax.fill(points[hull.vertices,0], points[hull.vertices,1], plt.cm.plasma(2*i+1), alpha=0.2)
        
    #     # poly = plt.Polygon(points[simplex], facecolor=cmap(values[i]), edgecolor='k')
    #     # ax.add_patch(poly)
    #     color = 'green'
    #     poly = plt.Polygon(points[simplex], alpha=0.5, edgecolor='k',linewidth=0.5)
    #     poly.set_fill('g')
    #     ax.add_patch(poly)
    #     # poly.set_capstyle('round')
    #     # plt.gca().add_patch(poly)
    
# plt.axis('equal')
# plt.xlim(-L, 2*L)
# plt.ylim(-L, 2*L)
plt.xlim(0, Ly)
plt.ylim(0, Lx)
plt.gca().set_aspect('equal')
plt.xlabel('Y')
plt.ylabel('X')
plt.grid()
plt.gca().invert_yaxis()

plt.plot([0,0],[0,Lx],linestyle='dashed', linewidth=0.6, color='black')
plt.plot([0,Ly],[0,0],linestyle='dashed', linewidth=0.6, color='black')
plt.plot([0,Ly],[Lx,Lx],linestyle='dashed', linewidth=0.6, color='black')
plt.plot([Ly,Ly],[0,Lx],linestyle='dashed', linewidth=0.6, color='black')
plt.savefig("vor_data.png", dpi=1000)
# plt.show()
