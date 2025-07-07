#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan  2 14:03:57 2024

@author: hossein
"""



import numpy as np
import math
from scipy.spatial import ConvexHull


class I_class:
    Ixx = 0.0
    Iyy = 0.0
    Ixy = 0.0
    center_x = 0.0
    center_y = 0.0
    
def inertia_calc_tri(points, center):
    
    
    
    I_info = I_class()
    
    I_info.center[0] = center[0]
    I_info.center[1] = center[1]
    
    centroid_x = (points[0][0]+points[1][0]+points[2][0])/3.0
    centroid_y = (points[0][1]+points[1][1]+points[2][1])/3.0
    
    u = [points[1][0]-points[0][0], points[1][1]-points[0][1]]
    v = [points[2][0]-points[0][0], points[2][1]-points[0][1]]
    
    area = 0.5*np.abs(u[0]*v[1] - u[1]*v[0])
    
    Ibar_xx = 0.
    Ibar_yy = 0.
    Ibar_xy = 0.
    
    I_info.Ixx = Ibar_xx + area * (centroid_x - I_info.center[0])**2
    I_info.Iyy = Ibar_yy + area * (centroid_y - I_info.center[1])**2
    I_info.Ixy = Ibar_xy + area * (centroid_x - I_info.center[0]) * (centroid_y - I_info.center[1])
    
    return I_info

def orientation(p, q, r):
    dx1 = q[0] - p[0]
    dy1 = q[1] - p[1]
    dx2 = r[0] - q[0]
    dy2 = r[1] - q[1]
    return dx1 * dy2 - dx2 * dy1

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

def I_calc_polygon(points, centroid):
    
    I_info = I_class()
    
    if np.abs(points[0][0]-centroid[0])<3.0:
        I_info.center_x = centroid[0].copy()
    elif points[0][0]-centroid[0] > Lx/2:
        I_info.center_x = (centroid[0]+Lx).copy()
    elif points[0][0]-centroid[0] < -Lx/2:
        I_info.center_x = (centroid[0]-Lx).copy()
    else:
        fff
    
    if np.abs(points[0][1]-centroid[1])<3.0:
        I_info.center_y = centroid[1].copy()
    elif points[0][1]-centroid[1] > Ly/2:
        I_info.center_y = (centroid[1]+Ly).copy()
    elif points[0][1]-centroid[1] < -Ly/2:
        I_info.center_y = (centroid[1]-Ly).copy()
    else:
        fff
        
    
    
    I_info.Ixx = 0.0
    I_info.Iyy = 0.0
    I_info.Ixy = 0.0
    
    N = len(points)
    
    for i in range(N):
        
        x_i = points[i][0] - I_info.center_x
        y_i = points[i][1] - I_info.center_y
        x_i_plus = points[(i+1)%N][0] - I_info.center_x
        y_i_plus = points[(i+1)%N][1] - I_info.center_y
        
        I_info.Ixx += (1/12)*(x_i*y_i_plus-x_i_plus*y_i)*(x_i**2 + x_i*x_i_plus + x_i_plus**2)
        I_info.Iyy += (1/12)*(x_i*y_i_plus-x_i_plus*y_i)*(y_i**2 + y_i*y_i_plus + y_i_plus**2)
        I_info.Ixy += (1/24)*(x_i*y_i_plus-x_i_plus*y_i)*(x_i*y_i_plus +2*x_i*y_i +2*x_i_plus*y_i_plus + x_i_plus*y_i)

    return I_info
    
# points = []
# points.append((2,0))
# points.append((1,1))
# points.append((-1,1))
# points.append((-1,-1))
# points.append((0.5,-1))

# convex_hull = graham_scan(points)

siteComXY = np.loadtxt('lattice/siteComXY.csv', delimiter=',')
LxLy =  np.loadtxt('lattice/LxLy.dat', delimiter=',')
Lx = LxLy[0]
Ly = LxLy[1]

with open('lattice/vertices.dat', 'r') as vertices_file:
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
            vertices_dict[index].append((float(line.split("\t")[0]), float(line.split("\t")[1])))
        else:
            continue
    
    
NSites = len(vertices_dict)
polygons = []
flag =0
I_list = []

I_info_array = np.zeros([NSites, 5])

for siteC in range(NSites):
    convex_hull = graham_scan(vertices_dict[siteC])
    if len(convex_hull) != len(vertices_dict[siteC]):
        flag = 1
    polygons.append(convex_hull)
    
    I_info_site = I_class()
    I_info_site = I_calc_polygon(convex_hull, (siteComXY[siteC, 0], siteComXY[siteC, 1]))
    I_list.append(I_info_site)
    
    I_info_array[siteC,:] = [I_info_site.center_x, I_info_site.center_y, I_info_site.Ixx, I_info_site.Iyy, I_info_site.Ixy]

np.savetxt('I_info_array.csv', X=I_info_array, fmt='%1.8f', delimiter=',')

