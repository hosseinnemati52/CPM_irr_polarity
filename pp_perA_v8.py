#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 20 08:45:30 2024

@author: hossein
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec  6 11:57:23 2023

@author: hossein
"""



import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as ps
import matplotlib.cm as cm
import math
from scipy.spatial import ConvexHull
import os

class distributionClass:
    N_tot = 0
    count = None
    bins = None
    dx = 0.0
    binCenters = None
    

class polygonInfoClass:
    
    vertices = []
    vertices.clear()
    vertices = []
    
    ordered_vertices = []
    ordered_vertices.clear()
    ordered_vertices = []
    
    center_x = 0.0
    center_y = 0.0
    
    area = 0.0
    
    Ixx = 0.0
    Iyy = 0.0
    Ixy = 0.0
    
    AR_Vor = 0.0
    circVor = 0.0

class fitClass:
  x = None
  y = None
  loglogFit = None
  linFit = None

class eqSamplingTimesClass:
    snapshotCList = []
    timesList = []
    sampleCList = []
    bunchCList = []
    

def PDF_creator(distObject, write_file_name, quantity_name):
    
    xlim_upper_c = len(distObject.bins) -1
    xlim_lower_c = 0
    xlim_upper = distObject.bins[xlim_upper_c]
    xlim_lower = distObject.bins[xlim_lower_c]
    
    margin = 3
    while np.sum(distObject.count[xlim_upper_c-margin:xlim_upper_c+1])==0:
        xlim_upper_c -= 1
        xlim_upper = distObject.bins[xlim_upper_c]
        
    while np.sum(distObject.count[xlim_lower_c:xlim_lower_c+margin+1])==0:
        xlim_lower_c += 1
        xlim_lower = distObject.bins[xlim_lower_c]
    
    xlim=[xlim_lower-2*distObject.dx, xlim_upper+2*distObject.dx]
    
    
    marker_sym_dict = dict()
    marker_sym_dict[3] = "^"
    marker_sym_dict[4] = "s"
    marker_sym_dict[5] = "p"
    marker_sym_dict[6] = "h"
    marker_sym_dict[7] = "$7$"
    marker_sym_dict[8] = "$8$"
    marker_sym_dict[9] = "$9$"
    
    polygons_flag = 0
    
    # if quantity_name == 'AR_Vor':
        
    #     polygons_flag =1
        
    #     polygons_val = dict()
    #     polygons_val[3] = 0.7776
    #     polygons_val[4] = 0.8862
    #     polygons_val[5] = 0.9299
    #     polygons_val[6] = 0.95231
    #     polygons_val[7] = 0.9654
    #     polygons_val[8] = 0.97368
    #     polygons_val[9] = 0.97931
        
        
    
    if quantity_name == 'circVor':
        polygons_flag =1
        polygons_val = polygons_circ.copy()
    
    if quantity_name == 'qVor':
        
        polygons_flag =1
        
        polygons_val = dict()
        polygons_val[3] = 0.7776
        polygons_val[4] = 0.8862
        polygons_val[5] = 0.9299
        polygons_val[6] = 0.95231
        polygons_val[7] = 0.9654
        polygons_val[8] = 0.97368
        polygons_val[9] = 0.97931
    
    if quantity_name == 'qInvVor':
        polygons_flag =1
        
        polygons_val = dict()
        polygons_val[3] = 2*np.sqrt(np.pi)/0.7776
        polygons_val[4] = 2*np.sqrt(np.pi)/0.8862
        polygons_val[5] = 2*np.sqrt(np.pi)/0.9299
        polygons_val[6] = 2*np.sqrt(np.pi)/0.95231
        polygons_val[7] = 2*np.sqrt(np.pi)/0.9654
        polygons_val[8] = 2*np.sqrt(np.pi)/0.97368
        polygons_val[9] = 2*np.sqrt(np.pi)/0.97931
        
    # if quantity_name == 'AR':
    #     xlim = [1, 4]
        
    if quantity_name == 'q':
        
        polygons_flag =1
        
        polygons_val = dict()
        polygons_val[3] = 0.7776
        polygons_val[4] = 0.8862
        polygons_val[5] = 0.9299
        polygons_val[6] = 0.95231
        polygons_val[7] = 0.9654
        polygons_val[8] = 0.97368
        polygons_val[9] = 0.97931
        
    if quantity_name == 'q_inv':
        
        polygons_flag =1
        
        polygons_val = dict()
        polygons_val[3] = 2*np.sqrt(np.pi)/0.7776
        polygons_val[4] = 2*np.sqrt(np.pi)/0.8862
        polygons_val[5] = 2*np.sqrt(np.pi)/0.9299
        polygons_val[6] = 2*np.sqrt(np.pi)/0.95231
        polygons_val[7] = 2*np.sqrt(np.pi)/0.9654
        polygons_val[8] = 2*np.sqrt(np.pi)/0.97368
        polygons_val[9] = 2*np.sqrt(np.pi)/0.97931
    
    if quantity_name == 'circ':
        polygons_flag =1
        polygons_val = polygons_circ.copy()
        
    
    # if quantity_name == 'hex':
    #     xlim = [0, 1]

    # old writing was not working correctly
    PDF_mean = 1.0*distObject.count/(distObject.N_tot*distObject.dx)
    PDF_err  = 1.0*np.sqrt(distObject.count)/(distObject.N_tot*distObject.dx)
    # f = open(write_file_name, "w")
    # f.write("mean = "+str(PDF_mean)+"\n")
    # f.write("----------------------------------------------\n")
    # f.write("std  = "+str(PDF_err)+"\n")
    # f.write("----------------------------------------------\n")
    # f.write("bins = "+str(distObject.bins))
    # f.close()
    
    # new writing
    PDF_data = np.zeros([3,len(distObject.binCenters)])
    PDF_data[0,:] = PDF_mean
    PDF_data[1,:] = PDF_err
    PDF_data[2,:] = distObject.binCenters
    np.savetxt(write_file_name, X = PDF_data, delimiter=',', fmt='%1.6f')

    x_plot = distObject.binCenters
    y_plot = PDF_mean
    menStd     = PDF_err
    width      = distObject.dx
    
    mean_val = np.sum(distObject.count * distObject.binCenters)/distObject.N_tot
    
    for yscale in ["linear", "log"]:
        plt.figure()
        plt.bar(x_plot, y_plot, width=width, zorder=1)
        plt.errorbar(x_plot, y_plot, yerr=menStd, fmt='none', capsize=0.2, ecolor='k', elinewidth=1.0, zorder=2)
        # plt.hist(AR_data, bins = bin_list)
        plt.xlim(xlim)
        plt.xlabel(quantity_name)
        plt.ylabel("PDF("+quantity_name+")")
        plt.axvline(x=mean_val, linestyle='dashed', color='k', linewidth=1, zorder=3)
        
        if polygons_flag:
            for poly in range(3,10):
                plt.axvline(x=polygons_val[poly], linestyle='dashed', color='k', linewidth=0.5, zorder=5)
                marker_sym = marker_sym_dict[poly]
                plt.scatter(polygons_val[poly]+0.00, 0.8*np.max(y_plot), marker=marker_sym, s=70, zorder=4)
        
        if yscale == "log":
                plt.ylim([(5e-6),np.max(y_plot)*5])
                
        # plt.axvline(x=min_val, linestyle='dashed', color='k', linewidth=0.5)
        # plt.axvline(x=max_val, linestyle='dashed', color='k', linewidth=0.5)
        plt.axhline(y=0, linestyle='dashed', color='k', linewidth=1)
        plt.yscale(yscale)
        title = quantity_name+", Avg: "+str(round(mean_val,4))
        plt.title(title)
        plt.savefig(quantity_name+"_PDF_"+yscale+".PNG", dpi = 300)
        # plt.show()
    return

def polygonCircAR_calc():
    
    polygons_AR = dict()
    polygons_circ = dict()
    r = 1.0
    for n in range(3, 10):
        points = []
        del(points)
        points = np.zeros([n, 2])
        
        for i in range(n):
            theta = i*(2*np.pi)/n
            points[i,0] = r * np.cos(theta)
            points[i,1] = r * np.sin(theta)
        
        A = 0.0
        for i in range(n):
            x_i = points[i,0]
            y_i = points[i,1]
            x_i_plus = points[(i+1)%n,0]
            y_i_plus = points[(i+1)%n,1]
            A += 0.5*(x_i*y_i_plus-x_i_plus*y_i)
        
        Ixx = 0.0
        Iyy = 0.0
        Ixy = 0.0
        for i in range(n):
            x_i = points[i,0]
            y_i = points[i,1]
            x_i_plus = points[(i+1)%n,0]
            y_i_plus = points[(i+1)%n,1]
            
            Ixx += (1/12)*(x_i*y_i_plus-x_i_plus*y_i)*(x_i**2 + x_i*x_i_plus + x_i_plus**2)
            Iyy += (1/12)*(x_i*y_i_plus-x_i_plus*y_i)*(y_i**2 + y_i*y_i_plus + y_i_plus**2)
            Ixy += (1/24)*(x_i*y_i_plus-x_i_plus*y_i)*(x_i*y_i_plus + 2*x_i*y_i + 2* x_i_plus*y_i_plus + x_i_plus*y_i)
        
        I_min = 0.5*((Ixx+Iyy)-np.sqrt((Ixx-Iyy)*(Ixx-Iyy)+4.0*Ixy*Ixy))
        I_max = 0.5*((Ixx+Iyy)+np.sqrt((Ixx-Iyy)*(Ixx-Iyy)+4.0*Ixy*Ixy))
        
        polygons_AR[n] = np.sqrt(I_max/I_min)
        polygons_circ[n] = ((A**2)/(2*np.pi))/(Ixx+Iyy)
    
    return polygons_AR, polygons_circ

def reverse_xy(array):
    
    rev_array=array.copy()
    n_row = np.shape(array)[0]
    
    for i in range(n_row):
        rev_array[i,0] = array[i,1]
        rev_array[i,1] = array[i,0]
    return rev_array
    
def VorCellsPlotter(VorCellInfoList, colorKeyVector, savePlotAddress, saveSwitch, cmap, com_switch, plot_title, xCom_instant, yCom_instant):
    
    N_Vor = len(VorCellInfoList)-1
    # xCom_instant =np.zeros(len(VorCellInfoList))
    # yCom_instant =np.zeros(len(VorCellInfoList))
    # for vorC in range(1,N_Vor+1):
    #     xCom_instant[vorC] = VorCellInfoList[vorC].center_x
    #     yCom_instant[vorC] = VorCellInfoList[vorC].center_y
        
    plt.figure()
    
    sigmaMat = colorKeyVector.copy()
    # NSites = len(colorKeyVector)
    
    # with open(verticesFileAddress, 'r') as vertices_file:
    #     vertices_dict=dict()
    #     # loop over the lines in the file
    #     for line in vertices_file:
    #         # do something with the line
    #         if ":" in line:
    #             index_str = line.split(":")[0]
    #             index = int(index_str)
    #             vertices_dict[index]=[]
    #             vertex_c=0
    #             continue
            
    #         if ("*" not in line):
    #             vertices_dict[index].append((float(line.split("\t")[1]), float(line.split("\t")[0])))
    #         else:
    #             continue
    
    polygons = [None]
    values = [None]
    for i in range(1,N_Vor+1):
        polygons.append(reverse_xy(np.array(VorCellInfoList[i].ordered_vertices)))
        values.append(sigmaMat[i])
        # values.append(i+1)
    
    # cmap = plt.get_cmap('rainbow')
    # if np.min(sigmaMat)==1 or np.min(sigmaMat)==0:
    #     vmin = 0
    #     vmax = NumCells
    #     lineColor = 'k'
    # else:
    #     vmin = 0.0
    #     vmax = 1.0
    #     lineColor = 'w'
    vmin = np.min(sigmaMat[1:])
    vmax = np.max(sigmaMat[1:])
    lineColor = 'k'
    norm = plt.Normalize(vmin, vmax)
    
    copy_direction_list = [[-1,-1], [-1, 0], [-1,+1], [0, +1], [+1, +1], [+1, 0], [+1, -1], [0,-1] ]
    
    length_scale = np.sqrt(Lx*Ly/N_Vor)
    margin = 2*length_scale
    edge_plot_check = [] #a list of sets containing the indices of sites
    for i in range(1, N_Vor+1):
        color = cmap(norm(values[i]))
        main_points = polygons[i]
        points = main_points.copy()
        hull = ConvexHull(points)
        convex_hull = graham_scan(points)
        plt.fill(*zip(*convex_hull), facecolor=color, alpha=1, edgecolor='black', linewidth=0.1, zorder=1)
        
        # edge_plot_list = []
        # for neigh in neigh_dict[i]:
        #     if (sigmaMat[i] != sigmaMat[neigh]) and (set({neigh, i}) not in edge_plot_check):
        #         edgePlot = [[],[]]  #[[y1, y2], [x1, x2]]
        #         for vertice_i in vertices_dict[i]:
        #             for vertice_neigh in vertices_dict[neigh]:
        #                 error = ((vertice_i[0]-vertice_neigh[0])**2+(vertice_i[1]-vertice_neigh[1])**2)**0.5 # r = pol(dx, dy)
        #                 if error< (1e-6)*length_scale:
        #                     edgePlot[0].append(vertice_i[0]) #this is y in the simulation
        #                     edgePlot[1].append(vertice_i[1]) #this is x in the simulation
                
        #         edge_plot_list.append(edgePlot)
        #     else:
        #         pass
        
        # for edgePlot in edge_plot_list:
        #     plt.plot(edgePlot[0], edgePlot[1], linewidth=0.1, color=lineColor, zorder=1)
        
        
        for copy_counter in range(0,8):
            
            copy_dir = copy_direction_list[copy_counter]
            
            points[:,0] = main_points[:,0].copy() + Ly*copy_dir[1]
            points[:,1] = main_points[:,1].copy() + Lx*copy_dir[0]
            
            condition = (abs(points[0,0]-Ly/2)<Ly/2+margin) and (abs(points[0,1]-Lx/2)<Lx/2+margin) 
            
            if condition:
                hull = ConvexHull(points)
                convex_hull = graham_scan(points)
                plt.fill(*zip(*convex_hull), facecolor=color, alpha=1, edgecolor='black', linewidth=0.1, zorder=1)
                
                # for edgePlot in edge_plot_list:
                #     plt.plot(np.array(edgePlot[0])+ Ly*copy_dir[1], np.array(edgePlot[1])+ Lx*copy_dir[0], linewidth=0.1, color=lineColor, zorder=1)
        
        # for j in range():
    
        
        
    # plt.figure()
    # cax = cm.ScalarMappable(cmap)
    # cax.set_array(sigmaMat)
    # plt.colorbar()
    plt.colorbar(cm.ScalarMappable(norm=norm, cmap=cmap))
    plt.xlim(0-yShift, Ly-yShift)
    plt.ylim(0-xShift, Lx-xShift)
    plt.gca().set_aspect('equal')
    plt.xlabel('Y')
    plt.ylabel('X')
    if lattice_switch =='irr':
        plt.xticks(np.arange(0, 175+1, 25.0))
    elif lattice_switch =='sq':
        plt.xticks(np.arange(0, 175+1, 25.0))
    plt.grid()
    plt.gca().invert_yaxis()
    plt.title(plot_title+", Avg: "+str(round(np.mean(sigmaMat[1:]),4)))
    
    plt.plot([0 -yShift,0 -yShift],[0 -xShift,Lx-xShift],linestyle='dashed', linewidth=0.6, color='black')
    plt.plot([0 -yShift,Ly-yShift],[0 -xShift,0 -xShift],linestyle='dashed', linewidth=0.6, color='black')
    plt.plot([0 -yShift,Ly-yShift],[Lx-xShift,Lx-xShift],linestyle='dashed', linewidth=0.6, color='black')
    plt.plot([Ly-yShift,Ly-yShift],[0 -xShift,Lx-xShift],linestyle='dashed', linewidth=0.6, color='black')
    
    
    if com_switch:
        plt.scatter(yCom_instant[1:]%Ly, xCom_instant[1:]%Lx, s=0.1, c=lineColor, zorder=2)
    
    plt.savefig(savePlotAddress, dpi=800)
    # plt.show()
    
    return

def verticesReader(fileAddress):
    vertices_dict=dict()
    vertices_dict.clear()
    vertices_dict=dict()
    
    with open(fileAddress, 'r') as vertices_file:
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
                vertices_dict[index].append((float(line.split("\t")[0]) -xShift , float(line.split("\t")[1])-yShift) )
            else:
                continue
        
    return vertices_dict

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

def VorCellInfoCalc(points):
    
    VorCellInfo = polygonInfoClass()
    
    VorCellInfo.vertices = points.copy()
    
    ################ ordering vertices ###############
    convex_hull = graham_scan(VorCellInfo.vertices.copy())
    VorCellInfo.ordered_vertices = convex_hull.copy()
    ################ ordering vertices ###############
    
    ################ area clac ###############
    VorCellInfo.area = 0.0
    N = len(VorCellInfo.ordered_vertices)
    for i in range(N):
        x_i = VorCellInfo.ordered_vertices[i][0]
        y_i = VorCellInfo.ordered_vertices[i][1]
        x_i_plus = VorCellInfo.ordered_vertices[(i+1)%N][0]
        y_i_plus = VorCellInfo.ordered_vertices[(i+1)%N][1]
        
        VorCellInfo.area += 0.5*(x_i*y_i_plus - x_i_plus*y_i)
    ################ area clac ###############
    
    ################ centroid clac ###############
    VorCellInfo.center_x = 0.0
    VorCellInfo.center_y = 0.0
    for i in range(N):
        x_i = VorCellInfo.ordered_vertices[i][0]
        y_i = VorCellInfo.ordered_vertices[i][1]
        x_i_plus = VorCellInfo.ordered_vertices[(i+1)%N][0]
        y_i_plus = VorCellInfo.ordered_vertices[(i+1)%N][1]
        
        VorCellInfo.center_x += (x_i + x_i_plus)*(x_i*y_i_plus - x_i_plus*y_i)
        VorCellInfo.center_y += (y_i + y_i_plus)*(x_i*y_i_plus - x_i_plus*y_i)
    
    VorCellInfo.center_x /= (6.0*VorCellInfo.area)
    VorCellInfo.center_y /= (6.0*VorCellInfo.area)
    ################ centroid clac ###############
    
    ################ inertia clac ###############
    VorCellInfo.Ixx = 0.0
    VorCellInfo.Iyy = 0.0
    VorCellInfo.Ixy = 0.0
    for i in range(N):
        x_i = VorCellInfo.ordered_vertices[i][0] - VorCellInfo.center_x
        y_i = VorCellInfo.ordered_vertices[i][1] - VorCellInfo.center_y
        x_i_plus = VorCellInfo.ordered_vertices[(i+1)%N][0] - VorCellInfo.center_x
        y_i_plus = VorCellInfo.ordered_vertices[(i+1)%N][1] - VorCellInfo.center_y
        
        VorCellInfo.Ixx += (x_i*y_i_plus - x_i_plus*y_i) * (x_i**2 + x_i*x_i_plus +x_i_plus**2)
        VorCellInfo.Iyy += (x_i*y_i_plus - x_i_plus*y_i) * (y_i**2 + y_i*y_i_plus +y_i_plus**2)
        VorCellInfo.Ixy += (x_i*y_i_plus - x_i_plus*y_i) * (x_i*y_i_plus + 2*x_i*y_i + 2*x_i_plus*y_i_plus +x_i_plus*y_i)
    
    VorCellInfo.Ixx /= 12.0
    VorCellInfo.Iyy /= 12.0
    VorCellInfo.Ixy /= 24.0
    ################ inertia clac ###############
    
    ################ AR clac ###############
    Ixx = VorCellInfo.Ixx
    Iyy = VorCellInfo.Iyy
    Ixy = VorCellInfo.Ixy
    I_min = 0.5*((Ixx+Iyy)-np.sqrt((Ixx-Iyy)*(Ixx-Iyy)+4.0*Ixy*Ixy));
    I_max = 0.5*((Ixx+Iyy)+np.sqrt((Ixx-Iyy)*(Ixx-Iyy)+4.0*Ixy*Ixy));
    VorCellInfo.AR_Vor = np.sqrt(I_max/I_min)
    ################ AR clac ###############
    
    ################ circ clac ###############
    VorCellInfo.circVor = (1/(2*np.pi))*(VorCellInfo.area**2)/(VorCellInfo.Ixx+VorCellInfo.Iyy)
    ################ circ clac ###############
    
    return VorCellInfo

def distributionFuncExtended(quantity_name):
    
    ppDataFolderName = pp_data_address
    ppPlotFolderName = pp_plot_address
    ppVorDataFolderName = pp_vor_address
    
    xlabel_description = ""
    if quantity_name == 'edgesVor':
        dx = .1
        min_val = 0.
        #max_val = 10.
        max_val = 1.5*np.sqrt(AvgCellArea)
        total_data = []
        f = open(pp_vor_address+"/eq_edges_data.txt", "r")
        for i in range(len(eq_indices)):
            edge_vor_data_str = f.readline()
            edge_vor_data = list(map(float, edge_vor_data_str.split('\t')))
            total_data = total_data + edge_vor_data.copy()
        f.close()
        write_file_name = ppVorDataFolderName+"/edgesVor_freq.csv"
        this_plot_switch = plot_switches_dict['q']
    
    if quantity_name == 'edgesVorLog':
        dx = .015 * 2
        min_val = -5.0
        #max_val = 1.0
        max_val = np.log10(1.5*np.sqrt(AvgCellArea))
        total_data = []
        f = open(pp_vor_address+"/eq_edges_data.txt", "r")
        for i in range(len(eq_indices)):
            edge_vor_data_str = f.readline()
            edge_vor_data = list(np.log10(np.array(list(map(float, edge_vor_data_str.split('\t'))))))
            total_data = total_data + edge_vor_data.copy()
        f.close()
        write_file_name = ppVorDataFolderName+"/edgesVorLog_freq.csv"
        this_plot_switch = plot_switches_dict['q']
        
    if quantity_name == 'distancesVor':
        dx = .1
        min_val = 0.
        #max_val = 12.
        max_val = 2.0*np.sqrt(AvgCellArea)
        total_data = []
        f = open(pp_vor_address+"/eq_distances_data.txt", "r")
        for i in range(len(eq_indices)):
            edge_vor_data_str = f.readline()
            edge_vor_data = list(map(float, edge_vor_data_str.split('\t')))
            total_data = total_data + edge_vor_data.copy()
        f.close()
        write_file_name = ppVorDataFolderName+"/distancesVor_freq.csv"
        this_plot_switch = plot_switches_dict['q']
        
    if quantity_name == 'AR_Vor':
        dist_eq = np.loadtxt(ppVorDataFolderName+"/AR_VorCellsEq.csv", delimiter=',')[1:, :]
        dx = 0.005 * 2
        min_val = 1.0
        # max_val = 3.6
        max_val = 4.0
        write_file_name = ppVorDataFolderName+"/AR_Vor_freq.csv"
        this_plot_switch = plot_switches_dict['AR']
    
    if quantity_name == 'circVor':
        dist_eq = np.loadtxt(ppVorDataFolderName+"/circVorCellsEq.csv", delimiter=',')[1:, :]
        dx = 0.001 *2
        min_val = 0.0
        max_val = 1.0
        write_file_name = ppVorDataFolderName+"/circVor_freq.csv"
        this_plot_switch = plot_switches_dict['circ']
    
    if quantity_name == 'qVor':
        dist_eq = np.loadtxt(ppVorDataFolderName+"/qCellEq.csv", delimiter=',')[1:, :]
        dx = 0.0005 *2
        min_val = 0.0
        # max_val = 3.6
        max_val = 1.0
        write_file_name = ppVorDataFolderName+"/qVor_freq.csv"
        this_plot_switch = plot_switches_dict['q']
    
    if quantity_name == 'qInvVor':
        dist_eq = np.loadtxt(ppVorDataFolderName+"/qCellEq.csv", delimiter=',')[1:, :]
        dist_eq = 2*(np.pi**0.5)/dist_eq
        dx = 0.001 *2
        min_val = 2*(np.pi**0.5)
        max_val = 10.0
        write_file_name = ppVorDataFolderName+"/qInvVor_freq.csv"
        this_plot_switch = plot_switches_dict['q']
        xlabel_description = '('+r'$P_{vor}/\sqrt{A_{vor}}$'+')'
        
    if quantity_name == 'AR':
        dist_eq = np.loadtxt(ppDataFolderName+"/AR_dist_eq.csv", delimiter=',')[1:, :]
        dx = 0.005 *2
        min_val = 1.0
        # max_val = 3.6
        max_val = 4.0
        write_file_name = ppDataFolderName+"/AR_freq.csv"
        this_plot_switch = plot_switches_dict[quantity_name]
        
    if quantity_name == 'q':
        dist_eq = np.loadtxt(ppDataFolderName+"/isoperi_dist_eq.csv", delimiter=',')[1:, :]
        dx = 0.002*2
        min_val = 0.0
        # max_val = 3.6
        max_val = 1.0
        write_file_name = ppDataFolderName+"/isoperi_freq.csv"
        this_plot_switch = plot_switches_dict[quantity_name]
        
    if quantity_name == 'q_inv':
        dist_eq = np.loadtxt(ppDataFolderName+"/isoperi_dist_eq.csv", delimiter=',')[1:, :]
        dist_eq = 2*(np.pi**0.5)/dist_eq
        dx = 0.01 *2
        min_val = 2*(np.pi**0.5)
        max_val = 10.0
        write_file_name = ppDataFolderName+"/isoperi_inv_freq.csv"
        this_plot_switch = plot_switches_dict['q']
        xlabel_description = '('+r'$P/\sqrt{A}$'+')'
    
    if quantity_name == 'circ':
        dist_eq = np.loadtxt(ppDataFolderName+"/circ_dist_eq.csv", delimiter=',')[1:, :]
        dx = 0.001 * 2
        min_val = 0.0
        max_val = 1.0
        write_file_name = ppDataFolderName+"/circ_freq.csv"
        this_plot_switch = plot_switches_dict['circ']
    
    if quantity_name == 'hex':
        dist_eq = np.loadtxt(ppDataFolderName+"/hex_dist_eq.csv", delimiter=',')[1:, :]
        dx = 0.005 * 2
        min_val = 0.0
        max_val = 1.0
        write_file_name = ppDataFolderName+"/hex_freq.csv"
        this_plot_switch = plot_switches_dict['order']
        
    
    try:
        len(total_data)
        title = r"$\alpha = $"+str(Alpha)+"\n Avg: "+str(round(np.mean(total_data), 3))+", N_tot = "+str(len(total_data))
    except:
        total_data = []
        for eqSampleC in range(np.shape(dist_eq)[1]):
            total_data = total_data + list(dist_eq[1:,eqSampleC])
        title = r"$\alpha = $"+str(Alpha)+"\n Avg: "+str(round(np.mean(total_data), 3))+", N_tot = "+str(np.shape(dist_eq)[1])+"*"+str(NumCells)+" = "+str(np.shape(dist_eq)[1]*NumCells)
    
    # q3, q1 = np.percentile(total_data, [75 ,25])
    # iqr = q3 - q1
    # dx_FD = 2*iqr/(len(total_data)**0.333)
    
    # if dx <0.75*dx_FD:
    #     print("bin size error!")
    #     print(quantity_name)
    #     ddddffff
    
    bin_list = []
    bin_list.append(min_val-dx/2)
    tracker = min_val + dx/2
    while (tracker <= max_val):
        bin_list.append(round(tracker, 4))
        tracker = tracker + dx
    bin_list.append(tracker)
    bin_list.append(tracker+dx)
    
    nBins = len(bin_list)-1
    
    
        
    freq_hist = np.histogram(total_data, bins = bin_list)

    freq_hist_err = (1.0*freq_hist[0])**0.5

    # f = open(write_file_name, "w")
    # f.write("mean = "+str(1.0*freq_hist[0]/(len(total_data)))+"\n")
    # f.write("----------------------------------------------\n")
    # f.write("std  = "+str(1.0*freq_hist_err/(len(total_data)))+"\n")
    # f.write("----------------------------------------------\n")
    # f.write("bins = "+str(bin_list))
    # f.close()
    
    freq_file = np.zeros([2,len(freq_hist[0])])
    freq_file[0,:] = freq_hist[0].copy()
    freq_file[1,:] = (np.array(bin_list[1:])-0.5*dx).copy()
    np.savetxt(write_file_name, X=freq_file, delimiter=',', fmt='%.6f')


    x_plot = np.array(bin_list[1:])-0.5*dx
    if this_plot_switch:
        y_plot = freq_hist[0]/len(total_data)
        menStd     = freq_hist_err/len(total_data)
        width      = dx
        
        for yscale in ["linear", "log"]:
            plt.figure()
            plt.bar(x_plot, y_plot, width=width)
            plt.errorbar(x_plot, y_plot, yerr=menStd, fmt='none', capsize=0.2, ecolor='k', elinewidth=1.0)
            # plt.hist(AR_data, bins = bin_list)
            plt.xlim([np.min(total_data)-3*dx,np.max(total_data)+3*dx])
            if yscale == "log":
                plt.ylim([(5e-6),np.max(y_plot)*5])
            plt.xlabel(quantity_name+xlabel_description)
            plt.ylabel("freq("+quantity_name+")")
            plt.axvline(x=np.mean(total_data), linestyle='dashed', color='k', linewidth=1)
            plt.axvline(x=min_val, linestyle='dashed', color='k', linewidth=0.5)
            plt.axvline(x=max_val, linestyle='dashed', color='k', linewidth=0.5)
            plt.axhline(y=0, linestyle='dashed', color='k', linewidth=1)
            plt.yscale(yscale)
            plt.title(title)
            plt.savefig(ppPlotFolderName+"/"+quantity_name+"_freq_eq_"+yscale+".PNG", dpi = 300)
            # plt.show()
    
    distObject = distributionClass()
    del(distObject)
    distObject = distributionClass()
    
    distObject.N_tot = len(total_data)
    distObject.count = freq_hist[0].copy()
    distObject.bins = np.array(bin_list).copy()
    distObject.dx = dx
    distObject.binCenters = x_plot.copy()
    
    plt.close()
    
    return distObject

def fancy_qVor_PDF_plotter():
    
    polygons_q = dict()
    polygons_q[3] = 0.7776
    polygons_q[4] = 0.8862
    polygons_q[5] = 0.9299
    polygons_q[6] = 0.95231
    polygons_q[7] = 0.9654
    polygons_q[8] = 0.97368
    polygons_q[9] = 0.97931
    
    marker_sym_dict = dict()
    marker_sym_dict[3] = "^"
    marker_sym_dict[4] = "s"
    marker_sym_dict[5] = "p"
    marker_sym_dict[6] = "h"
    marker_sym_dict[7] = "$7$"
    marker_sym_dict[8] = "$8$"
    marker_sym_dict[9] = "$9$"
    
    ################# q_Vor dist plot ####################
    plt.figure()
    x_plot = qVor_distObject.binCenters
    y_plot = 1.0*qVor_distObject.count/(qVor_distObject.N_tot*qVor_distObject.dx)
    y_err =  np.sqrt(1.0*qVor_distObject.count)/(qVor_distObject.N_tot*qVor_distObject.dx)
    plt.plot(x_plot, y_plot)
    plt.fill_between(x_plot, y_plot-y_err, y_plot+y_err, alpha=0.5)
    plt.xlim([0.85,1])
    plt.xlabel("qVor")
    plt.ylabel("PDF(qVor)")
    for poly in range(3,10):
        plt.axvline(x=polygons_q[poly], linestyle='dashed', color='k', linewidth=0.5)
        # plt.figure()
        # ps.RegularPolygon((polygons_q[poly], 0.8*np.max(y_plot)), poly, radius=0.05, orientation=0)
        marker_sym = marker_sym_dict[poly]
        plt.scatter(polygons_q[poly]+0.00, 0.8*np.max(y_plot), marker=marker_sym, s=70)
    # plt.grid()
    plt.axhline(y=0, linestyle='dashed', color='k', linewidth=1)
    plt.title(r"$\alpha = $"+str(Alpha))
    plt.savefig("q_vor_PDF_a"+str(int(np.floor(Alpha)))+"dot"+str(int(100*Alpha)%100)+".PNG", dpi = 300)
    ################# q_Vor dist plot ####################
	
    plt.close()
    
    return

def fittingFuncSmart(beta_threshold):
    
    try:
        fitting_begin_array = np.loadtxt('../fitting_begin_sq.csv', delimiter=',')
    except:
        try:
            fitting_begin_array = np.loadtxt('../fitting_begin_irr.csv', delimiter=',')
        except:
            fitting_begin_array = np.loadtxt('../fitting_begin_hex.csv', delimiter=',')
        
    
    fitting_begin_dict = dict()
    for alphaC in range(np.shape(fitting_begin_array)[1]):
        fitting_begin_dict[fitting_begin_array[0,alphaC]] = fitting_begin_array[1,alphaC]
        N_BLOCKS = int(fitting_begin_array[2,alphaC])
    
    jumping_runs = []
    
    # for n_blocks in [n_runs]: # This list MUST have just one element
    for n_blocks in [N_BLOCKS]: # This list MUST have just one element
    
        if n_runs%n_blocks ==0:
            block_size_list = n_blocks*[int(n_runs/n_blocks)]
        else:
            block_size_list = n_blocks*[0]
            counter = n_runs
            i=0
            while counter>0:
                block_size_list[i%n_blocks] += 1
                i += 1
                counter -= 1
            # block_size_typical = round(n_runs/n_blocks)
            # block_size_list = []
            # while counter >  block_size_typical+1:
            #     block_size_list.append(block_size_typical)
            #     counter = counter - block_size_typical
            # block_size_list.append(counter)
        
        exponents_list = []
        diff_list = []
        
        # if Alpha < 1.8 + (1e-10):
        #     t_begin_fit = t_w + 1
        # elif Alpha < 2.6 + (1e-10):
        #     t_begin_fit = t_w + 1
        # else:
        #     t_begin_fit = t_w + 1
            
        
        t_begin_fit = t_w + fitting_begin_dict[Alpha]
        
        t_begin_fit_ind = t_w_ind
        while time[t_begin_fit_ind] <= t_begin_fit:
            t_begin_fit_ind += 1
            
        
        all_fitting_slopes = []
        
        runC=1
        for blockC in range(len(block_size_list)):
            block_size = block_size_list[blockC]
            
            msd_block = np.zeros([block_size, totalSnapshotsNum])
            
            for runC_dummy in range(block_size):
                
                pp_data_address = "run_"+str(runC)+"/pp_data"
                msdAvgStd = np.loadtxt(pp_data_address +'/msdAvgStd.csv', delimiter=',')
                
                msd_block[runC_dummy,:] = msdAvgStd[0,:]
                
                runC = runC+1
            
            fitElement = fitClass()
            
            fitElement.x = time[t_begin_fit_ind:]-t_w
            fitElement.y = np.mean(msd_block, axis=0)[t_begin_fit_ind:]
            
            # fitElement.linFit = np.polyfit(fitElement.x, fitElement.y, 1)
            # fitElement.linSlope = (fitElement.linFit)[0]
            fitElement.loglogFit = np.polyfit(np.log10(fitElement.x[1:]), np.log10(fitElement.y[1:]), 1)
            fitElement.loglogSlope = (fitElement.loglogFit)[0]
            all_fitting_slopes.append(fitElement.loglogSlope)
            
            # plt.figure()
            # plt.scatter(np.log10(fitElement.x[1:]), np.log10(fitElement.y[1:]) )
            # plt.plot(np.log10(fitElement.x[1:]), np.log10(fitElement.x[1:])* (fitElement.loglogSlope) + ((fitElement.loglogFit)[1]), color = 'k')
            # plt.show()
            
            # ddd =1
            
            if fitElement.loglogSlope >= beta_threshold:
                jumping_runs.append(runC-1)
                continue
            else:
                exponents_list.append(fitElement.loglogSlope)
                # beta = fitElement.loglogSlope
                # if (fitElement.loglogSlope>1.01):
                #     plt.figure()
                #     plt.plot(fitElement.x, fitElement.y)
                #     plt.plot(fitElement.x, (fitElement.y[0]/fitElement.x[0])*fitElement.x, linestyle='dashed', color='k')
                # plt.xscale("log")
                # plt.yscale("log")
                
                fitElement = fitClass()
                d = 2 #dimension
                fitElement.x = 2*d*(1.0*time[t_begin_fit_ind:]-1.0*t_w)
                fitElement.y = np.mean(msd_block, axis=0)[t_begin_fit_ind:]
                fitElement.linlinFit = np.polyfit(fitElement.x, fitElement.y, 1)
                fitElement.linlinSlope = (fitElement.linlinFit)[0]
                diff_list.append(fitElement.linlinSlope)
            
            # plt.figure()
            # plt.plot(fitElement.x, fitElement.y)
            
        # print("block_size_list: "+str(block_size_list))
        # print("n_blocks: "+str(n_blocks)+", exp: "+str(np.sum(exponents_list*np.array(block_size_list)/np.sum(block_size_list))) +r"$ \pm $" + str(np.std(exponents_list)))
    
    
    np.savetxt('all_fitting_slopes.csv', X=np.array(all_fitting_slopes), delimiter=',', fmt='%.3f')
    np.savetxt('jumping_runs.csv', X=np.array(jumping_runs), delimiter=',', fmt='%d')
    
    return n_blocks-len(jumping_runs), block_size_list, exponents_list, diff_list

def fittingFunc():
    
    try:
        fitting_begin_array = np.loadtxt('../fitting_begin_sq.csv', delimiter=',')
    except:
        try:
            fitting_begin_array = np.loadtxt('../fitting_begin_irr.csv', delimiter=',')
        except:
            fitting_begin_array = np.loadtxt('../fitting_begin_hex.csv', delimiter=',')
        
    
    fitting_begin_dict = dict()
    for alphaC in range(np.shape(fitting_begin_array)[1]):
        fitting_begin_dict[fitting_begin_array[0,alphaC]] = fitting_begin_array[1,alphaC]
    
    
    for n_blocks in [5]: # This list MUST have just one element
    
        if n_runs%n_blocks ==0:
            block_size_list = n_blocks*[int(n_runs/n_blocks)]
        else:
            block_size_list = n_blocks*[0]
            counter = n_runs
            i=0
            while counter>0:
                block_size_list[i%n_blocks] += 1
                i += 1
                counter -= 1
            # block_size_typical = round(n_runs/n_blocks)
            # block_size_list = []
            # while counter >  block_size_typical+1:
            #     block_size_list.append(block_size_typical)
            #     counter = counter - block_size_typical
            # block_size_list.append(counter)
        
        exponents_list = []
        diff_list = []
        
        # if Alpha < 1.8 + (1e-10):
        #     t_begin_fit = t_w + 1
        # elif Alpha < 2.6 + (1e-10):
        #     t_begin_fit = t_w + 1
        # else:
        #     t_begin_fit = t_w + 1
            
        
        t_begin_fit = t_w + fitting_begin_dict[Alpha]
        
        t_begin_fit_ind = t_w_ind
        while time[t_begin_fit_ind] <= t_begin_fit:
            t_begin_fit_ind += 1
            
        
        
        runC=1
        for blockC in range(len(block_size_list)):
            block_size = block_size_list[blockC]
            
            msd_block = np.zeros([block_size, totalSnapshotsNum])
            
            for runC_dummy in range(block_size):
                
                pp_data_address = "run_"+str(runC)+"/pp_data"
                msdAvgStd = np.loadtxt(pp_data_address +'/msdAvgStd.csv', delimiter=',')
                
                msd_block[runC_dummy,:] = msdAvgStd[0,:]
                
                runC = runC+1
            
            fitElement = fitClass()
            
            fitElement.x = time[t_begin_fit_ind:]-t_w
            fitElement.y = np.mean(msd_block, axis=0)[t_begin_fit_ind:]
            
            # fitElement.linFit = np.polyfit(fitElement.x, fitElement.y, 1)
            # fitElement.linSlope = (fitElement.linFit)[0]
            fitElement.loglogFit = np.polyfit(np.log10(fitElement.x[1:]), np.log10(fitElement.y[1:]), 1)
            fitElement.loglogSlope = (fitElement.loglogFit)[0]
            exponents_list.append(fitElement.loglogSlope)
            beta = fitElement.loglogSlope
            # if (fitElement.loglogSlope>1.01):
            #     plt.figure()
            #     plt.plot(fitElement.x, fitElement.y)
            #     plt.plot(fitElement.x, (fitElement.y[0]/fitElement.x[0])*fitElement.x, linestyle='dashed', color='k')
            # plt.xscale("log")
            # plt.yscale("log")
            
            fitElement = fitClass()
            d = 2 #dimension
            fitElement.x = 2*d*(1.0*time[t_begin_fit_ind:]-1.0*t_w)
            fitElement.y = np.mean(msd_block, axis=0)[t_begin_fit_ind:]
            fitElement.linlinFit = np.polyfit(fitElement.x, fitElement.y, 1)
            fitElement.linlinSlope = (fitElement.linlinFit)[0]
            diff_list.append(fitElement.linlinSlope)
            
            # plt.figure()
            # plt.plot(fitElement.x, fitElement.y)
            
        print("block_size_list: "+str(block_size_list))
        print("n_blocks: "+str(n_blocks)+", exp: "+str(np.sum(exponents_list*np.array(block_size_list)/np.sum(block_size_list))) +r"$ \pm $" + str(np.std(exponents_list)))
    
    
    
    return n_blocks, block_size_list, exponents_list, diff_list

def get_ready_for_the_run():
    
    pp_data_address = "run_"+str(runC+1)+"/pp_data"
    pp_plot_address = "run_"+str(runC+1)+"/pp_plot"
    pp_vor_address = "run_"+str(runC+1)+"/pp_vor"
    
    plot_switches_dict = dict()
    with open("run_"+str(runC+1)+'/plot_switches.txt') as plot_switches:
        line = True
        while line:
            line = plot_switches.readline()
            line_list = line.split(": ")
            try:
                plot_switches_dict[line_list[0]]=int(line_list[1][0])
            except IndexError:
                break
            
    # eq_indices = np.loadtxt(pp_vor_address+"/"+"eq_snapshots_indices.csv", delimiter=',', dtype=int)
    
    eqSamplingTimesFile = np.loadtxt(pp_data_address+'/eqSamplingTimes.txt', dtype=str)
    eqSamplingTimes = eqSamplingTimesClass()
    
    eqSamplingTimes.snapshotCList = []
    eqSamplingTimes.timesList = []
    eqSamplingTimes.sampleCList = []
    eqSamplingTimes.bunchCList = []
            
    lineC= 1
    while(1):
        try:
            eqSamplingTimes.snapshotCList.append(int(eqSamplingTimesFile[lineC][0]))
            eqSamplingTimes.timesList.append(int(eqSamplingTimesFile[lineC][1]))
            eqSamplingTimes.sampleCList.append(int(eqSamplingTimesFile[lineC][2]))
            eqSamplingTimes.bunchCList.append(int(eqSamplingTimesFile[lineC][3]))
            
            lineC += 1
        except:
            break
    eq_indices = np.array(eqSamplingTimes.snapshotCList, dtype=int)
    N_eq_samples = len(eqSamplingTimes.snapshotCList)
    
    return pp_data_address, pp_plot_address, pp_vor_address, plot_switches_dict, eqSamplingTimes, eq_indices, N_eq_samples

n_runs = 0
while (os.path.isdir('run_'+str(n_runs+1))):
    n_runs += 1

simulationData = np.loadtxt('simulationData_vec.csv', dtype=str)

try:
    
    NSites                      = int(  simulationData[0][-1])
    Lx                          = float(  simulationData[1][-1])
    Ly                          = float(  simulationData[2][-1])
    Alpha                       = float(simulationData[3][-1])
    Kb                          = float(simulationData[4][-1])
    Tem                         = float(simulationData[5][-1])
    NumCells                    = int(  simulationData[6][-1])
    AvgCellArea                 = float(simulationData[7][-1])
    Lambda                      = float(simulationData[8][-1])
    # SweepLength                 = int(  simulationData[7][-1])
    SweepLength                 = NSites
    MaxMCSteps                  = int(  simulationData[10][-1])
    # MaxMCSteps                  = int(np.loadtxt('maxMCSteps.csv', dtype=int))
    samplesPerWrite             = int(  simulationData[11][-1])
    
    lattice_switch = 'irr'
    
    xShift = 0.0
    yShift = 0.0
    
except:
    
    L                           = int(  simulationData[0][-1])
    Alpha                       = float(simulationData[1][-1])
    Kb                          = float(simulationData[2][-1])
    Tem                         = float(simulationData[3][-1])
    NumCells                    = int(  simulationData[4][-1])
    AvgCellArea                 = float(simulationData[5][-1])
    Lambda                      = float(simulationData[6][-1])
    # SweepLength                 = int(  simulationData[7][-1])
    SweepLength                 = L*L
    MaxMCSteps                  = int(  simulationData[8][-1])
    # MaxMCSteps                  = int(np.loadtxt('maxMCSteps.csv', dtype=int))
    samplesPerWrite             = int(  simulationData[9][-1])
    
    Lx = L
    Ly = L
    
    xShift = 0.5
    yShift = 0.5
    
    lattice_switch = 'sq'

time = np.loadtxt('run_1/pp_data/time.csv', dtype=int, delimiter=',')
t_w = int(np.loadtxt('run_1/t_w.csv', dtype=int, delimiter=','))

t_w_ind = 0
while time[t_w_ind] != t_w:
    t_w_ind += 1

totalSnapshotsNum = len(time)



# ############### fitting MSD ##################
# # n_blocks, block_size_list, exponents_list, diff_list = fittingFunc()
# beta_threshold = 1.2
# n_blocks, block_size_list, exponents_list, diff_list = fittingFuncSmart(beta_threshold)
# ############### fitting MSD ##################


polygons_AR, polygons_circ = polygonCircAR_calc()

AR_list_alpha = [] #
AR_Vor_list_alpha = []
circ_list_alpha = [] #
circVor_list_alpha = []
hex_list_alpha = [] #
q_list_alpha = [] #
q_inv_list_alpha = [] #
q_vor_list_alpha = []
q_inv_vor_list_alpha = []
# q_vor_PDF_alpha = []
edge_vor_PDF_alpha = [] #
distances_vor_PDF_alpha = [] #

AR_distObject = distributionClass()
AR_Vor_distObject = distributionClass()
circ_distObject = distributionClass()
circVor_distObject = distributionClass()
q_distObject = distributionClass()
q_inv_distObject = distributionClass()
hex_distObject = distributionClass()
qVor_distObject = distributionClass()
qInvVor_distObject = distributionClass()
edge_vor_distObject = distributionClass()
edge_vor_log_distObject = distributionClass()
distances_vor_distObject = distributionClass()

for runC in range(n_runs):
    
    pp_data_address, pp_plot_address, pp_vor_address, plot_switches_dict, eqSamplingTimes, eq_indices, N_eq_samples = get_ready_for_the_run()
    
    ################# AR, circ, q, q_inv, hex Distributions ##########################
    distObject = distributionFuncExtended("AR")
    if runC==0:
        AR_distObject.N_tot = distObject.N_tot
        AR_distObject.count = distObject.count.copy()
        AR_distObject.bins  = distObject.bins.copy()
        AR_distObject.dx    = distObject.dx
        AR_distObject.binCenters = distObject.binCenters
    else:
        AR_distObject.N_tot += distObject.N_tot
        AR_distObject.count += distObject.count.copy()
    del(distObject)
    
    distObject = distributionFuncExtended("circ")
    if runC==0:
        circ_distObject.N_tot = distObject.N_tot
        circ_distObject.count = distObject.count.copy()
        circ_distObject.bins  = distObject.bins.copy()
        circ_distObject.dx    = distObject.dx
        circ_distObject.binCenters = distObject.binCenters
    else:
        circ_distObject.N_tot += distObject.N_tot
        circ_distObject.count += distObject.count.copy()
    del(distObject)
    
    distObject = distributionFuncExtended("q")
    if runC==0:
        q_distObject.N_tot = distObject.N_tot
        q_distObject.count = distObject.count.copy()
        q_distObject.bins  = distObject.bins.copy()
        q_distObject.dx    = distObject.dx
        q_distObject.binCenters = distObject.binCenters
    else:
        q_distObject.N_tot += distObject.N_tot
        q_distObject.count += distObject.count.copy()
    del(distObject)
    
    distObject = distributionFuncExtended("q_inv")
    if runC==0:
        q_inv_distObject.N_tot = distObject.N_tot
        q_inv_distObject.count = distObject.count.copy()
        q_inv_distObject.bins  = distObject.bins.copy()
        q_inv_distObject.dx    = distObject.dx
        q_inv_distObject.binCenters = distObject.binCenters
    else:
        q_inv_distObject.N_tot += distObject.N_tot
        q_inv_distObject.count += distObject.count.copy()
    del(distObject)
    
    distObject = distributionFuncExtended("hex")
    if runC==0:
        hex_distObject.N_tot = distObject.N_tot
        hex_distObject.count = distObject.count.copy()
        hex_distObject.bins  = distObject.bins.copy()
        hex_distObject.dx    = distObject.dx
        hex_distObject.binCenters = distObject.binCenters
    else:
        hex_distObject.N_tot += distObject.N_tot
        hex_distObject.count += distObject.count.copy()
    del(distObject)
    ################# AR, circ, q, q_inv, hex Distributions ##########################
    
    ############################## AR_Vor and circVor distributions #####################################
    AR_VorCellsEq = np.zeros([NumCells+1, N_eq_samples])
    circVorCellsEq = np.zeros([NumCells+1, N_eq_samples])
    for eqSampleC in range(N_eq_samples):
        verticesCellsEq = verticesReader(pp_vor_address+"/"+"verticesCellsEq_"+str(eqSampleC)+".dat")
        
        VorCellInfoList = []
        VorCellInfoList.append(None)
        for cellInd in range(1,NumCells+1):
            points = verticesCellsEq[cellInd-1].copy()
            VorCellInfo = VorCellInfoCalc(points)
            
            VorCellInfoList.append(VorCellInfo)
            AR_VorCellsEq[cellInd, eqSampleC]  = VorCellInfo.AR_Vor
            circVorCellsEq[cellInd, eqSampleC] = VorCellInfo.circVor
        
        # plt.figure()
        # plt.hist(AR_VorCellsEq[1:, eqSampleC], bins=30)
        # plt.axvline(x=np.mean(AR_VorCellsEq[1:, eqSampleC]), linestyle='dashed', color='k')
        # plt.title("AR_Vor , mean="+str(np.mean(AR_VorCellsEq[1:, eqSampleC])))
        # plt.show()
        # plt.figure()
        # plt.hist(circVorCellsEq[1:, eqSampleC], bins=30)
        # plt.axvline(x=np.mean(circVorCellsEq[1:, eqSampleC]), linestyle='dashed', color='k')
        # plt.title("circVor , mean="+str(np.mean(circVorCellsEq[1:, eqSampleC])))
        # plt.show()
    np.savetxt(pp_vor_address+"/"+"AR_VorCellsEq.csv", X=AR_VorCellsEq, delimiter=',', fmt='%1.5f')
    np.savetxt(pp_vor_address+"/"+"circVorCellsEq.csv", X=circVorCellsEq, delimiter=',', fmt='%1.5f')
    
    lastBunchIdx = eqSamplingTimes.bunchCList[-1]
    xCom_instant = np.loadtxt(pp_data_address+"/xComBunch_"+str(lastBunchIdx)+".csv", delimiter=',')[:,-1]
    yCom_instant = np.loadtxt(pp_data_address+"/yComBunch_"+str(lastBunchIdx)+".csv", delimiter=',')[:,-1]
    ### AR_Vor DISTRIBUTION
    distObject = distributionFuncExtended("AR_Vor")
    if runC==0:
        AR_Vor_distObject.N_tot = distObject.N_tot
        AR_Vor_distObject.count = distObject.count.copy()
        AR_Vor_distObject.bins  = distObject.bins.copy()
        AR_Vor_distObject.dx    = distObject.dx
        AR_Vor_distObject.binCenters = distObject.binCenters
    else:
        AR_Vor_distObject.N_tot += distObject.N_tot
        AR_Vor_distObject.count += distObject.count.copy()
    del(distObject)
    
    cmap = plt.get_cmap('summer')
    VorCellsPlotter(VorCellInfoList, AR_VorCellsEq[:,-1], pp_plot_address+"/AR_Vor_final.PNG", 1, cmap, 1, "AR_Vor", xCom_instant, yCom_instant)
    ### AR_Vor DISTRIBUTION
    ### circVor DISTRIBUTION
    distObject = distributionFuncExtended("circVor")
    if runC==0:
        circVor_distObject.N_tot = distObject.N_tot
        circVor_distObject.count = distObject.count.copy()
        circVor_distObject.bins  = distObject.bins.copy()
        circVor_distObject.dx    = distObject.dx
        circVor_distObject.binCenters = distObject.binCenters
    else:
        circVor_distObject.N_tot += distObject.N_tot
        circVor_distObject.count += distObject.count.copy()
    del(distObject)
    
    cmap = plt.get_cmap('cool')
    VorCellsPlotter(VorCellInfoList, circVorCellsEq[:,-1], pp_plot_address+"/circVor_final.PNG", 1, cmap, 1, "cirVor", xCom_instant, yCom_instant)
    ### circVor DISTRIBUTION
    cmap = plt.get_cmap('rainbow')
    VorCellsPlotter(VorCellInfoList, np.array(list(range(NumCells+1))), pp_plot_address+"/final_config_Vor.PNG", 1, cmap, 1, "cells Indices", xCom_instant, yCom_instant)
    ############################## AR_Vor and circVor distributions #####################################
    
    AR_list_run = []
    AR_Vor_list_run = np.mean(AR_VorCellsEq[1:,:], axis=0)
    circ_list_run = []
    circVor_list_run = np.mean(circVorCellsEq[1:,:], axis=0)
    hex_list_run = []
    q_list_run = []
    q_inv_list_run = []
    
    ARAvgStd = np.loadtxt(pp_data_address +'/ARAvgStd.csv', delimiter=',')
    circAvgStd = np.loadtxt(pp_data_address + '/circAvgStd.csv', delimiter=',')
    snapshotsHexAbsAvg = np.loadtxt(pp_data_address + '/snapshotsHexAbsAvg.csv', delimiter=',')
    isoperiAvgStd = np.loadtxt(pp_data_address +'/isoperiAvgStd.csv', delimiter=',')
    
    ind = t_w_ind
    while (ind<totalSnapshotsNum):
        
        if ind in eq_indices:
            AR_list_run.append(ARAvgStd[0, ind])
            circ_list_run.append(circAvgStd[0, ind])
            hex_list_run.append(snapshotsHexAbsAvg[ind])
            q_list_run.append(isoperiAvgStd[0, ind])
            q_inv_list_run.append(2.0*np.sqrt(np.pi)/isoperiAvgStd[0, ind])
        
        ind = ind+1
    
    ############################## qVor and qInvVor distributions #####################################
    lastBunchIdx = eqSamplingTimes.bunchCList[-1]
    xCom_instant = np.loadtxt(pp_data_address+"/xComBunch_"+str(lastBunchIdx)+".csv", delimiter=',')[:,-1]
    yCom_instant = np.loadtxt(pp_data_address+"/yComBunch_"+str(lastBunchIdx)+".csv", delimiter=',')[:,-1]
    qCellEq = np.loadtxt(pp_vor_address+'/qCellEq.csv', delimiter=',')
    ### qVor DISTRIBUTION
    
    distObject = distributionFuncExtended("qVor")
    if runC==0:
        qVor_distObject.N_tot = distObject.N_tot
        qVor_distObject.count = distObject.count.copy()
        qVor_distObject.bins  = distObject.bins.copy()
        qVor_distObject.dx    = distObject.dx
        qVor_distObject.binCenters = distObject.binCenters
    else:
        qVor_distObject.N_tot += distObject.N_tot
        qVor_distObject.count += distObject.count.copy()
    del(distObject)
    cmap = plt.get_cmap('hot')
    VorCellsPlotter(VorCellInfoList, qCellEq[:,-1], pp_plot_address+"/qVor_final.PNG", 1, cmap, 1, "qVor", xCom_instant, yCom_instant)
    ### AR_Vor DISTRIBUTION
    ### circVor DISTRIBUTION
    distObject = distributionFuncExtended("qInvVor")
    if runC==0:
        qInvVor_distObject.N_tot = distObject.N_tot
        qInvVor_distObject.count = distObject.count.copy()
        qInvVor_distObject.bins  = distObject.bins.copy()
        qInvVor_distObject.dx    = distObject.dx
        qInvVor_distObject.binCenters = distObject.binCenters
    else:
        qInvVor_distObject.N_tot += distObject.N_tot
        qInvVor_distObject.count += distObject.count.copy()
    del(distObject)
    cmap = plt.get_cmap('hot_r')
    VorCellsPlotter(VorCellInfoList, 2*(np.pi**0.5)/qCellEq[:,-1], pp_plot_address+"/qInvVor_final.PNG", 1, cmap, 1, "q_Inv_Vor", xCom_instant, yCom_instant)
    ### circVor DISTRIBUTION
    ############################## qVor and qInvVor distributions #####################################
    q_vor_list_run = np.mean(qCellEq[1:,:], axis=0)
    q_inv_vor_list_run = np.mean(2.0*np.sqrt(np.pi)/qCellEq[1:,:], axis=0)
    

    AR_list_alpha.append(np.mean(AR_list_run))
    AR_Vor_list_alpha.append(np.mean(AR_Vor_list_run))
    circ_list_alpha.append(np.mean(circ_list_run))
    circVor_list_alpha.append(np.mean(circVor_list_run))
    hex_list_alpha.append(np.mean(hex_list_run))
    q_list_alpha.append(np.mean(q_list_run))
    q_inv_list_alpha.append(np.mean(q_inv_list_run))
    q_vor_list_alpha.append(np.mean(q_vor_list_run))
    q_inv_vor_list_alpha.append(np.mean(q_inv_vor_list_run))
    
    ################## edges length and distance distribution ################################
    distObject = distributionFuncExtended("edgesVor")
    if runC==0:
        edge_vor_distObject.N_tot = distObject.N_tot
        edge_vor_distObject.count = distObject.count.copy()
        edge_vor_distObject.bins  = distObject.bins.copy()
        edge_vor_distObject.dx    = distObject.dx
        edge_vor_distObject.binCenters = distObject.binCenters
    else:
        edge_vor_distObject.N_tot += distObject.N_tot
        edge_vor_distObject.count += distObject.count.copy()
    del(distObject)
    
    distObject = distributionFuncExtended("edgesVorLog")
    if runC==0:
        edge_vor_log_distObject.N_tot = distObject.N_tot
        edge_vor_log_distObject.count = distObject.count.copy()
        edge_vor_log_distObject.bins  = distObject.bins.copy()
        edge_vor_log_distObject.dx    = distObject.dx
        edge_vor_log_distObject.binCenters = distObject.binCenters
    else:
        edge_vor_log_distObject.N_tot += distObject.N_tot
        edge_vor_log_distObject.count += distObject.count.copy()
    del(distObject)
    
    distObject = distributionFuncExtended("distancesVor")
    if runC==0:
        distances_vor_distObject.N_tot = distObject.N_tot
        distances_vor_distObject.count = distObject.count.copy()
        distances_vor_distObject.bins  = distObject.bins.copy()
        distances_vor_distObject.dx    = distObject.dx
        distances_vor_distObject.binCenters = distObject.binCenters
    else:
        distances_vor_distObject.N_tot += distObject.N_tot
        distances_vor_distObject.count += distObject.count.copy()
    del(distObject)
    ################## edges length and distance distribution ################################
    
    
    print(runC)



####################### PDF's plot ##############################
PDF_creator(AR_distObject, "AR_PDF.csv", 'AR')
PDF_creator(circ_distObject, "circ_PDF.csv", 'circ')
PDF_creator(q_distObject, "q_PDF.csv", 'q')
PDF_creator(q_inv_distObject, "q_inv_PDF.csv", 'q_inv')
PDF_creator(hex_distObject, "hex_PDF.csv", 'hex')

PDF_creator(AR_Vor_distObject, "AR_Vor_PDF.csv", 'AR_Vor')
PDF_creator(circVor_distObject, "circVor_PDF.csv", 'circVor')
PDF_creator(qVor_distObject, "qVor_PDF.csv", 'qVor')
PDF_creator(qInvVor_distObject, "qInvVor_PDF.csv", 'qInvVor')
PDF_creator(edge_vor_distObject, "edges_vor_PDF.csv", 'edgesVor')
PDF_creator(edge_vor_log_distObject, "edges_vor_log_PDF.csv", 'edgesVorLog')
PDF_creator(distances_vor_distObject, "distances_vor_PDF.csv", 'distancesVor')
####################### PDF's plot ##############################


fancy_qVor_PDF_plotter() #showing perfect polygons q






####################### writing summary file ##############################
f = open("pp_summary.txt", "w")
f.write("Alpha = "+str(Alpha)+"\n")
f.write("n_runs = "+str(n_runs)+"\n")
# f.write("n_blocks = "+str(n_blocks)+"\n")
# f.write("block_size_list = "+str(block_size_list)+"\n")
f.write("----------------------------------------------\n")

# circ_list_alpha.append(np.mean(circ_list_run))
# hex_list_alpha.append(np.mean(hex_list_run))
# q_list_alpha.append(np.mean(q_list_run))
# q_vor_list_alpha.append(np.mean(q_vor_list_run))

f.write("equilibrium AR        =\t"+str(np.mean(AR_list_alpha)) +      " \pm " +str(np.std(AR_list_alpha))     + "\n")
f.write("equilibrium AR_Vor    =\t"+str(np.mean(AR_Vor_list_alpha)) +  " \pm " +str(np.std(AR_Vor_list_alpha)) + "\n")
f.write("equilibrium circ      =\t"+str(np.mean(circ_list_alpha)) +    " \pm " +str(np.std(circ_list_alpha)) + "\n")
f.write("equilibrium circVor   =\t"+str(np.mean(circVor_list_alpha)) + " \pm " +str(np.std(circVor_list_alpha)) + "\n")
f.write("equilibrium hex   =\t"+str(np.mean(hex_list_alpha)) + " \pm " +str(np.std(hex_list_alpha)) + "\n")
f.write("equilibrium q     =\t"+str(np.mean(q_list_alpha)) + " \pm " +str(np.std(q_list_alpha)) + "\n")
f.write("equilibrium q_inv     =\t"+str(np.mean(q_inv_list_alpha)) + " \pm " +str(np.std(q_inv_list_alpha)) + "\n")
f.write("equilibrium q_vor =\t"+str(np.mean(q_vor_list_alpha)) + " \pm " +str(np.std(q_vor_list_alpha)) + "\n")
f.write("equilibrium q_inv_vor =\t"+str(np.mean(q_inv_vor_list_alpha)) + " \pm " +str(np.std(q_inv_vor_list_alpha)) + "\n")
f.write("----------------------------------------------\n")

# f.write("exponents_list  = "+str(exponents_list)+"\n")
# f.write("diff_coef_list  = "+str(diff_list)+"\n")

f.close()
####################### writing summary file ##############################





