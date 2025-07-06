#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 26 12:10:16 2023

@author: hossein
"""


import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as tri
import matplotlib.patches as patches
import freud
import time as TIME
import zipfile
import shutil
import sys
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection
from scipy.spatial import ConvexHull
import math


class eqSamplingTimesClass:
    snapshotCList = []
    timesList = []
    sampleCList = []
    bunchCList = []
    

def read_comparts_info(filepath):
    """
    Reads a file defining a vector and two symmetric matrices in a block format.

    The file format must have three blocks:
      - AvgAreaFrac: a column of floats ending with a line containing '#'
      - J_int: lower-triangular rows of a symmetric matrix ending with '#'
      - J_ext: same as J_int for another matrix

    Parameters
    ----------
    filepath : str
        Path to the text file.

    Returns
    -------
    tuple
        - vector : np.ndarray of shape (n,)
        - J_int  : np.ndarray of shape (n, n)
        - J_ext  : np.ndarray of shape (n, n)
    """
    # Parse blocks
    blocks = {}
    with open(filepath, 'r') as f:
        lines = [line.strip() for line in f if line.strip()]

    i = 0
    while i < len(lines):
        name = lines[i]
        i += 1
        data = []
        while i < len(lines) and lines[i] != '#':
            parts = lines[i].split()
            data.append([float(p) for p in parts])
            i += 1
        # Skip the '#' delimiter
        i += 1
        blocks[name] = data

    # Extract vector
    vector = np.array([row[0] for row in blocks['AvgAreaFrac']], dtype=float)

    # Helper to build full symmetric matrix from lower-triangular rows
    def build_symmetric(lower_rows):
        n = len(lower_rows)
        M = np.zeros((n, n), dtype=float)
        for r, row in enumerate(lower_rows):
            for c, val in enumerate(row):
                M[r, c] = val
                M[c, r] = val
        return M

    J_int = build_symmetric(blocks['J_int'])
    J_ext = build_symmetric(blocks['J_ext'])
    
    avgAreaFrac = vector.copy()
    
    return avgAreaFrac, J_int, J_ext

ppDataFolderName = 'pp_data'
ppPlotFolderName = 'pp_plot'

simulationData = np.loadtxt('simulationData_vec.csv', dtype=str)

try:
    
    NSites                      = int(  simulationData[0][-1])
    Lx                          = float(  simulationData[1][-1])
    Ly                          = float(  simulationData[2][-1])
    # Alpha                       = float(simulationData[3][-1])
    Kb                          = float(simulationData[3][-1])
    Tem                         = float(simulationData[4][-1])
    NumCells                    = int(  simulationData[5][-1])
    AvgCellArea                 = float(simulationData[6][-1])
    Lambda                      = float(simulationData[7][-1])
    # SweepLength                 = int(  simulationData[7][-1])
    SweepLength                 = NSites
    MaxMCSteps                  = int(  simulationData[9][-1])
    # MaxMCSteps                  = int(np.loadtxt('maxMCSteps.csv', dtype=int))
    samplesPerWrite             = int(  simulationData[10][-1])
    
    lattice_switch = 'irr'
    
    avgAreaFrac, J_int, J_ext = read_comparts_info('compart_params.txt')
    
    avgArea = AvgCellArea * avgAreaFrac
    
except:
    
    fff = 1
    # L                           = int(  simulationData[0][-1])
    # Alpha                       = float(simulationData[1][-1])
    # Kb                          = float(simulationData[2][-1])
    # Tem                         = float(simulationData[3][-1])
    # NumCells                    = int(  simulationData[4][-1])
    # AvgCellArea                 = float(simulationData[5][-1])
    # Lambda                      = float(simulationData[6][-1])
    # # SweepLength                 = int(  simulationData[7][-1])
    # SweepLength                 = L*L
    # MaxMCSteps                  = int(  simulationData[8][-1])
    # # MaxMCSteps                  = int(np.loadtxt('maxMCSteps.csv', dtype=int))
    # samplesPerWrite             = int(  simulationData[9][-1])
    
    # Lx = L
    # Ly = L
    
    # lattice_switch = 'sq'


# equilibration sampling times
eqSamplingTimesFile = np.loadtxt(ppDataFolderName+'/eqSamplingTimes.txt', dtype=str)
eqSamplingTimes = eqSamplingTimesClass()
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





numberOfSamplesPerWindow = 10





try:
    os.mkdir(ppPlotFolderName) 
except FileExistsError:
    pass

mainResumeFolderName = 'main_resume'
backupResumeFolderName = 'backup_resume'

initFolderName   = 'init' 
dataFolderName   = 'data' 

n_digits = 3



totalSnapshotsNum = np.loadtxt(ppDataFolderName+'/totalSnapshotsNum.csv', dtype=int, delimiter=',')
t_w = np.loadtxt('t_w.csv', delimiter=',', dtype=int)
WaitingMCSteps = t_w

final_sampling_numbers = int(0.2*totalSnapshotsNum)

def distributionFunc(quantity_name):
    
    xlabel_description = None
    if quantity_name == 'AR':
        dist_eq = np.loadtxt(ppDataFolderName+"/AR_dist_eq.csv", delimiter=',')[1:, :]
        dx = 0.01
        min_val = 1.0
        # max_val = 3.6
        max_val = 4.0
        
        write_file_name = ppDataFolderName+"/AR_freq.txt"
        this_plot_switch = plot_switches_dict[quantity_name]
        
    elif quantity_name == 'q':
        dist_eq = np.loadtxt(ppDataFolderName+"/isoperi_dist_eq.csv", delimiter=',')[1:, :]
        dx = 0.01
        min_val = 0.0
        # max_val = 3.6
        max_val = 1.0
        write_file_name = ppDataFolderName+"/isoperi_freq.txt"
        this_plot_switch = plot_switches_dict[quantity_name]
        
    elif quantity_name == 'q_inv':
        dist_eq = np.loadtxt(ppDataFolderName+"/isoperi_dist_eq.csv", delimiter=',')[1:, :]
        dist_eq = 2*(np.pi**0.5)/dist_eq
        dx = 0.01
        min_val = 2*(np.pi**0.5)
        max_val = 10.0
        write_file_name = ppDataFolderName+"/isoperi_inv_freq.txt"
        this_plot_switch = plot_switches_dict['q']
        xlabel_description = '('+r'$P/\sqrt(A)$'+')'
    
    elif quantity_name == 'circ':
        dist_eq = np.loadtxt(ppDataFolderName+"/circ_dist_eq.csv", delimiter=',')[1:, :]
        dx = 0.01
        min_val = 0.0
        max_val = 1.0
        write_file_name = ppDataFolderName+"/circ_freq.txt"
        this_plot_switch = plot_switches_dict['circ']
    
    elif quantity_name == 'hex':
        dist_eq = np.loadtxt(ppDataFolderName+"/hex_dist_eq.csv", delimiter=',')[1:, :]
        dx = 0.01
        min_val = 0.0
        max_val = 1.0
        write_file_name = ppDataFolderName+"/hex_freq.txt"
        this_plot_switch = plot_switches_dict['order']
        
    bin_list = []
    bin_list.append(min_val-dx/2)
    tracker = min_val + dx/2
    while (tracker <= max_val):
        bin_list.append(round(tracker, 4))
        tracker = tracker + dx
    bin_list.append(tracker)
    bin_list.append(tracker+dx)
    
    nBins = len(bin_list)-1
    
    total_data = []
    for eqSampleC in range(np.shape(dist_eq)[1]):
        total_data = total_data + list(dist_eq[1:,eqSampleC])
        
    freq_hist = np.histogram(total_data, bins = bin_list)
    # AR_PDF_all[eqSampleC,:] = AR_PDF_element[0]/(dAR * np.sum(AR_PDF_element[0]) )
    
        # plt.figure()
        # plt.hist(hex_data, bins=bin_list)
        # plt.axvline(x=np.mean(hex_data), linewidth=1, linestyle='dashed', c ='black')
        # plt.show()
    freq_hist_err = (1.0*freq_hist[0])**0.5
    
    # AR_normal_freq_avg = (1.0 * AR_PDF_element[0])/len(AR_data)
    # AR_normal_freq_err = (1.0 * AR_PDF_element[0]**0.5)/len(AR_data)
    
    # f = open(ppDataFolderName+"/AR_PDF.txt", "w")
    # f.write("mean = "+str(np.mean(AR_PDF_all, axis=0))+"\n")
    # f.write("----------------------------------------------\n")
    # f.write("std  = "+str(np.std(AR_PDF_all, axis=0))+"\n")
    # f.write("----------------------------------------------\n")
    # f.write("bins = "+str(bin_list))
    # f.close()
    
    # f = open(ppDataFolderName+"/AR_PDF.txt", "w")
    f = open(write_file_name, "w")
    f.write("mean = "+str(1.0*freq_hist[0]/(np.shape(dist_eq)[1]*NumCells))+"\n")
    f.write("----------------------------------------------\n")
    f.write("std  = "+str(1.0*freq_hist_err/np.shape(dist_eq)[1]*NumCells)+"\n")
    f.write("----------------------------------------------\n")
    f.write("bins = "+str(bin_list))
    f.close()

    if this_plot_switch:
        plt.figure()
        x_plot = np.array(bin_list[1:])-0.5*dx
        # # y_plot = np.mean(AR_PDF_all, axis=0)
        # # y_err = np.std(AR_PDF_all, axis=0)/np.sqrt(np.shape(AR_PDF_all)[0])
        # y_plot = AR_normal_freq_avg
        # y_err = AR_normal_freq_err
        # # y_err = np.std(AR_PDF_all, axis=0)
        # plt.plot(x_plot, y_plot)
        # plt.fill_between(x_plot, y_plot-y_err, y_plot+y_err, alpha=0.5)
        y_plot = freq_hist[0]/(np.shape(dist_eq)[1]*NumCells)
        menStd     = freq_hist_err/(np.shape(dist_eq)[1]*NumCells)
        width      = dx
        plt.bar(x_plot, y_plot, width=width, yerr=menStd)

        # plt.hist(AR_data, bins = bin_list)
        plt.xlim([np.min(total_data)-3*dx,np.max(total_data)+3*dx])
        plt.xlabel(quantity_name)
        plt.ylabel("freq("+quantity_name+")")
        plt.axvline(x=np.mean(total_data), linestyle='dashed', color='k', linewidth=1)
        plt.axvline(x=min_val, linestyle='dashed', color='k', linewidth=0.5)
        plt.axvline(x=max_val, linestyle='dashed', color='k', linewidth=0.5)
        plt.axhline(y=0, linestyle='dashed', color='k', linewidth=1)
        # plt.title(r"$\alpha = $"+str(Alpha)+"\n Avg: "+str(round(np.mean(total_data), 3))+", N_tot = "+str(np.shape(dist_eq)[1])+"*"+str(NumCells)+" = "+str(np.shape(dist_eq)[1]*NumCells))
        plt.savefig(ppPlotFolderName+"/"+quantity_name+"_freq_eq.PNG", dpi = 300)
        # plt.show()
    
    
    return

def g_of_r_calc(x_set, y_set, Lx, Ly, nBins, rCutOff, saving_name):
    
    dr = rCutOff / nBins
    N_cells = len(x_set) # not the same as NumClles. Because this function is going to be general.
    
    g_r_cells = np.zeros([N_cells, nBins])
    
    x_set = x_set % Lx
    y_set = y_set % Ly

    
    for cellC_1 in range(N_cells):
        
        copies = dict()
        copies['x'] = [0.0]
        copies['y'] = [0.0]
        
        if x_set[cellC_1] <= rCutOff:
            copies['x'].append(-Lx)
        if x_set[cellC_1] >= Lx-rCutOff:
            copies['x'].append(+Lx)
        if y_set[cellC_1] <= rCutOff:
            copies['y'].append(-Ly)
        if y_set[cellC_1] >= Ly-rCutOff:
            copies['y'].append(+Ly)
        
        r_data = []
        
        for cellC_2 in range(N_cells):
            
            if cellC_2 == cellC_1:
                continue
            
            for dx_copy in copies['x']:
                for dy_copy in copies['y']:
                    
                    Dx = x_set[cellC_1]-(x_set[cellC_2]+dx_copy)
                    Dy = y_set[cellC_1]-(y_set[cellC_2]+dy_copy)
                    
                    r = np.sqrt(Dx**2 + Dy**2)
                    
                    if r<=rCutOff:
                        r_data.append(r)
        
        hist_data = np.histogram(r_data, bins = np.arange(0, rCutOff+0.0001*dr, dr))
        areaElement = np.pi*((hist_data[1][1:])**2 - (hist_data[1][0:nBins])**2)
        plot_y = (hist_data[0] / areaElement) / (N_cells / (Lx*Ly))
        g_r_cells[cellC_1,:]= plot_y 
        # plot_x = hist_data[1][0:nBins]+dr/2
        # plt.plot(plot_x, plot_y)
        # plt.show()
        # plt.hist(r_data, bins = np.arange(0, rCutOff+0.0001*dr, dr))
        # plt.show()
        
        # ddd=5
    
    # g_r = np.mean(g_r_cells, axis=0)
    if saving_name:
        plt.figure()
        x_plot = np.arange(0+0.5*dr , rCutOff+0.0001*dr , dr)/ np.sqrt(AvgCellArea)
        y_plot = np.mean(g_r_cells, axis=0)
        y_plot_err = np.std(g_r_cells, axis=0)/np.sqrt(N_cells)
        plt.errorbar(x_plot, y_plot, yerr=y_plot_err, fmt="o", markersize =0.8, elinewidth=0.2, label='g(r)', zorder=1)
        plt.xlabel(r"$r/\sqrt{A_0}$")
        plt.ylabel(r"$g(r)$")
        plt.axhline(y=1, linewidth=1, linestyle='dashed', c ='black')
        plt.legend()
        plt.grid()
        # plt.savefig("g_of_r_"+str(eqSampleC)+".PNG", dpi=500)
        plt.savefig(saving_name, dpi=500)
        plt.show()
    
    g_r = np.mean(g_r_cells, axis=0)
    return g_r

def neigh_reader():
    address = '../lattice/neighs.dat'
    neigh_dict = dict()
    with open(address) as neigh_file:
        line = True
        while line:
            line = neigh_file.readline()
            line_list = line.split("\t")
            
            try:
                label = int(line_list[0])
                neigh_dict[label] = []
            except:
                return neigh_dict
            
            
            if line_list[-1] == '\n':
                ind_max = len(line_list)-2
            else:
                ind_max = len(line_list)-1
            
            ind = 3
            while (ind<=ind_max):
                neigh_dict[label].append(int(line_list[ind]))
                ind+=1
          
            # print(neigh_dict[label])
            # for line_element in line_list:
            #     if float(line_element) == int(line_element)
    return neigh_dict

if lattice_switch=='irr':
    neigh_dict = neigh_reader()

def xyComMaker():
    
    xCom = np.loadtxt(ppDataFolderName+'/'+'xComInit.csv', delimiter=',')
    yCom = np.loadtxt(ppDataFolderName+'/'+'yComInit.csv', delimiter=',')
    
    xCom = xCom.reshape([NumCells+1,1])
    yCom = yCom.reshape([NumCells+1,1])
    
    bunchC=0
    while(1):
        try:
            xCom_concat = np.loadtxt(ppDataFolderName+'/'+'xComBunch_'+str(bunchC)+'.csv', delimiter=',')
            yCom_concat = np.loadtxt(ppDataFolderName+'/'+'yComBunch_'+str(bunchC)+'.csv', delimiter=',')
            
            xCom = np.concatenate([xCom,xCom_concat], 1)
            yCom = np.concatenate([yCom,yCom_concat], 1)
            bunchC+=1
        except:
            break
    
    return xCom, yCom

def cumulativeAveraging(x, y):
    
    yCumAvg = 0*y
    yCumAvg[0] = y[0]
    
    integral = 0
    
    for i in range(1,len(yCumAvg)):
        integral += (x[i]-x[i-1])*(y[i]+y[i-1])/2
        yCumAvg[i] = integral/(x[i]-x[0])
        
    
    return yCumAvg

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

def cellsPlotterCSV(sigmaMat, savePlotAddress, saveSwitch):
    
    dpiRes = 500
    
    #sigmaMat = np.loadtxt(loadDataAddress, dtype=int, delimiter=',')
    # levels = np.arange(np.min(sigmaMat)+0.5,np.max(sigmaMat)-0.5, 1)
    # plt.imshow(sigmaMat, cmap='tab20')
    # plt.contour(sigmaMat, levels, colors='k')
    
    # plt.imshow(sigmaMat, cmap='viridis')
    # plt.contour(sigmaMat, colors='black', linewidths=0.5)
    
    # cmap = plt.get_cmap('hsv', np.max(sigmaMat) - np.min(sigmaMat) + 1)
    # cmap = plt.get_cmap('terrain', np.max(sigmaMat) - np.min(sigmaMat) + 1)
    cmap = plt.get_cmap('rainbow', np.max(sigmaMat) - np.min(sigmaMat) + 1)
    
    plt.figure()
    # Plot the image with discrete integer values and filled regions
    plt.imshow(sigmaMat, cmap=cmap, interpolation='nearest')
    
    # Draw contour lines on top of the filled regions
    # plt.contour(sigmaMat, levels=np.unique(sigmaMat), colors='black', linewidths=0.05, antialiased=False, corner_mask=False)
    plt.colorbar()
    
    lines_list = []
    L_ver = np.shape(sigmaMat)[0]
    L_hor = np.shape(sigmaMat)[1]
    
                
    verBorders = np.zeros([L_ver  , L_hor+1])
    for i in range(L_ver):
        for j in range(L_hor):
            if (sigmaMat[i,j]!=sigmaMat[i,(j-1)%L_hor]):
                verBorders[i,j]=1
                x_beg = j-0.5; y_beg=i-0.5;
                x_end = x_beg; y_end=y_beg+1;
                lines_list.append(np.array([[x_beg, y_beg],[x_end, y_end]]))
                
        j = L_hor
        if (sigmaMat[i,-1]!=sigmaMat[i,0]):
            verBorders[i,j]=1
            x_beg = (L_hor-1)+0.5; y_beg=i-0.5;
            x_end = x_beg;         y_end=y_beg+1;
            lines_list.append(np.array([[x_beg, y_beg],[x_end, y_end]]))


    horBorders = np.zeros([L_ver+1, L_hor  ])
    for i in range(L_ver):
        for j in range(L_hor):
            if (sigmaMat[i,j]!=sigmaMat[(i-1)%L_ver,j]):
                horBorders[i,j]=1
                x_beg = j-0.5;   y_beg=i-0.5;
                x_end = x_beg+1; y_end=y_beg;
                lines_list.append(np.array([[x_beg, y_beg],[x_end, y_end]]))
                
    i = L_ver
    for j in range(L_hor):
        if (sigmaMat[-1,j]!=sigmaMat[0,j]):
            horBorders[i,j]=1
            x_beg = j-0.5;   y_beg=(L_ver-1)+0.5;
            x_end = x_beg+1; y_end=y_beg;
            lines_list.append(np.array([[x_beg, y_beg],[x_end, y_end]]))

    # lines_list = [np.array([[20, 10],[40, 50]])]
    
    for line_c in range(len(lines_list)):
        line_beg = lines_list[line_c][0,:].copy()
        line_end = lines_list[line_c][1,:].copy()
        x_plot = [line_beg[0], line_end[0]]
        y_plot = [line_beg[1], line_end[1]]
        plt.plot(x_plot, y_plot, 'k', linewidth=0.2)
    
    plt.xlabel('Y')
    plt.ylabel('X')
    plt.title('color indicates cell index')
    
    if saveSwitch:
        #file_name_pic = fileName+".png"
        plt.savefig(savePlotAddress, dpi=dpiRes)
        
    return

def cellsPlotterIrr(sigmaMat, compartMat,  savePlotAddress, saveSwitch, cmap, com_switch, xCom_instant, yCom_instant):
    
    plt.figure()
    
    with open('../lattice/vertices.dat', 'r') as vertices_file:
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
                vertices_dict[index].append((float(line.split("\t")[1]), float(line.split("\t")[0])))
            else:
                continue
    
    polygons = []
    values = []
    for i in range(NSites):
        polygons.append(np.array(vertices_dict[i]))
        values.append(sigmaMat[i])
        # values.append(i+1)
    
    # cmap = plt.get_cmap('rainbow')
    if np.min(sigmaMat)==1 or np.min(sigmaMat)==0:
        vmin = 0
        vmax = NumCells
        lineColor = 'k'
    else:
        vmin = 0.0
        vmax = 1.0
        lineColor = 'w'
    norm = plt.Normalize(vmin, vmax)
    
    copy_direction_list = [[-1,-1], [-1, 0], [-1,+1], [0, +1], [+1, +1], [+1, 0], [+1, -1], [0,-1] ]
    
    length_scale = np.sqrt(Lx*Ly/NSites)
    margin = 2*length_scale
    edge_plot_check = [] #a list of sets containing the indices of sites
    for i in range(NSites):
        if compartMat[i] == 0:
            color = cmap(norm(values[i]))
        elif compartMat[i] == 1:
            color = 'grey'
        elif compartMat[i] == 2:
            color = 'yellow'
        main_points = polygons[i]
        points = main_points.copy()
        hull = ConvexHull(points)    
        convex_hull = graham_scan(points)
        plt.fill(*zip(*convex_hull), color=color, alpha=1, edgecolor='black', linewidth=0.0, zorder=1)
        
        edge_plot_list = []
        for neigh in neigh_dict[i]:
            if (sigmaMat[i] != sigmaMat[neigh]) and (set({neigh, i}) not in edge_plot_check):
                edgePlot = [[],[]]  #[[y1, y2], [x1, x2]]
                for vertice_i in vertices_dict[i]:
                    for vertice_neigh in vertices_dict[neigh]:
                        error = ((vertice_i[0]-vertice_neigh[0])**2+(vertice_i[1]-vertice_neigh[1])**2)**0.5 # r = pol(dx, dy)
                        if error< (1e-6)*length_scale:
                            edgePlot[0].append(vertice_i[0]) #this is y in the simulation
                            edgePlot[1].append(vertice_i[1]) #this is x in the simulation
                
                edge_plot_list.append(edgePlot)
            else:
                pass
        
        for edgePlot in edge_plot_list:
            plt.plot(edgePlot[0], edgePlot[1], linewidth=0.1, color=lineColor, zorder=1)
        
        
        for copy_counter in range(0,8):
            
            copy_dir = copy_direction_list[copy_counter]
            
            points[:,0] = main_points[:,0].copy() + Ly*copy_dir[1]
            points[:,1] = main_points[:,1].copy() + Lx*copy_dir[0]
            
            condition = (abs(points[0,0]-Ly/2)<Ly/2+margin) and (abs(points[0,1]-Lx/2)<Lx/2+margin) 
            
            if condition:
                hull = ConvexHull(points)
                convex_hull = graham_scan(points)
                plt.fill(*zip(*convex_hull), color=color, alpha=1, edgecolor='black', linewidth=0.0, zorder=1)
                
                for edgePlot in edge_plot_list:
                    plt.plot(np.array(edgePlot[0])+ Ly*copy_dir[1], np.array(edgePlot[1])+ Lx*copy_dir[0], linewidth=0.1, color=lineColor, zorder=1)
        
        # for j in range():
    
        
        
    # plt.figure()
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
    
    if com_switch:
        plt.scatter(yCom_instant[1:]%Ly, xCom_instant[1:]%Lx, s=0.1, c=lineColor, zorder=2)
    
    plt.savefig(savePlotAddress, dpi=800)
    # plt.show()
    
    return

def temporalPlotsV2(saveSwitch):
    
    function_time_start = TIME.time()
    
    
    
    dpiRes = 500
    
    
    ############################ ENERGY PLOT ##########################
    
    # time = data.t
    # energy = data.e
    
    time = np.loadtxt(ppDataFolderName+'/time.csv',  delimiter=',', dtype=int)
    energy = np.loadtxt(ppDataFolderName+'/energy.csv',  delimiter=',', dtype=float)
    
    eTemporalCumAvg = cumulativeAveraging(time, energy)
    np.savetxt(ppDataFolderName+'/'+'energyCumAvg.csv', eTemporalCumAvg, fmt='%.2f')
        
    # np.savetxt(ppDataFolderName+'/'+'time.csv', time, fmt='%d')
    
    # x_plot = time[1:]
    # y_plot = energy[1:]
    
    ind=0
    t = time[ind]
    while(t!=t_w):
        t = time[ind]
        ind+=1
    ind-=1
    t_w_ind = ind
    
    x_plot = time
    y_plot = energy
    
    if plot_switches_dict['temporals']:

        plt.figure()
        plt.plot(x_plot, y_plot, linewidth=0.6, label='instantaneous')
        plt.axvline(x=WaitingMCSteps, linewidth=1, linestyle='dashed', c ='black')
        plt.xlabel('MCS')
        plt.ylabel('Energy')
        
        # eTemporalCumAvg = 0*eTemporal
        # for t_c in range(len(eTemporal)):
        #     eTemporalCumAvg[t_c] = np.mean(eTemporal[0:t_c+1])
        
        # eTemporalCumAvg = cumulativeAveraging(time, energy)
        eTemporalCumAvg = cumulativeAveraging(x_plot, y_plot)
        
        # np.savetxt(ppDataFolderName+'/'+'energy.csv', energy, fmt='%.2f')
        # np.savetxt(ppDataFolderName+'/'+'energyCumAvg.csv', eTemporalCumAvg, fmt='%.2f')
        
        x_plot = time
        y_plot = eTemporalCumAvg
        plt.plot(x_plot, y_plot, linewidth=0.6, label='cumulative avg')
        plt.legend()
        plt.xscale("log")
        plt.grid(linewidth=0.3)
        
        if saveSwitch:
            plt.savefig(ppPlotFolderName+'/Energy.png',dpi=dpiRes)
        
        
        x_plot = time[t_w_ind:]
        y_plot = energy[t_w_ind:]
        
        plt.figure()
        plt.plot(x_plot, y_plot, linewidth=0.6, label='instantaneous')
        plt.axvline(x=WaitingMCSteps, linewidth=1, linestyle='dashed', c ='black')
        plt.xlabel('MCS')
        plt.ylabel('Energy')
        
        # eTemporalCumAvg = 0*eTemporal
        # for t_c in range(len(eTemporal)):
        #     eTemporalCumAvg[t_c] = np.mean(eTemporal[0:t_c+1])
        
        # eTemporalCumAvg = cumulativeAveraging(time, energy)
        eTemporalCumAvg = cumulativeAveraging(x_plot-WaitingMCSteps, y_plot)
        
        # np.savetxt(ppDataFolderName+'/'+'energy.csv', energy, fmt='%.2f')
        # np.savetxt(ppDataFolderName+'/'+'energyCumAvg.csv', eTemporalCumAvg, fmt='%.2f')
        
        x_plot = time[t_w_ind:]
        y_plot = eTemporalCumAvg[0:]
        plt.plot(x_plot, y_plot, linewidth=0.6, label='cumulative avg')
        plt.legend()
        plt.xscale("log")
        plt.grid(linewidth=0.3)
        
        if saveSwitch:
            plt.savefig(ppPlotFolderName+'/Energy_eq.png',dpi=dpiRes)
        
        
        
    
    ############################ ENERGY PLOT ##########################    
    
    ############################ KURTOSIS PLOT ##########################
    
    # time = data.t
    # energy = data.e
    
    time = np.loadtxt(ppDataFolderName+'/time.csv',  delimiter=',', dtype=int)
    kurt = np.loadtxt(ppDataFolderName+'/kurt.csv',  delimiter=',', dtype=float)
    
    # np.savetxt(ppDataFolderName+'/'+'time.csv', time, fmt='%d')
    
    if plot_switches_dict['temporals']:
        
        x_plot = time[1:]
        y_plot = kurt[1:]
        
        plt.figure()
        plt.plot(x_plot, y_plot, linewidth=0.6)
        plt.axvline(x=WaitingMCSteps, linewidth=1, linestyle='dashed', c ='black')
        plt.xlabel('MCS')
        plt.ylabel(r'$\alpha_2$')
        plt.title('Kurtosis '+r'$(\alpha_2)$')
        
        # eTemporalCumAvg = 0*eTemporal
        # for t_c in range(len(eTemporal)):
        #     eTemporalCumAvg[t_c] = np.mean(eTemporal[0:t_c+1])
        
        # eTemporalCumAvg = cumulativeAveraging(time, energy)
        
        # np.savetxt(ppDataFolderName+'/'+'energy.csv', energy, fmt='%.2f')
        # np.savetxt(ppDataFolderName+'/'+'energyCumAvg.csv', eTemporalCumAvg, fmt='%.2f')
        
        # x_plot = time[1:]
        # y_plot = eTemporalCumAvg[1:]
        # plt.plot(x_plot, y_plot, linewidth=0.6, label='cumulative avg')
        plt.legend()
        plt.xscale("log")
        plt.grid(linewidth=0.3)
        
        if saveSwitch:
            plt.savefig(ppPlotFolderName+'/Kurt.png',dpi=dpiRes)
    ############################ KURTOSIS PLOT ##########################

    ############################ totCom PLOT ##########################
    
    # time = data.t
    # energy = data.e
    
    time = np.loadtxt(ppDataFolderName+'/time.csv',  delimiter=',', dtype=int)
    totXComVec = np.loadtxt(ppDataFolderName+'/xComTotVec.csv',  delimiter=',', dtype=float)
    totYComVec = np.loadtxt(ppDataFolderName+'/yComTotVec.csv',  delimiter=',', dtype=float)
    totRComVec = np.loadtxt(ppDataFolderName+'/RComTotVec.csv',  delimiter=',', dtype=float)
    
    # np.savetxt(ppDataFolderName+'/'+'time.csv', time, fmt='%d')
    if plot_switches_dict['temporals']:
        
        x_plot = time[0:]
        y1_plot = totXComVec[0:]
        y2_plot = totYComVec[0:]
        y3_plot = totRComVec[0:]
        
        plt.figure()
        plt.plot(x_plot, y1_plot, linewidth=0.6, label='totXCom')
        plt.plot(x_plot, y2_plot, linewidth=0.6, label='totYCom')
        plt.axvline(x=WaitingMCSteps, linewidth=1, linestyle='dashed', c ='black')
        plt.xlabel('MCS')
        plt.ylabel(r'$\bar{x}, \bar{y}$')
        plt.title(r'$\bar{x}, \bar{y}$')
        plt.legend()
        plt.xscale("log")
        plt.grid(linewidth=0.3)
        
        if saveSwitch:
            plt.savefig(ppPlotFolderName+'/xBaryBar.png',dpi=dpiRes)
        
    
        plt.figure()
        plt.plot(x_plot, y3_plot, linewidth=0.6, label='totRCom')
        plt.axvline(x=WaitingMCSteps, linewidth=1, linestyle='dashed', c ='black')
        plt.xlabel('MCS')
        plt.ylabel(r'$\bar{R}$')
        plt.title(r'$\bar{R}$')
        plt.legend()
        plt.xscale("log")
        plt.grid(linewidth=0.3)
        
        if saveSwitch:
            plt.savefig(ppPlotFolderName+'/rBar.png',dpi=dpiRes)
    ############################ totCom PLOT ##########################

    
    # ############################ ENERGY AUTO-CORRELATION PLOT ##########
    # # max_lag = int(np.floor(0.1*len(eTemporal)))
    # # max_lag = 2*WaitingMCSteps
    # try:
    #     energyAutoCorrData = np.loadtxt(ppDataFolderName+'/'+'energyAutocorr.csv', delimiter=',')
    #     x_plot = energyAutoCorrData[0,:]
    #     y_plot = energyAutoCorrData[1,:]
        
    # except:
            
    #     max_lag = int(np.floor(0.25*len(energy)))
    #     energyAutoCorr = np.zeros(1+max_lag)
    #     energyAutoCorr[0] = 1
    #     for lag in range(1,max_lag+1):
    #         # energySet1 = eTemporal[lag:]
    #         # energySet2 = eTemporal[0:-lag]
    #         energySet1 = energy[lag:]
    #         energySet2 = energy[0:-lag]
    #         energyAutoCorr[lag] = np.corrcoef(energySet1, energySet2)[0,1]
        
    #     x_plot = list(range(max_lag+1))
    #     y_plot = energyAutoCorr
        
    #     energyAutoCorrData = np.zeros([2,len(x_plot)])
    #     energyAutoCorrData[0,:] = x_plot
    #     energyAutoCorrData[1,:] = y_plot
    #     np.savetxt(ppDataFolderName+'/'+'energyAutocorr.csv', energyAutoCorrData, fmt='%6.2f', delimiter=',')
        
    # plt.figure()
    # plt.plot(x_plot, y_plot)
    # # plt.axvline(x=WaitingMCSteps, linewidth=1, linestyle='dashed', c ='black')
    # plt.xlabel('lag (number of samples)')
    # plt.ylabel('Energy autocorrelation')
    # plt.xscale("log")
    # plt.grid(linewidth=0.4)
    
    # if saveSwitch:
    #     plt.savefig(ppPlotFolderName+'/energyAutocorr.png',dpi=dpiRes)
    
    
    # ############################ ENERGY AUTO-CORRELATIO PLOT ##########
    
    
    # ############################ ENERGY STD PLOT ######################    
    
    
    # np.savetxt(ppDataFolderName+'/'+'samplePerWindow.csv', [numberOfSamplesPerWindow], fmt='%d', delimiter=',')
    
    # try:
    #     energyWinStdData = np.loadtxt(ppDataFolderName+'/'+'energyWinStd.csv', delimiter=',')
    #     energyStdPlottingCum = np.loadtxt(ppDataFolderName+'/'+'energyStdCum.csv', delimiter=',')
        
    #     timeStdPlotting = (energyWinStdData[0,:]+energyWinStdData[1,:])/2.
    #     x_plot = timeStdPlotting
    #     energyStdPlottingWin = energyWinStdData[2,:]
    #     y_plot = energyStdPlottingWin
    #     y_plot_cum = energyStdPlottingCum
            
    # except:
    #     # i = WaitingMCSteps
    #     i=0
    #     timeStdPlotting = []
    #     energyStdPlottingWin = []
    #     energyStdPlottingCum = []
        
    #     winBegData = []
    #     winEndData = []
    #     while( (i + numberOfSamplesPerWindow-1) <= len(time)-1 ):
                
    #         windowBegIdx = i
    #         # windowBegIdx = WaitingMCSteps
    #         windowEndIdx = i + numberOfSamplesPerWindow-1   #inclusive
    #         windowBegTime = time[windowBegIdx]
    #         winBegData.append(windowBegTime)
    #         windowEndTime = time[windowEndIdx]
    #         winEndData.append(windowEndTime)
    #         windowAvgTime = (windowBegTime + windowEndTime)/2
            
    #         timeStdPlotting.append(windowAvgTime)
    #         energyStdPlottingWin.append(np.std(energy[windowBegIdx:windowEndIdx+1]))
    #         energyStdPlottingCum.append(np.std(energy[0:windowEndIdx+1]))
            
            
    #         i+=1
        
    #     timeStdPlotting = np.array(timeStdPlotting)
    #     energyStdPlottingWin = np.array(energyStdPlottingWin)
    #     energyStdPlottingCum = np.array(energyStdPlottingCum)
        
    #     winBegData = np.array(winBegData)
    #     winEndData = np.array(winEndData)
    #     energyWinStdData = np.zeros([3,len(timeStdPlotting)])
    #     energyWinStdData[0,:] = winBegData
    #     energyWinStdData[1,:] = winEndData
    #     energyWinStdData[2,:] = energyStdPlottingWin
            
    #     np.savetxt(ppDataFolderName+'/'+'energyWinStd.csv', energyWinStdData, fmt='%10.2f', delimiter=',')
    #     np.savetxt(ppDataFolderName+'/'+'energyStdCum.csv', energyStdPlottingCum, fmt='%.2f', delimiter=',')
        
    #     x_plot = timeStdPlotting
    #     y_plot = energyStdPlottingWin
    #     y_plot_cum = energyStdPlottingCum
    
    # plt.figure()
    # plt.plot(x_plot, y_plot)
    
    # energyStdPlottingWinCumAvg =cumulativeAveraging(timeStdPlotting, energyStdPlottingWin)
    # y_plot_cumAvg = energyStdPlottingWinCumAvg
    # plt.plot(x_plot, y_plot_cumAvg, label='cum avg')
    # np.savetxt(ppDataFolderName+'/'+'energyWinStdCumAvg.csv', energyStdPlottingWinCumAvg, fmt='%.2f', delimiter=',')
    
    # title = 'Moving window Std (window width : '+str(numberOfSamplesPerWindow)+' samples)'
    # plt.title(title)
    # # plt.axvline(x=WaitingMCSteps, linewidth=1, linestyle='dashed', c ='black')
    # plt.xlabel('MC step')
    # plt.ylabel('Energy standard deviation')
    # plt.xscale("log")
    # plt.legend()
    # plt.grid(linewidth=0.4)
    
    # if saveSwitch:
    #     plt.savefig(ppPlotFolderName+'/energyStdWin.png',dpi=dpiRes)
    
    # # x_plot = timeStdPlotting
    # # y_plot = energyStdPlottingCum
    
    # plt.figure()
    # plt.plot(x_plot, y_plot_cum)
    # title = 'Cumulative window Std'
    # plt.title(title)
    # # plt.axvline(x=WaitingMCSteps, linewidth=1, linestyle='dashed', c ='black')
    # plt.xlabel('MC step')
    # plt.ylabel('Energy standard deviation')
    # plt.xscale("log")
    # plt.grid(linewidth=0.4)
    
    # if saveSwitch:
    #     plt.savefig(ppPlotFolderName+'/energyStdCum.png',dpi=dpiRes)
    
    # ############################ ENERGY STD PLOT ######################
    
    function_time_end = TIME.time()
    delta_t_function = function_time_end-function_time_start
    print('Temporal plots function was done: dt = '+str(delta_t_function))
    
    return


def trajectoriesV2(saveSwitch):
    
    function_time_start = TIME.time()
    
    time = np.loadtxt(ppDataFolderName+'/time.csv',  delimiter=',', dtype=int)
    
    dpiRes = 500
    
    ind=0
    t = time[ind]
    while(t!=t_w):
        t = time[ind]
        ind+=1
    ind-=1
    t_w_ind = ind
    
    ############################ TRAJECTORIES PLOT ##########################
    
    # snapshotsX = np.loadtxt('snapshotsX.csv',  delimiter=',')
    # snapshotsY = np.loadtxt('snapshotsY.csv',  delimiter=',')
    # xCom = np.loadtxt(ppDataFolderName+'/'+'xCom.csv', delimiter=',')
    # yCom = np.loadtxt(ppDataFolderName+'/'+'yCom.csv', delimiter=',')
    
    # xCom = np.loadtxt(ppDataFolderName+'/'+'xComInit.csv', delimiter=',')
    # yCom = np.loadtxt(ppDataFolderName+'/'+'yComInit.csv', delimiter=',')
    
    # xCom = xCom.reshape([NumCells+1,1])
    # yCom = yCom.reshape([NumCells+1,1])
    
    # bunchC=0
    # while(1):
    #     try:
    #         xCom_concat = np.loadtxt(ppDataFolderName+'/'+'xComBunch_'+str(bunchC)+'.csv', delimiter=',')
    #         yCom_concat = np.loadtxt(ppDataFolderName+'/'+'yComBunch_'+str(bunchC)+'.csv', delimiter=',')
            
    #         xCom = np.concatenate([xCom,xCom_concat], 1)
    #         yCom = np.concatenate([yCom,yCom_concat], 1)
    #         bunchC+=1
    #     except:
    #         break

    plt.figure()
    for cellIndex in list(range(1, NumCells+1)):
        
        x_plot = xCom[cellIndex, :]
        y_plot = yCom[cellIndex, :]
        
        plt.plot(y_plot, x_plot, linewidth=0.2)
    
    # plt.plot( [-2,0,10],[1,4,7], linewidth=0.2)
        
    plt.xlabel('Y')
    plt.gca().invert_yaxis()
    plt.ylabel('X')
    plt.grid(linewidth=0.2)
    # plt.axis('equal', xlim=(np.min(snapshotsY)-10, np.max(snapshotsY)+10), ylim=(np.min(snapshotsX)-10, np.max(snapshotsX)+10))
    plt.axis('equal')
    plt.title('unfolded trajectories')
    # plt.plot([-0.5,L+0.5],[-0.5,-0.5],linestyle='dashed', linewidth=0.6, color='black')
    # plt.plot([-0.5,-0.5],[-0.5,L+0.5],linestyle='dashed', linewidth=0.6, color='black')
    # plt.plot([L+0.5,-0.5],[L+0.5,L+0.5],linestyle='dashed', linewidth=0.6, color='black')
    # plt.plot([L+0.5,L+0.5],[L+0.5,-0.5],linestyle='dashed', linewidth=0.6, color='black')
    plt.plot([0,Ly],[0,0],linestyle='dashed', linewidth=0.6, color='black')
    plt.plot([0,0],[0,Lx],linestyle='dashed', linewidth=0.6, color='black')
    plt.plot([Ly,0],[Lx,Lx],linestyle='dashed', linewidth=0.6, color='black')
    plt.plot([Ly,Ly],[Lx,0],linestyle='dashed', linewidth=0.6, color='black')
    
    if saveSwitch:
        plt.savefig(ppPlotFolderName+'/'+'trajectories.png',dpi=dpiRes)
    ############################ TRAJECTORIES PLOT ##########################
    
    ############################ TRAJECTORIES PLOT T_W ##########################
    plt.figure()
    for cellIndex in list(range(1, NumCells+1)):
        
        x_plot = xCom[cellIndex, t_w_ind:]
        y_plot = yCom[cellIndex, t_w_ind:]
        
        plt.plot(y_plot, x_plot, linewidth=0.2)
    
    # plt.plot( [-2,0,10],[1,4,7], linewidth=0.2)
        
    plt.xlabel('Y')
    plt.gca().invert_yaxis()
    plt.ylabel('X')
    plt.grid(linewidth=0.2)
    # plt.axis('equal', xlim=(np.min(snapshotsY)-10, np.max(snapshotsY)+10), ylim=(np.min(snapshotsX)-10, np.max(snapshotsX)+10))
    plt.axis('equal')
    plt.title('unfolded trajectories')
    # plt.plot([-0.5,L+0.5],[-0.5,-0.5],linestyle='dashed', linewidth=0.6, color='black')
    # plt.plot([-0.5,-0.5],[-0.5,L+0.5],linestyle='dashed', linewidth=0.6, color='black')
    # plt.plot([L+0.5,-0.5],[L+0.5,L+0.5],linestyle='dashed', linewidth=0.6, color='black')
    # plt.plot([L+0.5,L+0.5],[L+0.5,-0.5],linestyle='dashed', linewidth=0.6, color='black')
    plt.plot([0,Ly],[0,0],linestyle='dashed', linewidth=0.6, color='black')
    plt.plot([0,0],[0,Lx],linestyle='dashed', linewidth=0.6, color='black')
    plt.plot([Ly,0],[Lx,Lx],linestyle='dashed', linewidth=0.6, color='black')
    plt.plot([Ly,Ly],[Lx,0],linestyle='dashed', linewidth=0.6, color='black')
    
    if saveSwitch:
        plt.savefig(ppPlotFolderName+'/'+'trajectories_t_w.png',dpi=dpiRes)
    ############################ TRAJECTORIES PLOT T_W ##########################

        
    
    
    
    # ############################ AVGR2 PLOT ##########################    
    # # x_w = xCom[:, t_w_ind]
    # # y_w = yCom[:, t_w_ind]
    
    # # deltaX = 0.0*xCom
    # # deltaY = 0.0*yCom
    # # deltaR2 = 0.0*xCom
    # deltaT = 0*(time)
    
    # for snapshotC in range(len(time)):
        
    #     t = time[snapshotC]
        
    #     if (t>=t_w):
            
    #         deltaT[snapshotC]=t-t_w
            
    #         # deltaX[:,snapshotC] = xCom[:,snapshotC]-x_w[:]
    #         # deltaY[:,snapshotC] = yCom[:,snapshotC]-y_w[:]
            
    #         # deltaR2[:,snapshotC] = (deltaX[:,snapshotC])**2 + (deltaY[:,snapshotC])**2
    
    # msdAvgStd = np.loadtxt(ppDataFolderName+'/msdAvgStd.csv',  delimiter=',')
    # # deltaR2Avg = np.mean(deltaR2[1:,:], axis = 0)
    # # deltaR2Std = np.std( deltaR2[1:,:], axis = 0)
    # deltaR2Avg = msdAvgStd[0,:]
    # deltaR2Std = msdAvgStd[1,:]
    
    # # np.savetxt(ppDataFolderName+'/'+'deltaT.csv', deltaT, fmt='%9d', delimiter=',')
    # # np.savetxt(ppDataFolderName+'/'+'deltaX.csv', deltaX, fmt='%9.2f', delimiter=',')
    # # np.savetxt(ppDataFolderName+'/'+'deltaY.csv', deltaY, fmt='%9.2f', delimiter=',')
    # # np.savetxt(ppDataFolderName+'/'+'deltaR2.csv', deltaR2, fmt='%9.2f', delimiter=',')
    # # np.savetxt(ppDataFolderName+'/'+'deltaR2Avg.csv', deltaR2Avg, fmt='%9.2f', delimiter=',')
    # # np.savetxt(ppDataFolderName+'/'+'deltaR2Std.csv', deltaR2Std, fmt='%9.2f', delimiter=',')
    
    
    # plt.figure()
    # plt.xlabel(r'$t-t_w$'+' (MCS)')
    # plt.ylabel(r'$\langle r^2 \rangle$')
    # plt.grid()

    # x_plot = deltaT[t_w_ind+1:]
    # y_plot = deltaR2Avg[t_w_ind+1:]
    # y_plot_err = deltaR2Std[t_w_ind+1:]/np.sqrt(NumCells)
    # plt.errorbar(x_plot, y_plot, yerr=y_plot_err, fmt="o", markersize =0.8, elinewidth=0.2, label='MSD', zorder=1)
    
    # b = deltaR2Avg[t_w_ind+1]/deltaT[t_w_ind+1]
    # y_plot_supposed = b*x_plot
    # plt.plot(x_plot, y_plot_supposed, zorder=2, color='k', linestyle='dashed')
    
    # # plt.scatter(x_plot, y_plot, label='MSD', c='black', s=0.4)
    # plt.xscale("log")
    # plt.yscale("log")
    # plt.legend()
    
    # if saveSwitch:
    #     plt.savefig(ppPlotFolderName+'/'+'avgR2.png',dpi=dpiRes)
        
    # ############################ AVGR2 PLOT ##########################
    
    
    ############################ newMSD PLOT ##########################    
    
    x_plot = []
    y_plot = []
    y_plot_err = []
    
    with open(ppDataFolderName+'/newMSD.csv', 'r') as newMSD:
        # Iterate through each line in the file
        for line in newMSD:
            # Strip any leading/trailing whitespace (including newline characters)
            line = line.strip()
            # Split the line using ':' as the delimiter
            parts = line.split(':')
            dt = int(parts[0])
            num_list = parts[1].split(',')
            num_list = [float(num) for num in num_list]
            datapoints = np.array(num_list)
            
            x_plot.append(dt)
            y_plot.append(np.mean(datapoints))
            y_plot_err.append(np.std(datapoints)/np.sqrt(len(datapoints)))
            #y_plot_err.append(np.std(datapoints))
            
    
    
    plt.figure()
    plt.xlabel(r'$t-t_w$'+' (MCS)')
    plt.ylabel(r'$\langle r^2 \rangle$')
    plt.grid()

    x_plot = np.array(x_plot)
    y_plot = np.array(y_plot)
    y_plot_err = np.array(y_plot_err)
    plt.errorbar(x_plot, y_plot, yerr=y_plot_err, fmt="o", markersize =0.8, elinewidth=0.4, label='MSD', zorder=1, ecolor='c')
    
    # b = deltaR2Avg[t_w_ind+1]/deltaT[t_w_ind+1]
    b = y_plot[0]/x_plot[0]
    y_plot_supposed = b*x_plot
    plt.plot(x_plot, y_plot_supposed, zorder=2, color='k', linestyle='dashed')
    
    # plt.scatter(x_plot, y_plot, label='MSD', c='black', s=0.4)
    plt.xscale("log")
    plt.yscale("log")
    plt.legend()
    
    if saveSwitch:
        plt.savefig(ppPlotFolderName+'/'+'avgR2.png',dpi=dpiRes)
    
    ############################ newMSD PLOT ##########################    
    
    
    ############################ AVGR2_T_0 PLOT ##########################    
    
    xComInit = np.loadtxt(ppDataFolderName+"/xComInit.csv", delimiter=',')
    yComInit = np.loadtxt(ppDataFolderName+"/yComInit.csv", delimiter=',')
    
    x_plot= time
    y_plot = 0.0*time
    y_plot_err  = 0.0*time
    for t_c in range(len(time)):
        # r2_vec = ( (xCom[1:, t_c]) - (xComInit[1:]) )**2 + ( (yCom[1:, t_c]) - (yComInit[1:]) )**2
        r2_vec = ( (xCom[1:, t_c] - np.mean(xCom[1:, t_c]) ) - (xComInit[1:] - np.mean(xComInit[1:])) )**2 + ( (yCom[1:, t_c] - np.mean(yCom[1:, t_c])) - (yComInit[1:] - np.mean(yComInit[1:])) )**2
        y_plot[t_c] = np.mean(r2_vec)
        y_plot_err[t_c] = np.std(r2_vec)/np.sqrt(NumCells - 1)
    
    msdAvgStd = np.loadtxt(ppDataFolderName+'/msdAvgStd.csv',  delimiter=',')

    plt.figure()
    plt.xlabel(r'$t$'+' (MCS)')
    plt.ylabel(r'$\langle r^2 \rangle$')
    plt.grid()

    plt.errorbar(x_plot, y_plot, yerr=y_plot_err, fmt="o", markersize =0.8, elinewidth=0.2, label='MSD', zorder=1)
    
    b = y_plot[1]/time[1]
    y_plot_supposed = b*x_plot
    plt.plot(x_plot, y_plot_supposed, zorder=2, color='k', linestyle='dashed')
    
    # plt.scatter(x_plot, y_plot, label='MSD', c='black', s=0.4)
    plt.xscale("log")
    plt.yscale("log")
    plt.legend()
    
    if saveSwitch:
        plt.savefig(ppPlotFolderName+'/'+'avgR2_T_0.png',dpi=dpiRes)
        
############################ AVGR2_T_0 PLOT ##########################    
    
    function_time_end = TIME.time()
    delta_t_function = function_time_end-function_time_start
    print('Trajectories function was done: dt = '+str(delta_t_function))
    
    return

def qFunc(saveSwitch):
    
    time = np.loadtxt(ppDataFolderName+'/time.csv',  delimiter=',', dtype=int)
    
    function_time_start = TIME.time()
    
    ind=0
    t = time[ind]
    while(t!=t_w):
        t = time[ind]
        ind+=1
    ind-=1
    t_w_ind = ind
    
    dpiRes=500
    
    # q =    np.loadtxt(ppDataFolderName+'/'+'isoperi.csv', delimiter=',', dtype=float)
    qAvgStd = np.loadtxt(ppDataFolderName+'/'+'isoperiAvgStd.csv', delimiter=',', dtype=float)
    qAvg = qAvgStd[0,:]
    qStd = qAvgStd[1,:]
    qAvgCumAvg = cumulativeAveraging(time, qAvg)
    np.savetxt(ppDataFolderName+'/'+'isoperiAvgCumAvg.csv', qAvgCumAvg, delimiter=',')
    
    
    qAvg_final_avg = np.mean(qAvg[-final_sampling_numbers:])
    qAvg_final_std = np.std( qAvg[-final_sampling_numbers:])
    
    np.savetxt(ppDataFolderName+'/'+'isoperiAvgFinal_AvgStdNum.csv', [qAvg_final_avg, qAvg_final_std, final_sampling_numbers], fmt='%10.4f', delimiter=',')
    
    if plot_switches_dict['q']:
        
        plt.figure()
        plt.xlabel(r'$t$'+' (MCS)')
        plt.ylabel(r'$\langle Q \rangle$')
        plt.title(r'$\langle Q \rangle=2\sqrt{\pi}\langle { \sqrt{A} / P} \rangle$'+'\n'+r'$\langle Q \rangle_{final}=$'\
                  +str(np.floor((10**n_digits)*qAvg_final_avg)/(10**n_digits))+r'$\pm$'\
                  +str(np.floor((10**n_digits)*qAvg_final_std/np.sqrt(final_sampling_numbers))/(10**n_digits))    )
        plt.grid()
    
        x_plot = time[1:]
        y_plot = qAvg[1:]
        y_plot_err = qStd[1:]/np.sqrt(NumCells-1)
        plt.errorbar(x_plot, y_plot, yerr=y_plot_err, fmt="o", markersize =0.8, elinewidth=0.2, label=r'$\langle Q \rangle$', zorder=1)
        plt.plot(x_plot, qAvgCumAvg[1:], label='cum avg', zorder=2)
        # plt.scatter(x_plot, y_plot, label='MSD', c='black', s=0.4)
        plt.xscale("log")
        # plt.yscale("log")
        plt.legend()
        plt.axvline(x=WaitingMCSteps, linewidth=1, linestyle='dashed', c ='black')
    
        if saveSwitch:
            plt.savefig(ppPlotFolderName+'/'+'isoperi.png',dpi=dpiRes)
            
        
        qAvg = qAvgStd[0,t_w_ind:]
        qStd = qAvgStd[1,t_w_ind:]
        qAvgCumAvg = cumulativeAveraging(time[t_w_ind:], qAvg)
        
        plt.figure()
        plt.xlabel(r'$t$'+' (MCS)')
        plt.ylabel(r'$\langle Q \rangle$')
        plt.title(r'$\langle Q \rangle=2\sqrt{\pi}\langle { \sqrt{A} / P} \rangle$'+'\n'+r'$\langle Q \rangle_{final}=$'\
                  +str(np.floor((10**n_digits)*qAvg_final_avg)/(10**n_digits))+r'$\pm$'\
                  +str(np.floor((10**n_digits)*qAvg_final_std/np.sqrt(final_sampling_numbers))/(10**n_digits))    )
        plt.grid()
    
        x_plot = time[t_w_ind:]
        y_plot = qAvg
        y_plot_err = qStd/np.sqrt(NumCells-1)
        plt.errorbar(x_plot, y_plot, yerr=y_plot_err, fmt="o", markersize =0.8, elinewidth=0.2, label=r'$\langle Q \rangle$', zorder=1)
        plt.plot(x_plot, qAvgCumAvg, label='cum avg', zorder=2)
        # plt.scatter(x_plot, y_plot, label='MSD', c='black', s=0.4)
        plt.xscale("log")
        # plt.yscale("log")
        plt.legend()
        plt.axvline(x=WaitingMCSteps, linewidth=1, linestyle='dashed', c ='black')
    
        if saveSwitch:
            plt.savefig(ppPlotFolderName+'/'+'isoperi_eq.png',dpi=dpiRes)
        
    
    ################ EQ DISTRIBUTION ################
    distributionFunc("q")
    distributionFunc("q_inv")
    ################ EQ DISTRIBUTION ################
    
    function_time_end = TIME.time()
    delta_t_function = function_time_end-function_time_start
    
    print('q function was done: dt = '+str(delta_t_function))
    
    return
    
def ARFunc(saveSwitch):
    
    time = np.loadtxt(ppDataFolderName+'/time.csv',  delimiter=',', dtype=int)
    
    function_time_start = TIME.time()
    
    ind=0
    t = time[ind]
    while(t!=t_w):
        t = time[ind]
        ind+=1
    ind-=1
    t_w_ind = ind
    
    dpiRes=500
    
    # q =    np.loadtxt(ppDataFolderName+'/'+'isoperi.csv', delimiter=',', dtype=float)
    qAvgStd = np.loadtxt(ppDataFolderName+'/'+'ARAvgStd.csv', delimiter=',', dtype=float)
    qAvg = qAvgStd[0,:]
    qStd = qAvgStd[1,:]
    qAvgCumAvg = cumulativeAveraging(time, qAvg)
    np.savetxt(ppDataFolderName+'/'+'ARAvgCumAvg.csv', qAvgCumAvg, delimiter=',')
    
    
    qAvg_final_avg = np.mean(qAvg[-final_sampling_numbers:])
    qAvg_final_std = np.std( qAvg[-final_sampling_numbers:])
    
    np.savetxt(ppDataFolderName+'/'+'ARAvgFinal_AvgStdNum.csv', [qAvg_final_avg, qAvg_final_std, final_sampling_numbers], fmt='%10.4f', delimiter=',')
    
    
    if plot_switches_dict['AR']:
        
        plt.figure()
        plt.xlabel(r'$t$'+' (MCS)')
        plt.ylabel(r'$\langle AR \rangle$')
        plt.title(r'$\langle AR \rangle=\langle { \sqrt{I_{max}/I_{min}} } \rangle$'+'\n'+r'$\langle AR \rangle_{final}=$'\
                  +str(np.floor((10**n_digits)*qAvg_final_avg)/(10**n_digits))+r'$\pm$'\
                  +str(np.floor((10**n_digits)*qAvg_final_std/np.sqrt(final_sampling_numbers))/(10**n_digits))    )
        plt.grid()
    
        x_plot = time
        y_plot = qAvg
        y_plot_err = qStd/np.sqrt(NumCells-1)
        plt.errorbar(x_plot, y_plot, yerr=y_plot_err, fmt="o", markersize =0.8, elinewidth=0.2, label=r'$\langle AR \rangle$', zorder=1)
        plt.plot(x_plot, qAvgCumAvg, label='cum avg', zorder=2)
        plt.axvline(x=WaitingMCSteps, linewidth=1, linestyle='dashed', c ='black')
        # plt.scatter(x_plot, y_plot, label='MSD', c='black', s=0.4)
        plt.xscale("log")
        # plt.yscale("log")
        plt.legend()
        
        
        if saveSwitch:
            plt.savefig(ppPlotFolderName+'/'+'AR.png',dpi=dpiRes)
            
            # plt.xlim([1e6,MaxMCSteps])
            # plt.ylim([1.285,1.325])
            
            # plt.xscale("linear")
            # plt.savefig(ppPlotFolderName+'/'+'AR_zoom.png',dpi=dpiRes)
        
            
        
    
        qAvg = qAvgStd[0,t_w_ind:]
        qStd = qAvgStd[1,t_w_ind:]
        qAvgCumAvg = cumulativeAveraging(time[t_w_ind:], qAvg)
        
        plt.figure()
        plt.xlabel(r'$t$'+' (MCS)')
        plt.ylabel(r'$\langle AR \rangle$')
        plt.title(r'$\langle AR \rangle=\langle { \sqrt{I_{max}/I_{min}} } \rangle$'+'\n'+r'$\langle AR \rangle_{final}=$'\
                  +str(np.floor((10**n_digits)*qAvg_final_avg)/(10**n_digits))+r'$\pm$'\
                  +str(np.floor((10**n_digits)*qAvg_final_std/np.sqrt(final_sampling_numbers))/(10**n_digits))    )
        plt.grid()
    
        x_plot = time[t_w_ind:]
        y_plot = qAvg
        y_plot_err = qStd/np.sqrt(NumCells-1)
        plt.errorbar(x_plot, y_plot, yerr=y_plot_err, fmt="o", markersize =0.8, elinewidth=0.2, label=r'$\langle AR \rangle$', zorder=1)
        plt.plot(x_plot, qAvgCumAvg, label='cum avg', zorder=2)
        plt.axvline(x=WaitingMCSteps, linewidth=1, linestyle='dashed', c ='black')
        # plt.scatter(x_plot, y_plot, label='MSD', c='black', s=0.4)
        plt.xscale("log")
        # plt.yscale("log")
        plt.legend()
        
        if saveSwitch:
            plt.savefig(ppPlotFolderName+'/'+'AR_eq.png',dpi=dpiRes)
        
    
    ################ EQ DISTRIBUTION ################
    distributionFunc("AR")
    ################ EQ DISTRIBUTION ################
    
    function_time_end = TIME.time()
    delta_t_function = function_time_end-function_time_start
    
    print('AR function was done: dt = '+str(delta_t_function))
    
    return

def circFunc(saveSwitch):
    
    time = np.loadtxt(ppDataFolderName+'/time.csv',  delimiter=',', dtype=int)
    
    function_time_start = TIME.time()
    
    ind=0
    t = time[ind]
    while(t!=t_w):
        t = time[ind]
        ind+=1
    ind-=1
    t_w_ind = ind
    
    dpiRes=500
    
    # q =    np.loadtxt(ppDataFolderName+'/'+'isoperi.csv', delimiter=',', dtype=float)
    qAvgStd = np.loadtxt(ppDataFolderName+'/'+'circAvgStd.csv', delimiter=',', dtype=float)
    qAvg = qAvgStd[0,:]
    qStd = qAvgStd[1,:]
    qAvgCumAvg = cumulativeAveraging(time, qAvg)
    np.savetxt(ppDataFolderName+'/'+'circAvgCumAvg.csv', qAvgCumAvg, delimiter=',')
    
    
    # qAvg_final_avg = np.mean(qAvg[-final_sampling_numbers:])
    # qAvg_final_std = np.std( qAvg[-final_sampling_numbers:])
    
    # np.savetxt(ppDataFolderName+'/'+'ARAvgFinal_AvgStdNum.csv', [qAvg_final_avg, qAvg_final_std, final_sampling_numbers], fmt='%10.4f', delimiter=',')
    
    
    if plot_switches_dict['circ']:
        
        plt.figure()
        plt.xlabel(r'$t$'+' (MCS)')
        plt.ylabel(r'$\langle circularity \rangle$')
        plt.title(r'$\langle circularity \rangle=\langle { {A^2}/{(2\pi (I_{xx}+I_{yy})})} \rangle$')
        plt.grid()
    
        x_plot = time
        y_plot = qAvg
        y_plot_err = qStd/np.sqrt(NumCells-1)
        plt.errorbar(x_plot, y_plot, yerr=y_plot_err, fmt="o", markersize =0.8, elinewidth=0.2, label=r'$\langle circ \rangle$', zorder=1)
        plt.plot(x_plot, qAvgCumAvg, label='cum avg', zorder=2)
        plt.axvline(x=WaitingMCSteps, linewidth=1, linestyle='dashed', c ='black')
        # plt.scatter(x_plot, y_plot, label='MSD', c='black', s=0.4)
        plt.xscale("log")
        # plt.yscale("log")
        plt.legend()
        
        
        if saveSwitch:
            plt.savefig(ppPlotFolderName+'/'+'circ.png',dpi=dpiRes)
            
            # plt.xlim([1e6,MaxMCSteps])
            # plt.ylim([1.285,1.325])
            
            # plt.xscale("linear")
            # plt.savefig(ppPlotFolderName+'/'+'AR_zoom.png',dpi=dpiRes)
        
            
        
    
        qAvg = qAvgStd[0,t_w_ind:]
        qStd = qAvgStd[1,t_w_ind:]
        qAvgCumAvg = cumulativeAveraging(time[t_w_ind:], qAvg)
        
        plt.figure()
        plt.xlabel(r'$t$'+' (MCS)')
        plt.ylabel(r'$\langle circularity \rangle$')
        plt.title(r'$\langle circularity \rangle=\langle { {A^2}/{(2\pi (I_{xx}+I_{yy})})} \rangle$')
        plt.grid()
    
        x_plot = time[t_w_ind:]
        y_plot = qAvg
        y_plot_err = qStd/np.sqrt(NumCells-1)
        plt.errorbar(x_plot, y_plot, yerr=y_plot_err, fmt="o", markersize =0.8, elinewidth=0.2, label=r'$\langle circ \rangle$', zorder=1)
        plt.plot(x_plot, qAvgCumAvg, label='cum avg', zorder=2)
        plt.axvline(x=WaitingMCSteps, linewidth=1, linestyle='dashed', c ='black')
        # plt.scatter(x_plot, y_plot, label='MSD', c='black', s=0.4)
        plt.xscale("log")
        # plt.yscale("log")
        plt.legend()
        
        if saveSwitch:
            plt.savefig(ppPlotFolderName+'/'+'circ_eq.png',dpi=dpiRes)
        
    
    ################ EQ DISTRIBUTION ################
    distributionFunc("circ")
    ################ EQ DISTRIBUTION ################
    
    function_time_end = TIME.time()
    delta_t_function = function_time_end-function_time_start
    
    print('circ function was done: dt = '+str(delta_t_function))
    
    return

def orderParametersFunc(saveSwitch):
    
    function_time_start = TIME.time()
    
    freudcmap = 'plasma_r'
        

    ngridx=1000
    ngridy=1000
    xi = np.linspace(-0.5, L-0.5, ngridx)
    yi = np.linspace(-0.5, L-0.5, ngridy)
    Xi, Yi = np.meshgrid(xi, yi)
    
    dpi_res=500
    
    # tTemporalSnapshots = np.loadtxt('tTemporalSnapshots.csv',  delimiter=',', dtype=int)
    # snapshotsX = np.loadtxt('snapshotsX.csv',  delimiter=',')
    # snapshotsY = np.loadtxt('snapshotsY.csv',  delimiter=',')
    
    tTemporalSnapshots = np.loadtxt(ppDataFolderName+"/"+'time.csv',  delimiter=',', dtype=int)
    # snapshotsX = np.loadtxt(ppDataFolderName+"/"+'xCom.csv',  delimiter=',')
    # snapshotsY = np.loadtxt(ppDataFolderName+"/"+'yCom.csv',  delimiter=',')
    # snapshotsX = xCom
    # snapshotsY = yCom
    
    # Wrapper (for periodicity)
    periodicSnapshotsX = ((xCom+0.5)%L)-0.5
    periodicSnapshotsY = ((yCom+0.5)%L)-0.5
    
    # Coordinates origin changing
    periodicSnapshotsX = periodicSnapshotsX - (L/2-0.5)
    periodicSnapshotsY = periodicSnapshotsY - (L/2-0.5)
    
    periodicSnapshotsX[0,:] = 0.*periodicSnapshotsX[0,:]
    periodicSnapshotsY[0,:] = 0.*periodicSnapshotsY[0,:]
    
    snapshotsHexAbs = 0.0*xCom
    snapshotsHexAbsAvg = np.zeros(np.shape(xCom)[1])
    snapshotsHexAbsStd = np.zeros(np.shape(xCom)[1])
    snapshotsHexAng = 0.0*xCom
    
    box = freud.box.Box(Lx=L, Ly=L, Lz=0, xy=0, xz=0, yz=0, is2D=True)
    
    eqSampleC=0
    hexAllCells = np.zeros([NumCells, len(eqSamplingTimes.timesList)])
    
    for t_c in range(0, len(tTemporalSnapshots)):
            
        t = tTemporalSnapshots[t_c]
        
        points = []
        
        for cellIndex in range(1,NumCells+1):
            
            point = []
            point.append(periodicSnapshotsX[cellIndex, t_c])
            point.append(periodicSnapshotsY[cellIndex, t_c])
            point.append(0.)
            
            points.append(point)
        
        
        
        hex_order = freud.order.Hexatic(k=6)
        hex_order.compute(system=(box, points))
        snapshotsHexAbs[1:,t_c]=np.abs(hex_order.particle_order)
        snapshotsHexAbsAvg[t_c]=np.mean(snapshotsHexAbs[1:,t_c])
        snapshotsHexAbsStd[t_c]=np.std( snapshotsHexAbs[1:,t_c])
        snapshotsHexAng[1:,t_c]=np.angle(hex_order.particle_order, deg=False)
        
        # print(t_c)
        
        if t in eqSamplingTimes.timesList:
            hexAllCells[:,eqSampleC] = np.abs(hex_order.particle_order)
            eqSampleC += 1
        
        if t_c == len(tTemporalSnapshots)-1:
        # if t_c == 0:
            
            if plot_switches_dict['order']:
            
                x_plot = ((xCom[1:,t_c]+0.5)%L)-0.5
                y_plot = ((yCom[1:,t_c]+0.5)%L)-0.5
                z_plot = snapshotsHexAbs[1:,t_c]
                
                # x_plot = [12, 4]
                # y_plot = [8, -1]
                # z_plot = [0.5, 0.1]
                
                plt.figure()
                plt.scatter(y_plot, x_plot, c=z_plot, cmap=freudcmap, vmin=0., vmax=1., s=10)
                plt.xlabel('Y')
                plt.gca().invert_yaxis()
                plt.ylabel('X')
                plt.grid(linewidth=0.2)
                plt.axis('equal')
                title = 'Hexatic order; '+r'$\alpha = $'+str(Alpha)+'\n' \
                +r'$t = $'+str(t)+' , '\
                +r'$\langle \psi_6 \rangle=$'+str(np.floor((10**n_digits)*snapshotsHexAbsAvg[t_c])/(10**n_digits))
                plt.title(title)
                plt.colorbar()
                plt.plot([-0.5,L-0.5],[-0.5,-0.5],linestyle='dashed', linewidth=0.6, color='black')
                plt.plot([-0.5,-0.5],[-0.5,L-0.5],linestyle='dashed', linewidth=0.6, color='black')
                plt.plot([L-0.5,-0.5],[L-0.5,L-0.5],linestyle='dashed', linewidth=0.6, color='black')
                plt.plot([L-0.5,L-0.5],[L-0.5,-0.5],linestyle='dashed', linewidth=0.6, color='black')
                if saveSwitch:
                    plt.savefig(ppPlotFolderName+"/"+'finalScatterHex.PNG', dpi=dpi_res)
                # plt.show()

            
            
                plt.figure()
                triang = tri.Triangulation(((xCom[1:,t_c]+0.5)%L)-0.5, ((yCom[1:,t_c]+0.5)%L)-0.5)
                interpolator = tri.LinearTriInterpolator(triang, z_plot)
                zi = interpolator(Xi, Yi)
                
                plt.contour(Yi, Xi, zi, levels=8, linewidths=0.5, colors='k',  vmin=0., vmax=1.)
                cntr1 = plt.contourf(Yi, Xi, zi, levels=8, cmap=freudcmap)
                plt.xlabel('Y')
                plt.gca().invert_yaxis()
                plt.ylabel('X')
                plt.grid(linewidth=0.2)
                plt.axis('equal')
                title = 'Hexatic order; '+r'$\alpha = $'+str(Alpha)+'\n' \
                +r'$t = $'+str(t)+' , '\
                +r'$\langle \psi_6 \rangle=$'+str(np.floor((10**n_digits)*snapshotsHexAbsAvg[t_c])/(10**n_digits))
                plt.title(title)
                plt.colorbar()
                plt.plot([-0.5,L-0.5],[-0.5,-0.5],linestyle='dashed', linewidth=0.6, color='black')
                plt.plot([-0.5,-0.5],[-0.5,L-0.5],linestyle='dashed', linewidth=0.6, color='black')
                plt.plot([L-0.5,-0.5],[L-0.5,L-0.5],linestyle='dashed', linewidth=0.6, color='black')
                plt.plot([L-0.5,L-0.5],[L-0.5,-0.5],linestyle='dashed', linewidth=0.6, color='black')
                # plt.savefig('aaa.png', dpi=1000)
                
                # fig.colorbar(cntr1, ax=ax1)
                # ax1.plot(x, y, 'ko', ms=3)
                # ax1.set(xlim=(-2, 2), ylim=(-2, 2))
                # ax1.set_title('grid and contour (%d points, %d grid points)' %
                #   (npts, ngridx * ngridy))
                if saveSwitch:
                    plt.savefig(ppPlotFolderName+"/"+'finalContourHex.PNG', dpi=dpi_res)
            # plt.show()        
        # break
    
    del periodicSnapshotsX
    del periodicSnapshotsY
        
    # np.savetxt(ppDataFolderName+"/"+'snapshotsHexAbs.csv', X=snapshotsHexAbs, fmt='%6.3f', delimiter=',')
    np.savetxt(ppDataFolderName+"/"+'snapshotsHexAbsAvg.csv', X=snapshotsHexAbsAvg, fmt='%1.3f', delimiter=',')
    np.savetxt(ppDataFolderName+"/"+'snapshotsHexAbsStd.csv', X=snapshotsHexAbsStd, fmt='%1.3f', delimiter=',')
    # np.savetxt(ppDataFolderName+"/"+'snapshotsHexAng.csv', X=snapshotsHexAng, fmt='%6.3f', delimiter=',')
    
    hexAllCells = np.concatenate((np.zeros([1,np.shape(hexAllCells)[1]]), hexAllCells), 0)
    np.savetxt(ppDataFolderName+"/"+'hex_dist_eq.csv', X=hexAllCells, fmt='%1.4f', delimiter=',')
    
    #######################EQ HEX DISTRIBUTION#########################
    ################ EQ DISTRIBUTION ################
    distributionFunc("hex")
    ################ EQ DISTRIBUTION ################
    #######################EQ HEX DISTRIBUTION#########################
    
    
    ######################### ADVANCED PLOT ################################
    if plot_switches_dict['adv_hex']:
    
        sigmaHex = np.zeros([L,L], dtype=float)
        finalTimeIndex = len(tTemporalSnapshots)-1
        # finalTimeIndex = 0
        final_config = np.loadtxt(mainResumeFolderName+'/sigmaMatLS.csv',  delimiter=',', dtype=int)
        for row in range(L):
            for col in range(L):
                
                cellIndex = final_config[row, col]
                sigmaHex[row, col] = snapshotsHexAbs[cellIndex, finalTimeIndex]
        
        
        del snapshotsHexAbs
        del snapshotsHexAbsAvg
        del snapshotsHexAbsStd
        del snapshotsHexAng
        
        # cmap = plt.get_cmap(freudcmap, np.max(sigmaHex) - np.min(sigmaHex) + 1)
        cmap = freudcmap
        
        plt.figure()
        # Plot the image with discrete integer values and filled regions
        plt.imshow(sigmaHex, cmap=cmap, interpolation='nearest', vmin=0., vmax=1.)
        
        # Draw contour lines on top of the filled regions
        # plt.contour(sigmaMat, levels=np.unique(sigmaMat), colors='black', linewidths=0.05, antialiased=False, corner_mask=False)
        plt.colorbar()
        
        lines_list = []
        L_ver = np.shape(sigmaHex)[0]
        L_hor = np.shape(sigmaHex)[1]
        
                    
        verBorders = np.zeros([L_ver  , L_hor+1])
        for i in range(L_ver):
            for j in range(L_hor):
                if (sigmaHex[i,j]!=sigmaHex[i,(j-1)%L_hor]):
                    verBorders[i,j]=1
                    x_beg = j-0.5; y_beg=i-0.5;
                    x_end = x_beg; y_end=y_beg+1;
                    lines_list.append(np.array([[x_beg, y_beg],[x_end, y_end]]))
                    
            j = L_hor
            if (sigmaHex[i,-1]!=sigmaHex[i,0]):
                verBorders[i,j]=1
                x_beg = (L_hor-1)+0.5; y_beg=i-0.5;
                x_end = x_beg;         y_end=y_beg+1;
                lines_list.append(np.array([[x_beg, y_beg],[x_end, y_end]]))
    
    
        horBorders = np.zeros([L_ver+1, L_hor  ])
        for i in range(L_ver):
            for j in range(L_hor):
                if (sigmaHex[i,j]!=sigmaHex[(i-1)%L_ver,j]):
                    horBorders[i,j]=1
                    x_beg = j-0.5;   y_beg=i-0.5;
                    x_end = x_beg+1; y_end=y_beg;
                    lines_list.append(np.array([[x_beg, y_beg],[x_end, y_end]]))
                    
        i = L_ver
        for j in range(L_hor):
            if (sigmaHex[-1,j]!=sigmaHex[0,j]):
                horBorders[i,j]=1
                x_beg = j-0.5;   y_beg=(L_ver-1)+0.5;
                x_end = x_beg+1; y_end=y_beg;
                lines_list.append(np.array([[x_beg, y_beg],[x_end, y_end]]))
    
        # lines_list = [np.array([[20, 10],[40, 50]])]
        
        for line_c in range(len(lines_list)):
            line_beg = lines_list[line_c][0,:].copy()
            line_end = lines_list[line_c][1,:].copy()
            x_plot = [line_beg[0], line_end[0]]
            y_plot = [line_beg[1], line_end[1]]
            plt.plot(x_plot, y_plot, 'k', linewidth=0.2)
        
        plt.xlabel('Y')
        plt.ylabel('X')
        plt.title('color indicates Hexatic order')
        
        if saveSwitch:
            #file_name_pic = fileName+".png"
            plt.savefig(ppPlotFolderName+"/"+'advancedHex.PNG', dpi=dpi_res)
        
        function_time_end = TIME.time()
        delta_t_function = function_time_end-function_time_start
        print('order parameter function was done: dt = '+str(delta_t_function))
    
    return

def orderParametersFunc_Irr(saveSwitch):
    
    function_time_start = TIME.time()
    
    freudcmap = 'plasma_r'
        

    ngridx=1000
    ngridy=1000
    xi = np.linspace(0, Lx, ngridx)
    yi = np.linspace(0, Ly, ngridy)
    Xi, Yi = np.meshgrid(xi, yi)
    
    dpi_res=500
    
    # tTemporalSnapshots = np.loadtxt('tTemporalSnapshots.csv',  delimiter=',', dtype=int)
    # snapshotsX = np.loadtxt('snapshotsX.csv',  delimiter=',')
    # snapshotsY = np.loadtxt('snapshotsY.csv',  delimiter=',')
    
    tTemporalSnapshots = np.loadtxt(ppDataFolderName+"/"+'time.csv',  delimiter=',', dtype=int)
    # snapshotsX = np.loadtxt(ppDataFolderName+"/"+'xCom.csv',  delimiter=',')
    # snapshotsY = np.loadtxt(ppDataFolderName+"/"+'yCom.csv',  delimiter=',')
    # snapshotsX = xCom
    # snapshotsY = yCom
    
    # Wrapper (for periodicity)
    periodicSnapshotsX = ((xCom)%Lx)
    periodicSnapshotsY = ((yCom)%Ly)
    
    # Coordinates origin changing
    # periodicSnapshotsX = periodicSnapshotsX - (L/2-0.5)
    # periodicSnapshotsY = periodicSnapshotsY - (L/2-0.5)
    
    periodicSnapshotsX[0,:] = 0.*periodicSnapshotsX[0,:]
    periodicSnapshotsY[0,:] = 0.*periodicSnapshotsY[0,:]
    
    snapshotsHexAbs = 0.0*xCom
    snapshotsHexAbsAvg = np.zeros(np.shape(xCom)[1])
    snapshotsHexAbsStd = np.zeros(np.shape(xCom)[1])
    snapshotsHexAng = 0.0*xCom
    
    box = freud.box.Box(Lx=Lx, Ly=Ly, Lz=0, xy=0, xz=0, yz=0, is2D=True)
    
    eqSampleC=0
    hexAllCells = np.zeros([NumCells, len(eqSamplingTimes.timesList)])
    
    for t_c in range(0, len(tTemporalSnapshots)):
            
        t = tTemporalSnapshots[t_c]
        
        points = []
        
        for cellIndex in range(1,NumCells+1):
            
            point = []
            point.append(periodicSnapshotsX[cellIndex, t_c])
            point.append(periodicSnapshotsY[cellIndex, t_c])
            point.append(0.)
            
            points.append(point)
        
        
        
        hex_order = freud.order.Hexatic(k=6)
        hex_order.compute(system=(box, points))
        snapshotsHexAbs[1:,t_c]=np.abs(hex_order.particle_order)
        snapshotsHexAbsAvg[t_c]=np.mean(snapshotsHexAbs[1:,t_c])
        snapshotsHexAbsStd[t_c]=np.std( snapshotsHexAbs[1:,t_c])
        snapshotsHexAng[1:,t_c]=np.angle(hex_order.particle_order, deg=False)
        
        # print(t_c)
        
        if t in eqSamplingTimes.timesList:
            hexAllCells[:,eqSampleC] = np.abs(hex_order.particle_order)
            eqSampleC += 1
        
        
        if t_c == len(tTemporalSnapshots)-1:
        # if t_c == 0:
            
            if plot_switches_dict['order']:

                # x_plot = ((xCom[1:,t_c]+0.5)%L)-0.5
                # y_plot = ((yCom[1:,t_c]+0.5)%L)-0.5
                x_plot = ((xCom[1:,t_c])%Lx)
                y_plot = ((yCom[1:,t_c])%Ly)
                z_plot = snapshotsHexAbs[1:,t_c]
                
                # x_plot = [12, 4]
                # y_plot = [8, -1]
                # z_plot = [0.5, 0.1]
                
                plt.figure()
                plt.scatter(y_plot, x_plot, c=z_plot, cmap=freudcmap, vmin=0., vmax=1., s=10)
                plt.xlabel('Y')
                plt.gca().invert_yaxis()
                plt.ylabel('X')
                plt.grid(linewidth=0.2)
                plt.axis('equal')
                # title = 'Hexatic order; '+r'$\alpha = $'+str(Alpha)+'\n' \
                # +r'$t = $'+str(t)+' , '\
                # +r'$\langle \psi_6 \rangle=$'+str(np.floor((10**n_digits)*snapshotsHexAbsAvg[t_c])/(10**n_digits))
                title = 'Hexatic order; '+'\n' \
                +r'$t = $'+str(t)+' , '\
                +r'$\langle \psi_6 \rangle=$'+str(np.floor((10**n_digits)*snapshotsHexAbsAvg[t_c])/(10**n_digits))
                plt.title(title)
                plt.colorbar()
                plt.plot([0,Ly],[0,0],linestyle='dashed', linewidth=0.6, color='black')
                plt.plot([0,0],[0,Lx],linestyle='dashed', linewidth=0.6, color='black')
                plt.plot([Ly,0],[Lx,Lx],linestyle='dashed', linewidth=0.6, color='black')
                plt.plot([Ly,Ly],[Lx,0],linestyle='dashed', linewidth=0.6, color='black')
                if saveSwitch:
                    plt.savefig(ppPlotFolderName+"/"+'finalScatterHex.PNG', dpi=dpi_res)
                # plt.show()

            
            
                plt.figure()
                # triang = tri.Triangulation(((xCom[1:,t_c]+0.5)%L)-0.5, ((yCom[1:,t_c]+0.5)%L)-0.5)
                triang = tri.Triangulation(((xCom[1:,t_c])%Lx), ((yCom[1:,t_c])%Ly))
                interpolator = tri.LinearTriInterpolator(triang, z_plot)
                zi = interpolator(Xi, Yi)
                
                plt.contour(Yi, Xi, zi, levels=8, linewidths=0.5, colors='k',  vmin=0., vmax=1.)
                cntr1 = plt.contourf(Yi, Xi, zi, levels=8, cmap=freudcmap)
                plt.xlabel('Y')
                plt.gca().invert_yaxis()
                plt.ylabel('X')
                plt.grid(linewidth=0.2)
                plt.axis('equal')
                # title = 'Hexatic order; '+r'$\alpha = $'+str(Alpha)+'\n' \
                # +r'$t = $'+str(t)+' , '\
                # +r'$\langle \psi_6 \rangle=$'+str(np.floor((10**n_digits)*snapshotsHexAbsAvg[t_c])/(10**n_digits))
                title = 'Hexatic order; '+'\n' \
                +r'$t = $'+str(t)+' , '\
                +r'$\langle \psi_6 \rangle=$'+str(np.floor((10**n_digits)*snapshotsHexAbsAvg[t_c])/(10**n_digits))
                plt.title(title)
                plt.colorbar()
                plt.plot([0,Ly],[0,0],linestyle='dashed', linewidth=0.6, color='black')
                plt.plot([0,0],[0,Lx],linestyle='dashed', linewidth=0.6, color='black')
                plt.plot([Ly,0],[Lx,Lx],linestyle='dashed', linewidth=0.6, color='black')
                plt.plot([Ly,Ly],[Lx,0],linestyle='dashed', linewidth=0.6, color='black')
                # plt.savefig('aaa.png', dpi=1000)
            
                # fig.colorbar(cntr1, ax=ax1)
                # ax1.plot(x, y, 'ko', ms=3)
                # ax1.set(xlim=(-2, 2), ylim=(-2, 2))
                # ax1.set_title('grid and contour (%d points, %d grid points)' %
                #   (npts, ngridx * ngridy))
                if saveSwitch:
                    plt.savefig(ppPlotFolderName+"/"+'finalContourHex.PNG', dpi=dpi_res)
                # plt.show()        
            # break
    
    del periodicSnapshotsX
    del periodicSnapshotsY
    
    np.savetxt(ppDataFolderName+"/"+'hexAllCellsEq.csv', X=hexAllCells, fmt='%1.4f', delimiter=',')
    
    # np.savetxt(ppDataFolderName+"/"+'snapshotsHexAbs.csv', X=snapshotsHexAbs, fmt='%6.3f', delimiter=',')
    np.savetxt(ppDataFolderName+"/"+'snapshotsHexAbsAvg.csv', X=snapshotsHexAbsAvg, fmt='%1.3f', delimiter=',')
    np.savetxt(ppDataFolderName+"/"+'snapshotsHexAbsStd.csv', X=snapshotsHexAbsStd, fmt='%1.3f', delimiter=',')
    # np.savetxt(ppDataFolderName+"/"+'snapshotsHexAng.csv', X=snapshotsHexAng, fmt='%6.3f', delimiter=',')
    
    
    hexAllCells = np.concatenate((np.zeros([1,np.shape(hexAllCells)[1]]), hexAllCells), 0)
    np.savetxt(ppDataFolderName+"/"+'hex_dist_eq.csv', X=hexAllCells, fmt='%1.4f', delimiter=',')
    
    #######################EQ HEX DISTRIBUTION#########################
    ################ EQ DISTRIBUTION ################
    distributionFunc("hex")
    ################ EQ DISTRIBUTION ################
    #######################EQ HEX DISTRIBUTION#########################
    
    
    ######################### ADVANCED PLOT ################################
    sigmaHex = np.zeros(NSites, dtype=float)
    finalTimeIndex = len(tTemporalSnapshots)-1
    # finalTimeIndex = 0
    final_config = np.loadtxt(mainResumeFolderName+'/sigmaMatLS.csv',  delimiter=',', dtype=int)
 
    for siteC in range(NSites):
        cellIndex = final_config[siteC]
        sigmaHex[siteC] = snapshotsHexAbs[cellIndex, finalTimeIndex]
    
    
    del snapshotsHexAbs
    del snapshotsHexAbsAvg
    del snapshotsHexAbsStd
    del snapshotsHexAng
    
    # cmap = plt.get_cmap(freudcmap, np.max(sigmaHex) - np.min(sigmaHex) + 1)
    # cmap = freudcmap
    
    
    

    if plot_switches_dict['adv_hex']:
        
        cmap = plt.get_cmap(freudcmap)
        # cmap = 'freudcmap'
        plt.figure()
        cellsPlotterIrr(sigmaHex, 0*sigmaHex, ppPlotFolderName+'/advancedHex.png' , 1, cmap, 1, xCom[:,-1], yCom[:,-1])

    
    # plt.figure()
    # # Plot the image with discrete integer values and filled regions
    # plt.imshow(sigmaHex, cmap=cmap, interpolation='nearest', vmin=0., vmax=1.)
    
    # # Draw contour lines on top of the filled regions
    # # plt.contour(sigmaMat, levels=np.unique(sigmaMat), colors='black', linewidths=0.05, antialiased=False, corner_mask=False)
    # plt.colorbar()
    
    # lines_list = []
    # L_ver = np.shape(sigmaHex)[0]
    # L_hor = np.shape(sigmaHex)[1]
    
                
    # verBorders = np.zeros([L_ver  , L_hor+1])
    # for i in range(L_ver):
    #     for j in range(L_hor):
    #         if (sigmaHex[i,j]!=sigmaHex[i,(j-1)%L_hor]):
    #             verBorders[i,j]=1
    #             x_beg = j-0.5; y_beg=i-0.5;
    #             x_end = x_beg; y_end=y_beg+1;
    #             lines_list.append(np.array([[x_beg, y_beg],[x_end, y_end]]))
                
    #     j = L_hor
    #     if (sigmaHex[i,-1]!=sigmaHex[i,0]):
    #         verBorders[i,j]=1
    #         x_beg = (L_hor-1)+0.5; y_beg=i-0.5;
    #         x_end = x_beg;         y_end=y_beg+1;
    #         lines_list.append(np.array([[x_beg, y_beg],[x_end, y_end]]))


    # horBorders = np.zeros([L_ver+1, L_hor  ])
    # for i in range(L_ver):
    #     for j in range(L_hor):
    #         if (sigmaHex[i,j]!=sigmaHex[(i-1)%L_ver,j]):
    #             horBorders[i,j]=1
    #             x_beg = j-0.5;   y_beg=i-0.5;
    #             x_end = x_beg+1; y_end=y_beg;
    #             lines_list.append(np.array([[x_beg, y_beg],[x_end, y_end]]))
                
    # i = L_ver
    # for j in range(L_hor):
    #     if (sigmaHex[-1,j]!=sigmaHex[0,j]):
    #         horBorders[i,j]=1
    #         x_beg = j-0.5;   y_beg=(L_ver-1)+0.5;
    #         x_end = x_beg+1; y_end=y_beg;
    #         lines_list.append(np.array([[x_beg, y_beg],[x_end, y_end]]))

    # # lines_list = [np.array([[20, 10],[40, 50]])]
    
    # for line_c in range(len(lines_list)):
    #     line_beg = lines_list[line_c][0,:].copy()
    #     line_end = lines_list[line_c][1,:].copy()
    #     x_plot = [line_beg[0], line_end[0]]
    #     y_plot = [line_beg[1], line_end[1]]
    #     plt.plot(x_plot, y_plot, 'k', linewidth=0.2)
    
    # plt.xlabel('Y')
    # plt.ylabel('X')
    # plt.title('color indicates Hexatic order')
    
    # if saveSwitch:
    #     #file_name_pic = fileName+".png"
    #     plt.savefig(ppPlotFolderName+"/"+'advancedHex.PNG', dpi=dpi_res)
    
    # function_time_end = TIME.time()
    # delta_t_function = function_time_end-function_time_start
    # print('order parameter function was done: dt = '+str(delta_t_function))
    
    return


def temporalHex(saveSwitch):
    
    time = np.loadtxt(ppDataFolderName+'/time.csv',  delimiter=',', dtype=int)
    
    function_time_start = TIME.time()
    
    # t_w = 8e6
    
    ind=0
    t = time[ind]
    while(t!=t_w):
        t = time[ind]
        ind+=1
    ind-=1
    t_w_ind = ind
    
    dpiRes=500
    

    qAvg = np.loadtxt(ppDataFolderName+'/'+'snapshotsHexAbsAvg.csv', delimiter=',', dtype=float)
    qStd = np.loadtxt(ppDataFolderName+'/'+'snapshotsHexAbsStd.csv', delimiter=',', dtype=float)
    qAvgCumAvg = cumulativeAveraging(time, qAvg)
    np.savetxt(ppDataFolderName+'/'+'HexAvgCumAvg.csv', qAvgCumAvg, delimiter=',')
    
    
    qAvg_final_avg = np.mean(qAvg[-final_sampling_numbers:])
    qAvg_final_std = np.std( qAvg[-final_sampling_numbers:])
    
    np.savetxt(ppDataFolderName+'/'+'HexAvgFinal_AvgStdNum.csv', [qAvg_final_avg, qAvg_final_std, final_sampling_numbers], fmt='%10.4f', delimiter=',')
    
    if plot_switches_dict['temporals'] and plot_switches_dict['order']:
        
        plt.figure()
        plt.xlabel(r'$t$'+' (MCS)')
        plt.ylabel(r'$\langle |\psi_6| \rangle$')
        plt.title(r'$\langle |\psi_6| \rangle$'+'\n'+r'$\langle |\psi_6| \rangle_{final}=$'\
                  +str(np.floor((10**n_digits)*qAvg_final_avg)/(10**n_digits))+r'$\pm$'\
                  +str(np.floor((10**n_digits)*qAvg_final_std/np.sqrt(final_sampling_numbers))/(10**n_digits))    )
        plt.grid()
    
        x_plot = time
        y_plot = qAvg
        y_plot_err = qStd/np.sqrt(NumCells-1)
        plt.errorbar(x_plot, y_plot, yerr=y_plot_err, fmt="o", markersize =0.8, elinewidth=0.2, label=r'$\langle AR \rangle$', zorder=1)
        plt.plot(x_plot, qAvgCumAvg, label='cum avg', zorder=2)
        plt.axvline(x=WaitingMCSteps, linewidth=1, linestyle='dashed', c ='black')
        # plt.scatter(x_plot, y_plot, label='MSD', c='black', s=0.4)
        plt.xscale("log")
        # plt.yscale("log")
        plt.legend()
        
        
        if saveSwitch:
            plt.savefig(ppPlotFolderName+'/'+'Hex.png',dpi=dpiRes)
        
        
        qAvg = qAvg[t_w_ind:]
        qStd = qStd[t_w_ind:]
        qAvgCumAvg = cumulativeAveraging(time[t_w_ind:], qAvg)    
        
        plt.figure()
        plt.xlabel(r'$t$'+' (MCS)')
        plt.ylabel(r'$\langle |\psi_6| \rangle$')
        plt.title(r'$\langle |\psi_6| \rangle$'+'\n'+r'$\langle |\psi_6| \rangle_{final}=$'\
                  +str(np.floor((10**n_digits)*qAvg_final_avg)/(10**n_digits))+r'$\pm$'\
                  +str(np.floor((10**n_digits)*qAvg_final_std/np.sqrt(final_sampling_numbers))/(10**n_digits))    )
        plt.grid()
    
        x_plot = time[t_w_ind:]
        y_plot = qAvg
        y_plot_err = qStd/np.sqrt(NumCells-1)
        plt.errorbar(x_plot, y_plot, yerr=y_plot_err, fmt="o", markersize =0.8, elinewidth=0.2, label=r'$\langle AR \rangle$', zorder=1)
        plt.plot(x_plot, qAvgCumAvg, label='cum avg', zorder=2)
        plt.axvline(x=WaitingMCSteps, linewidth=1, linestyle='dashed', c ='black')
        # plt.scatter(x_plot, y_plot, label='MSD', c='black', s=0.4)
        plt.xscale("log")
        # plt.yscale("log")
        plt.legend()
        
        
        if saveSwitch:
            plt.savefig(ppPlotFolderName+'/'+'Hex_eq.png',dpi=dpiRes)
        
    return 

def g_of_r(saveSwitch):
    
    nBins = 200
    rCutOff = 6*np.sqrt(AvgCellArea)
    dr = rCutOff/nBins
    
    N_eqSamples = len(eqSamplingTimes.snapshotCList)
    
    g_of_r_all = np.zeros([N_eqSamples, nBins])
    
    for eqSampleC in range(N_eqSamples):
        
        snapshotC = eqSamplingTimes.snapshotCList[eqSampleC]
        time = eqSamplingTimes.timesList[eqSampleC]
        sampleC = eqSamplingTimes.sampleCList[eqSampleC]
        bunchC = eqSamplingTimes.bunchCList[eqSampleC]
        
    
        x_file_name = ppDataFolderName+"/xComBunch_"+str(bunchC)+".csv"
        y_file_name = ppDataFolderName+"/yComBunch_"+str(bunchC)+".csv"
        
        x_set = (np.loadtxt(x_file_name, delimiter=',' , dtype=float)[1:, sampleC]) % Lx
        y_set = (np.loadtxt(y_file_name, delimiter=',' , dtype=float)[1:, sampleC]) % Ly
        
        # g_of_r_sample = g_of_r_calc(x_set, y_set, Lx, Ly, nBins, rCutOff, "g_of_r_"+str(eqSampleC)+".PNG")
        g_of_r_sample = g_of_r_calc(x_set, y_set, Lx, Ly, nBins, rCutOff, None)
        
        print(time)
        g_of_r_all[eqSampleC, :] = g_of_r_sample
        
        # plt.figure()
        # x_plot = np.arange(0+0.5*dr , rCutOff+0.0001*dr , dr)/ np.sqrt(AvgCellArea)
        # y_plot = g_of_r_sample
        # y_plot_err = g_of_r_sample*0
        # plt.errorbar(x_plot, y_plot, yerr=y_plot_err, fmt="o", markersize =0.8, elinewidth=0.2, label='g(r)', zorder=1)
        # plt.xlabel(r"$r/\sqrt{A_0}$")
        # plt.ylabel(r"$g(r)$")
        # plt.axhline(y=1, linewidth=1, linestyle='dashed', c ='black')
        # plt.legend()
        # plt.grid()
        # plt.savefig("g_of_r_"+str(eqSampleC)+".PNG", dpi=500)
        # plt.show()
        
    if saveSwitch:
        bins_vec = np.arange(0+0.5*dr , rCutOff+0.0001*dr , dr)/ np.sqrt(AvgCellArea)
        g_of_r_all_with_bins = np.zeros([N_eqSamples+1, len(bins_vec)])
        g_of_r_all_with_bins[0,:] = bins_vec
        g_of_r_all_with_bins[1:,:] = g_of_r_all
        np.savetxt(ppDataFolderName + "/g_of_r_eq.csv", g_of_r_all_with_bins, delimiter=',', fmt='%1.5f')
    
    if plot_switches_dict['g_of_r']:
        plt.figure()
        x_plot = bins_vec
        y_plot = np.mean(g_of_r_all, axis=0)
        y_plot_err = np.std(g_of_r_all, axis=0)/np.sqrt(N_eqSamples)
        # y_plot_err = np.std(g_of_r_all, axis=0)
        plt.errorbar(x_plot, y_plot, yerr=y_plot_err, fmt="o", markersize =0.8, elinewidth=0.2, label='g(r)', zorder=1)
        plt.xlabel(r"$r/\sqrt{A_0}$")
        plt.ylabel(r"$g(r)$")
        plt.axhline(y=1, linewidth=1, linestyle='dashed', c ='black')
        plt.legend()
        plt.grid()
        plt.savefig(ppPlotFolderName+"/g_of_r.PNG", dpi=500)

    return

# cellsPlotterCSV(initFolderName+'/sigma_init.csv',ppPlotFolderName+'/initial_config.png' , 1)
# cellsPlotterCSV(mainResumeFolderName+'/sigmaMatLS.csv',ppPlotFolderName+'/final_config.png' , 1)

plot_switches_dict = dict()

with open('plot_switches.txt') as plot_switches:
        line = True
        while line:
            line = plot_switches.readline()
            line_list = line.split(": ")
            try:
                plot_switches_dict[line_list[0]]=int(line_list[1][0])
            except IndexError:
                break
        


print('xyComMaker: running...')
xCom, yCom = xyComMaker()
print('xyComMaker: done')

# temporalSnapshots = snapshotsMaker()    
cmap = plt.get_cmap('rainbow')

if plot_switches_dict['init']:
    print('cellsPlotter (init) : running...')
    initial_config  = np.loadtxt(initFolderName+'/sigma_init.csv',  delimiter=',', dtype=int)
    initial_compart = np.loadtxt(initFolderName+'/compart_init.csv',  delimiter=',', dtype=int)
    if lattice_switch == 'irr':
        cellsPlotterIrr(initial_config, initial_compart, ppPlotFolderName+'/initial_config.png' , 1, cmap, 1, xCom[:,0], yCom[:,0])
    elif lattice_switch == 'sq':
        cellsPlotterCSV(initial_config, ppPlotFolderName+'/initial_config.png' , 1)
    print('cellsPlotter (init) : done.')

if plot_switches_dict['final']:
    print('cellsPlotter (final) : running...')
    final_config  = np.loadtxt(mainResumeFolderName+'/sigmaMatLS.csv',  delimiter=',', dtype=int)
    final_compart = np.loadtxt(mainResumeFolderName+'/compartMatLS.csv',  delimiter=',', dtype=int)
    if lattice_switch == 'irr':
        cellsPlotterIrr(final_config, final_compart, ppPlotFolderName+'/final_config.png' , 1, cmap, 1, xCom[:,-1], yCom[:,-1])
    elif lattice_switch == 'sq':
        cellsPlotterCSV(final_config, ppPlotFolderName+'/final_config.png' , 1)
    print('cellsPlotter (final) : done.')


print('temporalPlotsV2 : running...')
# temporalPlots(1)
temporalPlotsV2(1)
print('temporalPlotsV2 : done')
# trajectoriesPlot(1)

if plot_switches_dict['traj']:
    print('trajectoriesV2 : running...')
    trajectoriesV2(1)
    print('trajectoriesV2 : done')

print('qFunc : running...')
qFunc(1)
print('qFunc : done')


print('ARFunc : running...')
ARFunc(1)
print('ARFunc : done')

print('circFunc : running...')
circFunc(1)
print('circFunc : done')


print('orderParametersFunc: running...')
if lattice_switch == 'irr':
    orderParametersFunc_Irr(1)
elif lattice_switch == 'sq':
    orderParametersFunc(1)
print('orderParametersFunc: done')


print('temporalHex : running...')
# temporalPlots(1)
temporalHex(1)
print('temporalHex : done')
# trajectoriesPlot(1)

#
#print('g(r) : running...')
## temporalPlots(1)
#g_of_r(1)
#print('g(r) : done')
#
