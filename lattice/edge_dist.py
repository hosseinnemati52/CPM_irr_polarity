#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 19 17:19:06 2024

@author: hossein
"""








import numpy as np
import math
import matplotlib.pyplot as plt

edges = []

with open('walls.dat', 'r') as edges_file:

    for line in edges_file:
        # do something with the line
        edges = edges + [float(number) for number in line.split('\t') if number.strip()][3:]

np.savetxt("edges_mean.csv", X= np.array([np.mean(edges)]), fmt = "%.15f")