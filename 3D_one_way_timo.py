#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 14 13:36:38 2017

@author: nicholas
"""

import sys
import csv
import numpy as np
from timo2 import BEMT_FEA_mesh_interp
from timo2 import finite_element_1d
from matplotlib import pyplot as plt
path2points = '/home/nicholas/Documents/Summer_Proj/Simulations/Examples/3D/initial/system/'
working_dir = '/home/nicholas/Documents/Summer_Proj/Simulations/Examples/3D/initial/postProcessing/sets/5983/'
sys.path.insert(0,working_dir)

rho = 1000

N_probes = 300

pressure_span_probes = [7,23,40,53,70,83,97]

pressure_upper = np.zeros(((len(pressure_span_probes)),N_probes))

pressure_lower = np.zeros_like(pressure_upper)

force = np.zeros(len(pressure_span_probes))

# wing shape according to processResults.py file
chord = 1
wing_span = 0.667

for i,span in enumerate(pressure_span_probes):
    
    with open(working_dir+'surfaceProbesLower'+str(span)+'_p_k_omega.xy') as f:
        reader = csv.reader(f,delimiter='\t')        
        for j,row in enumerate(reader):
            pressure_upper[i,j] = float(row[3])*rho
            
    with open(working_dir+'surfaceProbesLower'+str(span)+'_p_k_omega.xy') as f:
        reader = csv.reader(f,delimiter='\t')        
        for row in reader:
            pressure_lower[i,j] = float(row[3])*rho




pressure_diff = np.abs(pressure_upper - pressure_lower)
up_points = []
low_points = []
keep_flag_up = 0
keep_flag_low = 0
with open(path2points+'sampleDict') as f:
    for i,line in enumerate(f):
        if keep_flag_up == 1 and '(' and ')' in line:
            up_points.append(line.strip())
        if keep_flag_low == 1 and '(' and ')' in line:
            low_points.append(line.strip())
        if 'surfaceProbesUpper23' in line:
            keep_flag_up = 1
            keep_flag_low = 0
        if 'surfaceProbesLower23' in line:
            keep_flag_low = 1
            keep_flag_up = 0
        if 'surfaceProbesUpper' in line:
            keep_flag_low = 0

print(len(low_points))

points_x_up = np.zeros(N_probes)
points_y_up = np.zeros_like(points_x_up)

points_x_low = np.zeros(N_probes)
points_y_low = np.zeros_like(points_x_up)

h = np.zeros(N_probes - 1)

for i,line in enumerate(up_points):
    newline = ''
    char_keep = ''
    for letter in line:
        if letter == '(' or letter == ')' or letter == ';' or letter == ');':
            pass
        else:
            char_keep = letter
        newline += str(char_keep) 
    if newline.split() != []:
        i -= 1
        points_x_up[i] = float(newline.split()[0])
        points_y_up[i] = float(newline.split()[1])

for i,line in enumerate(low_points):
    newline = ''
    char_keep = ''
    for letter in line:
        if letter == '(' or letter == ')' or letter == ';' or letter == ');':
            pass
        else:
            char_keep = letter
        newline += str(char_keep) 
    if newline.split() != []:
        i -= 1
        points_x_low[i] = float(newline.split()[0])
        points_y_low[i] = float(newline.split()[1])




for i in range(len(points_x_up)-1):
    
    h[i] = np.sqrt((points_x_up[i] - points_x_up[i+1])**2 + (points_x_up[i] - points_x_up[i+1])**2)
integrand = np.zeros((len(pressure_span_probes),len(h)))
# integrate pressue over surface
for j in range(len(pressure_span_probes)):
    for i in range(len(h)):
        integrand[j,i] = h[i]/2 * (pressure_upper[j,i] + pressure_upper[j,i+1])
    force[j] = np.sum(integrand[j,:])

plt.plot(force)    
thickness = abs(points_y_up - points_y_low)

    
    
    
    