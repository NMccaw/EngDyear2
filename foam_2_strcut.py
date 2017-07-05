#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul  5 13:50:24 2017

@author: nicholas
"""

import sys
import csv
import numpy as np
from timo2 import BEMT_FEA_mesh_interp
from timo2 import finite_element_1d

N_elements_FEA = 1000
L = 1

simulation_directory = 'NACA0020_00deg/' 
simulation_file_path = '/home/nicholas/Documents/Summer_Proj/Simulations/Examples/'
path_to_pressure_lower = 'postProcessing/sampleDict/2924/surfaceProbesLower_p_k_omega.xy'
path_to_pressure_upper = 'postProcessing/sampleDict/2924/surfaceProbesUpper_p_k_omega.xy'

path = simulation_file_path + simulation_directory
sys.path.insert(0,path)

pressure_upper_list = []
pressure_lower_list = []

with open(path+path_to_pressure_upper) as f:
    reader = csv.reader(f,delimiter='\t')        
    for row in reader:
        pressure_upper_list.append(float(row[3]))
        
with open(path+path_to_pressure_lower) as f:
    reader = csv.reader(f,delimiter='\t')        
    for row in reader:
        pressure_lower_list.append(float(row[3]))
        
pressure_upper = np.zeros(len(pressure_upper_list))
pressure_lower = np.zeros(len(pressure_lower_list))

N_elements_BEMT = len(pressure_upper)

for i in range(len(pressure_upper)):
    pressure_upper[i] = pressure_upper_list[i]
    pressure_lower[i] = pressure_lower_list[i]
    
q_upper = BEMT_FEA_mesh_interp(N_elements_FEA,N_elements_BEMT,pressure_upper,L)
q_lower = BEMT_FEA_mesh_interp(N_elements_FEA,N_elements_BEMT,pressure_lower,L)
q = q_upper+q_lower

E = 5e6
nu = 0.3
kappa = 5/6
h = 0.2
b = 1

#U,free_dof = finite_element_1d(N_elements_FEA,q,L,E,nu,kappa,h,b)
up_points = []
low_points = []
keep_flag_up = 0
keep_flag_low = 0
with open(path + 'system/sampleDict') as f:
    for i,line in enumerate(f):
        if keep_flag_up == 1 and '(' and ')' in line:
            up_points.append(line.strip())
        if keep_flag_low == 1 and '(' and ')' in line:
            low_points.append(line.strip())
        if 'surfaceProbesUpper' in line:
            keep_flag_up = 1
        if 'surfaceProbesLower' in line:
            keep_flag_low = 1
            keep_flag_up = 0

for line in up_points:
    if line == ');':
        up_points.remove(line)

for line in low_points:
    if line == ');':
        low_points.remove(line)
        




