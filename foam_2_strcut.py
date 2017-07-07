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
from matplotlib import pyplot as plt

N_elements_FEA = 1000
L = 1
time = '3178' # 3178 ,2924
simulation_directory = 'NACA0020_01deg/' 
simulation_file_path = '/home/nicholas/Documents/Summer_Proj/Simulations/Examples/'
path_to_pressure_lower = 'postProcessing/sampleDict/'+ time + '/surfaceProbesLower_p_k_omega.xy'
path_to_pressure_upper = 'postProcessing/sampleDict/'+ time +'/surfaceProbesUpper_p_k_omega.xy'

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
q = q_upper - q_lower


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

points_x_up = np.zeros_like(pressure_upper)
points_y_up = np.zeros_like(points_x_up)


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
        
points_x_low = np.zeros_like(pressure_upper)
points_y_low = np.zeros_like(points_x_up)

        
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
        
    
    
h = np.abs(points_y_low - points_y_up) + 0.000001

h_fea = BEMT_FEA_mesh_interp(N_elements_FEA,N_elements_BEMT,h,L) 

E = 67000e6
nu = 0.3
kappa = 5/6

b = np.ones_like(h_fea) 
U,free_dof = finite_element_1d(N_elements_FEA,-q,L,E,nu,kappa,h_fea,b)
plt.plot(points_y_up)
plt.plot(points_y_low) 
plt.plot(np.linspace(0,300,num = 1000),(-q_upper/max(q_upper))+1 )   
plt.plot(np.linspace(0,300,num = 1000),(q_lower/max(q_lower))-1 )  
plt.title('Geometry and pressure fields')
plt.savefig('geometry.eps')
plt.figure(2)   
plt.plot(free_dof[:N_elements_FEA],U[:N_elements_FEA])
plt.title('Deflection of foil')
plt.ylabel('Deflection (m)')
plt.savefig('deflection.pdf')
plt.figure(3)
plt.plot(-q_upper,label = 'Upper surface')
plt.plot(-q_lower,label = 'Lower surface')
plt.legend()
plt.title('-Pressure along foil')
plt.ylabel('-Pressure (Pa)')
plt.savefig('pressure_up_down.pdf')
plt.figure(4)
plt.plot(q)
plt.title('combined pressure acting on foil')
plt.ylabel('Pressure (Pa)')
plt.savefig('cob_pressure.pdf')
plt.figure(5)
plt.plot(free_dof[N_elements_FEA:],U[N_elements_FEA:])
plt.title('Rotation of foil')
plt.ylabel('rotation (o)')
plt.savefig('rotation.pdf')
