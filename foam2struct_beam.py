#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 11 10:23:29 2017

@author: nicholas
"""

import sys
import csv
import numpy as np
from timo2 import BEMT_FEA_mesh_interp
from timo2 import finite_element_1d
from matplotlib import pyplot as plt

N_elements_FEA = 1000
N_elements_CFD = 10
L = 1
rho = 1000
Vref = 1 
b = 0.2 * np.ones(N_elements_CFD)
h = 0.006 * np.ones(N_elements_CFD)

E = 207e9 #67000e6
nu = 0.3
kappa = 5/6



time = '3178' # 3178 ,2924
simulation_directory = 'NACA0020_01deg/' 
simulation_file_path = '/home/nicholas/Documents/Summer_Proj/Simulations/Examples/'
forcesCoeff_dat = 'postProcessing/forcesCoeffs/0/forceCoeffs.dat'

path = simulation_file_path + simulation_directory
sys.path.insert(0,path)

with open(path + forcesCoeff_dat) as f:
    reader = csv.reader(f,delimiter='\t')
    for row in reader:       
        if 'Cl           ' in row:
            for i,value in enumerate(row):
                if value.strip() == 'Cl':
                    j = i
                    
cl = float(row[j])
 
lift_force = (0.5 * rho * Vref**2 * cl * b * L) / (L/N_elements_CFD)

lift = BEMT_FEA_mesh_interp(N_elements_FEA,N_elements_CFD,lift_force,L)

h = BEMT_FEA_mesh_interp(N_elements_FEA,N_elements_CFD,h,L)
b = BEMT_FEA_mesh_interp(N_elements_FEA,N_elements_CFD,b,L)


U,free_dof = finite_element_1d(N_elements_FEA,lift,L,E,nu,kappa,h,b)


plt.plot(free_dof[:N_elements_FEA],U[:N_elements_FEA]/L)
plt.title('Deflection of foil')
plt.ylabel('Deflection (m)')

plt.figure(2)
plt.plot(free_dof[:N_elements_FEA],U[N_elements_FEA:]/L)