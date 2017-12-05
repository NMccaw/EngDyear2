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
from timo_explicit_3D import timo_3D_explicit
# define paths to CFD data
path2points = '/home/nicholas/Documents/Summer_Proj/Simulations/Examples/3D/initial/system/'
working_dir = '/home/nicholas/Documents/Summer_Proj/Simulations/Examples/3D/initial/postProcessing/sampleDict/5983/'
sys.path.insert(0,working_dir)

rho = 1
thickness = 0.006
N_probes = 300
N_elements = 1000
Uref = 10
AoA = 9
E_axis = 0.5
# structural data
kappa = 5/6
nu = 0.3 
E = 207e9
# wing shape according to processResults.py file
chord = 0.667
wing_span = 1
# define the span points where the pressure probes are
pressure_span_probes = [7,23,40,53,70]
for span in  np.linspace(80,100,21):
    pressure_span_probes.append(int(span))

pressure_upper = np.zeros(((len(pressure_span_probes)),N_probes))

pressure_lower = np.zeros_like(pressure_upper)

force_upper = np.zeros(len(pressure_span_probes))

force_lower = np.zeros_like(force_upper)
points_x_up = np.zeros(N_probes)
points_y_up = np.zeros_like(points_x_up)

points_x_low = np.zeros(N_probes)
points_y_low = np.zeros_like(points_x_up)

h_upper = np.zeros(N_probes - 1)
h_lower = np.zeros(N_probes - 1)
integrand_upper = np.zeros((len(pressure_span_probes),len(h_upper)))
integrand_lower = np.zeros_like(integrand_upper)
integrand_upper_cop = np.zeros_like(integrand_upper)
integrand_lower_cop = np.zeros_like(integrand_upper)


# read in pressure data from CFD results
for i,span in enumerate(pressure_span_probes):
    
    with open(working_dir+'surfaceProbesUpper'+str(span)+'_p_k_omega.xy') as f:
        reader = csv.reader(f,delimiter='\t')        
        for j,row in enumerate(reader):
            pressure_upper[i,j] = float(row[3])*rho
            
    with open(working_dir+'surfaceProbesLower'+str(span)+'_p_k_omega.xy') as f:
        reader = csv.reader(f,delimiter='\t')        
        for j,row in enumerate(reader):
            pressure_lower[i,j] = float(row[3])*rho
            
            

pressure_diff = np.abs(pressure_upper - pressure_lower)
up_points = []
low_points = []
keep_flag_up = 0
keep_flag_low = 0
# extract the data points which defines the upper and lower pressure probes
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
# extract coordinates for upper probes
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
# extract coordinates for lower points
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

# compute the distance between each probe point. For integration             
for i in range(len(points_x_up)-1):
    
    h_upper[i] = np.sqrt((points_x_up[i] - points_x_up[i+1])**2 + (points_y_up[i] - points_y_up[i+1])**2)
    h_lower[i] = np.sqrt((points_x_low[i] - points_x_low[i+1])**2 + (points_y_low[i] - points_y_low[i+1])**2)

    
# integrate pressue over line
integral_x_p_up = np.zeros_like(force_upper)
integral_x_p_low = np.zeros_like(force_lower)
# set points along x_axis to find c.o.p
cop_points_finder = np.linspace(0,1,len(points_x_up))
for j in range(len(pressure_span_probes)):
    for i in range(len(h_upper)):
        integrand_upper[j,i] = h_upper[i]/2 * (pressure_upper[j,i] + pressure_upper[j,i+1])
        integrand_lower[j,i] = h_lower[i]/2 * (pressure_lower[j,i] + pressure_lower[j,i+1])
        
        integrand_upper_cop[j,i] = h_upper[i]/2 * (pressure_upper[j,i]*cop_points_finder[i] + pressure_upper[j,i+1]*cop_points_finder[i+1])
        integrand_lower_cop[j,i] = h_lower[i]/2 * (pressure_lower[j,i]*cop_points_finder[i]  + pressure_lower[j,i+1]*cop_points_finder[i+1] )
        
    force_upper[j] = np.sum(integrand_upper[j,:])
    force_lower[j] = np.sum(integrand_lower[j,:])
    
    integral_x_p_up[j] = np.sum(integrand_upper_cop[j,:])
    integral_x_p_low[j] = np.sum(integrand_lower_cop[j,:])
# calculate centre of pressure    
cop = (integral_x_p_up+integral_x_p_low)/(force_upper+force_lower)

x_bar = (E_axis - cop) *chord # distance from elastic axis to c.o.p

h_span = np.zeros(len(pressure_span_probes))
for i in range(len(pressure_span_probes)-1):
    h_span[i+1] = np.abs(pressure_span_probes[i] - pressure_span_probes[i+1]) + h_span[i]


lines = ['k--','r--','k-','r-','b--','b-','g-']

plt.figure(2)
#plt.plot(pressure_upper[0,:])
for i in range(pressure_upper.shape[0]//4):
    plt.plot(np.linspace(0,1,pressure_upper.shape[1]),-pressure_upper[i,:]/(0.5*rho*Uref**2), lines[i],label = 'span = '+str(pressure_span_probes[i])+'%')
    plt.plot(np.linspace(0,1,pressure_upper.shape[1]),-pressure_lower[i,:]/(0.5*rho*Uref**2), lines[i])

plt.ylabel('-Cp')
plt.xlabel('x/c')
plt.legend()
plt.savefig('pressure_dist1.pdf')




plt.figure(20)
#plt.plot(pressure_upper[0,:])
for i in range(pressure_upper.shape[0]//4,2*pressure_upper.shape[0]//4):
    plt.plot(np.linspace(0,1,pressure_upper.shape[1]),-pressure_upper[i,:]/(0.5*rho*Uref**2), lines[i-pressure_upper.shape[0]//4],label = 'span = '+str(pressure_span_probes[i])+'%')
    plt.plot(np.linspace(0,1,pressure_upper.shape[1]),-pressure_lower[i,:]/(0.5*rho*Uref**2), lines[i-pressure_upper.shape[0]//4])

plt.ylabel('-Cp')
plt.xlabel('x/c')
plt.legend()
plt.savefig('pressure_dist2.pdf')

plt.figure(200)
#plt.plot(pressure_upper[0,:])
for i in range(2*pressure_upper.shape[0]//4,3*pressure_upper.shape[0]//4):
    plt.plot(np.linspace(0,1,pressure_upper.shape[1]),-pressure_upper[i,:]/(0.5*rho*Uref**2), lines[i-2*pressure_upper.shape[0]//4],label = 'span = '+str(pressure_span_probes[i])+'%')
    plt.plot(np.linspace(0,1,pressure_upper.shape[1]),-pressure_lower[i,:]/(0.5*rho*Uref**2), lines[i-2*pressure_upper.shape[0]//4])

plt.ylabel('-Cp')
plt.xlabel('x/c')
plt.legend()
plt.savefig('pressure_dist3.pdf')

plt.figure(2000)
#plt.plot(pressure_upper[0,:])
for i in range(3*pressure_upper.shape[0]//4,4*pressure_upper.shape[0]//4-1):
    plt.plot(np.linspace(0,1,pressure_upper.shape[1]),-pressure_upper[i,:]/(0.5*rho*Uref**2), lines[i-3*pressure_upper.shape[0]//4],label = 'span = '+str(pressure_span_probes[i])+'%')
    plt.plot(np.linspace(0,1,pressure_upper.shape[1]),-pressure_lower[i,:]/(0.5*rho*Uref**2), lines[i-3*pressure_upper.shape[0]//4])

plt.ylabel('-Cp')
plt.xlabel('x/c')
plt.legend()

plt.savefig('pressure_dist4.pdf')



plt.figure(10)
#plt.plot(force_upper)
#plt.plot(force_lower)
force = force_upper - force_lower
plt.plot(pressure_span_probes,force*np.cos(AoA*np.pi/180)/(0.5*rho*Uref**2))
plt.xlabel('% span')
plt.ylabel('Normal Coefficient Cn')
plt.ylim([0,0.5])
plt.savefig('normal_dist.pdf')





def FEA_mesh_interp(N_elements_FEA,N_elements_BEMT,Cl,L,h_span):
    '''
    Interpolates between the FEA mesh and BEMT mesh. FEA mesh will be much finer
    than BEMT mesh so the force vectors will differ between the two. This function
    maps the force vectors of the BEMT to the FEA elements.
    
    Parameters
    ----------
    N_elements_FEA: int
        number of FEA elements
    N_elements_BEMT: int
        number of BEMT elements
    Cl: array of floats
        lift coefficient or force vector at each blade element
    L: float
        length of beam/blade.
        
    Returns
    --------
    q: array of floats
        force on each FEA element
    '''
    
    LM_BEMT = np.zeros((2,N_elements_BEMT))
    LM_FEA = np.zeros((2,N_elements_FEA))
    for e,h in enumerate(h_span):
        
        LM_BEMT[0,e] = e*L*h/100
        LM_BEMT[1,e] = (e+1)*L*h/100
        if e > 0:
            LM_BEMT[0,e] = L*h_span[e-1]/100
            LM_BEMT[1,e] = L*h_span[e]/100
            
    for e in range(N_elements_FEA):
        LM_FEA[0,e] = e*L/N_elements_FEA
        LM_FEA[1,e] = (e+1)*L/N_elements_FEA
    
    q = np.zeros(N_elements_FEA)
    residual = 0
    for j in range(N_elements_FEA):   
        k = j 
        for i in range(N_elements_BEMT):
            if LM_FEA[0,k] >= LM_BEMT[0,i] and LM_FEA[1,k] <= LM_BEMT[1,i]:
                q[k] = Cl[i]
                k+=1
            elif LM_FEA[0,k] < LM_BEMT[1,i] and LM_FEA[1,k] > LM_BEMT[1,i]:
                residual = abs(LM_FEA[1,j] - LM_BEMT[1,i])
                q[k] = residual*Cl[i] + (1-residual)*Cl[i+1]
                k+=1
    return q

moment = x_bar*force
Cm = moment/(0.5*rho*Uref**2)
Cn = force*np.cos(AoA*np.pi/180)/(0.5*rho*Uref**2)
print('Cn = {}'.format(Cn[3]))
print('Cmx = {}'.format(Cm[3]))
print('C.O.P = {}'.format(cop[3]))
moment_fea = FEA_mesh_interp(N_elements,pressure_upper.shape[0],moment,wing_span,pressure_span_probes)
force_fea = FEA_mesh_interp(N_elements,pressure_upper.shape[0],force,wing_span,pressure_span_probes)

b = np.ones_like(force_fea)
q = np.zeros((N_elements,6))
t = np.ones_like(b) * thickness
q[:,1] = -force_fea
q[:,3] = -moment_fea
#U,free_dof = finite_element_1d(N_elements,force_fea,wing_span,E,nu,kappa,t,b)
U = timo_3D_explicit(N_elements,q,E,nu,kappa,t,b,wing_span)

plt.figure(3)
plt.plot(np.linspace(0,1,len(U[1::6])),U[1::6])
plt.ylabel('Deflection (m) ')
plt.xlabel('y/span')
plt.tight_layout()
plt.savefig('deflection.pdf')
plt.figure(4)
plt.plot(np.linspace(0,1,len(U[3::6])),U[3::6])
plt.ylabel(r'Change in Angle of attack $\alpha$ (rad)')
plt.xlabel('y/span')
plt.tight_layout()
plt.savefig('AoA.pdf')


#plt.figure(4)
#plt.plot(np.linspace(0,1,pressure_upper.shape[0]),-force)
#plt.xlabel('y/span')
#plt.ylabel('Lift Force (N)')
#plt.savefig('lift_dist.pdf')