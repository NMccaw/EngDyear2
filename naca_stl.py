#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug  2 16:11:19 2017

@author: nicholas
"""

from stl import mesh
from mpl_toolkits import mplot3d
from matplotlib import pyplot as plt
import numpy as np 
import stl


N_elements = 10

path2points = '/home/nicholas/Documents/Summer_Proj/Simulations/Examples/3D/initial/constant/triSurface/'
file = 'naca_ascii.stl'

naca = mesh.Mesh.from_file(path2points+file)
foil_section_fraction = 1/N_elements
 
naca.vectors[:,:,2] *= foil_section_fraction


indice = naca.vectors.shape[0]

to_add = np.zeros((N_elements,naca.vectors.shape[0],3))

added = naca.vectors


angles_all = U[5::6]
angles = angles_all[0::len(U[3::6])//N_elements]

displacements_all = U[1::6]

displacements = displacements_all[0::len(U[3::6])//N_elements] * 10000

aoa_all = U[3::6]
AoA = aoa_all[0::len(U[3::6])//N_elements] *10000

print(displacements)
c = 0.667

  
vertices = np.zeros(((N_elements)*naca.vectors.shape[0],3,3))
for i in range(N_elements):  
    to_add[i,:,:] = naca.vectors[:,:,2] + (foil_section_fraction*i)
    #naca.rotate([0,0,1],np.deg2rad(00))
    naca = mesh.Mesh.from_file(path2points+file)
    naca.translate([-0.2*c,0,0])
    naca.rotate([0,0,1],-AoA[i])
    naca.translate([0.2*c,0,0])
    naca.rotate([1,0,0],angles[i])
    naca.translate([0,displacements[i],0])
    naca.vectors[:,:,2] = ((foil_section_fraction)*naca.vectors[:,:,2]) + (foil_section_fraction*i)
    vertices[i*indice:(i+1)*indice,:,:] = naca.vectors
    
init_AoA = np.deg2rad(9.6)

figure = plt.figure()
axes = mplot3d.Axes3D(figure)

axes.add_collection3d(mplot3d.art3d.Poly3DCollection(vertices))
axes.scatter(-0.3288,displacements[0]+(0.7*0.667)*np.sin(AoA[0])-(0.667*0.5)*np.sin(init_AoA),0,c='r',marker='o')
axes.scatter(np.cos(init_AoA)*(-0.5*c) +np.sin(AoA[9])*(c*0.7)*np.sin(AoA[9]),displacements[9]-(0.667*0.5)*np.sin(init_AoA)-(0.7*c)*np.sin(AoA[9]),0.9,c='r',marker='x')
#scale = naca.points.flatten(-1)
#axes.auto_scale_xyz(scale, scale, scale)
axes.view_init(90,0)
plt.savefig('stl_AoA.pdf')

figure = plt.figure(2)
axes = mplot3d.Axes3D(figure)

axes.add_collection3d(mplot3d.art3d.Poly3DCollection(vertices))

#scale = naca.points.flatten(-1)
#axes.auto_scale_xyz(scale, scale, scale)
axes.scatter(-0.3288,displacements[0]+(0.7*0.667)*np.sin(AoA[0])+(0.5*0.667)*np.sin(init_AoA),0,c='r',marker='o')
axes.scatter(-0.3288,displacements[9],0.9,c='r',marker='o')
#axes.view_init(0,0)
plt.savefig('stl_deflection.pdf')


naca_deformed = mesh.Mesh(np.zeros(vertices.shape[0], dtype = mesh.Mesh.dtype))

naca_deformed.vectors = vertices
    
    
    
naca_deformed.save('test.stl',mode=stl.Mode.ASCII)

'''

naca_new = mesh.Mesh.from_file('test.stl')
figure = plt.figure(2)
axes = mplot3d.Axes3D(figure)

axes.add_collection3d(mplot3d.art3d.Poly3DCollection(naca_new.vectors))

#scale = naca.points.flatten(-1)
#axes.auto_scale_xyz(scale, scale, scale)
axes.view_init(30,30)

'''