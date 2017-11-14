#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 20 14:47:55 2017

Creates structural mesh for NACA 00XX foil in OpenFOAM.  

@author: nicholas
"""
import numpy as np
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

c = 1 # chord 
t = 0.2 # max thickness
wing_span = 0.67


x = np.linspace(0,1,100)
z_points = np.linspace(0,wing_span,100)


y = 5*t*(0.2969*np.sqrt(x/c) - 0.1260*(x/c) -0.3516*(x/c)**2 + 0.2843*(x/c)**3 - 0.1015*(x/c)**4)

y_upper = y 
y_lower = -y 

points_mat = np.zeros((3,100))
points_mat[0,:] = x
points_mat[1,:] = y
points_mat[2,:] = z_points
def array_checker(dictionary):
    new_value = '('
    for i in dictionary:
        if type(dictionary[i]) == np.ndarray:
            for j in range(len(dictionary[i])):
                if j == len(dictionary[i])-1:
                    new_value+= str(dictionary[i][j])
                else:
                    new_value+= str(dictionary[i][j]) + ' '
            dictionary[i] = new_value
            dictionary[i] += ')'
class make_file:
    
    def __init__(self, points_mat,file,value_dict,block_mesh_dict,snappy_dict):
        
        self.x = points_mat[0,:]
        self.y = points_mat[1,:]
        self.z = points_mat[2,:]
        self.f = open(file, "w")

        
        self.rho = rho
        self.visc = visc
        self.flowVelocity = flowVelocity
        self.magSpeed = magSpeed
        self.pressure = pressure
        self.turbulentKE = turbulentKE
        self.turbulentOmega = turbulentOmega
        self.turbulentOmegaWall = turbulentOmegaWall
        self.refArea = refArea
        self.refL = refL
    def new_line(self):
        self.f.write('\n')
        
    def tab(self):
        self.f.write('\t')
    def end(self):
        self.new_line()
        
        self.f.write('// ************************************************************************* //')
        self.f.close()
    def make_header(self):
        self.f.write('/*--------------------------------*- C++ -*----------------------------------*\ ' + '\n')
        self.f.write('| =========                 |                                                 |' + '\n')
        self.f.write('| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |' + '\n')
        self.f.write('|  \\    /   O peration     | Version:  4.1                                   |' + '\n')
        self.f.write('|   \\  /    A nd           | Web:      www.OpenFOAM.org                      |' + '\n')
        self.f.write('|    \\/     M anipulation  |                                                 |' + '\n')
        self.f.write('\*---------------------------------------------------------------------------*/' + '\n')
        
    def make_bloackMeshDict(self):
        
        self.f.write('FoamFile' + '\n')
        self.f.write('{' + '\n')
        self.f.write('version' + '\t' + '2.0;'+' \n')
        self.f.write('format' + '\t' + 'ascii;' +'\n')
        self.f.write('class' + '\t' + 'dictionary;' +'\n')
        self.f.write('object' + '\t' + 'blockMeshDict;' + '\n')
        self.f.write('convertToMeters 1;'+'\n')
        self.f.write('}'+'\n')
        self.f.write('vertices'+'\n')
        self.f.write('(' +'\n')
        for k in range(points_mat.shape[1]):
                for i in range(points_mat.shape[1]):
                    self.f.write('\t'+'('+str(self.x[i]) + ' ' + str(self.y[i]) + ' ' + str(self.z[k])+')' + '\n')
                    self.f.write('\t'+'('+str(self.x[i]) + ' ' + str(-self.y[i]) + ' ' + str(self.z[k])+')' + '\n')
        self.f.write(');' + '\n')
        self.f.write('blocks' + '\n' + '(' + '\n' )
        for j in range(0,points_mat.shape[1]**2,2):
            self.f.write('\t' +'hex '+'(' + str(j) +' '+ str(j+1) +' '+ str(j+3) +' '+ str(j+2)+' ' +str(j+200)+' '+str(j+201)+' ' +str(j+203)+' ' +str(j+202)+') ' + '(1 1 1)' +' simpleGrading (1 1 1)')
        
        self.f.write( ');' + '\n')
        self.f.write('boundary' + '\n' + '(' + '\n')
        self.f.write('root' +'\n' + '{' +'\n' +'type patch;' + '\n' + 'faces'+'\n' +'(' +'\n')
        for i in range(0,2*points_mat.shape[1],2):
            self.f.write('('+str(i) +' '+str(i+1) + ' ' +str(i+3)+ ' ' + str(i+2) +')' + '\n')
            
        self.f.write(');'+'\n' + '}'+'\n')
            
        self.f.write('upper_surface' +'\n' + '{' +'\n' +'type patch;' + '\n' + 'faces'+'\n' +'(' +'\n')
        for i in range(0,points_mat.shape[1]**2,2):
            self.f.write('('+str(i) +' '+str(i+2) + ' ' +str(i+202)+ ' ' + str(i+200) +')' + '\n')
        self.f.write(');'+'\n' + '}'+'\n')
            
        self.f.write(');' +'\n')
            
        self.f.write('mergePatchPairs'+'\n' + '(' +'\n' +');' +'\n')
        self.f.write('\n' + '//' + '\n' + '*************************************************************************' + ' // ')
        self.f.close()
        
    def make_run_conditions(self):
        
        self.f.write('// FLOW PROPERTIES')
        for i in value_dict:
            self.new_line()
            self.f.write(i)
            self.tab()
            self.f.write(str(value_dict[i])+';')

        self.f.write('\n')
        self.f.write('// BlOCK MESH SETTINGS')
        for i in block_mesh_dict:
            self.new_line()
            self.f.write(i)
            self.tab()
            self.f.write(str(block_mesh_dict[i])+';')
        
        
        self.f.write('\n')
        self.f.write('// SNAPPY SETTINGS')
        for i in snappy_dict:
            if 'Title' in i:
                self.new_line()
                self.f.write(snappy_dict[i])
            else:
                self.new_line()
                self.f.write(i)
                self.tab()
                self.f.write(str(snappy_dict[i])+';')
        
        
        
if __name__ == '__main__':
    rho = 1
    visc = 1.667e-5
    flowVelocity = np.array([-10,0,0])
    magSpeed = 10
    pressure = 0 
    turbulentKE = 0.667
    turbulentOmega = 0.2236
    turbulentOmegaWall = 3.0e-014
    refArea = 0.667
    refL = 0.667
    
    # black mesh settings
    nx = 50
    ny = 12
    nz = 5
    xu = 4
    xl = -8
    yu = 1.75
    yl = -1.75
    zl = 0
    zu = 2.5
    
    
    
    # SNAPPY SETTINGS
    snappy_dict = {}
    #  general
    featAngleRef = 65
    nCellBetLvl =  3
    locInMesh = np.array([1.1, 0.011, 0.011]) 
    nSnapIter = 1#//5;
    yTE = -0.05554 #negative y for positive angle of yaw
    yLE = 0.05554
    
    snappy_dict['featAngleRef'] = featAngleRef
    snappy_dict['nCellBetLvl'] = nCellBetLvl
    snappy_dict['locInMesh'] = locInMesh
    snappy_dict['nSnapIter'] = nSnapIter
    snappy_dict['yTE'] = yTE
    snappy_dict['yLE'] = yLE
    
    # layers
    expRat =  1.1
    finalThick = 0.2
    snappy_dict['Title layers'] = '// layers'

    snappy_dict['expRat'] = expRat
    snappy_dict['finalThick'] = finalThick

    # refinement cylinder - Tip refinement
    xuCyl = 0.35
    xlCyl = -0.35
    zCyl = 1
    rCyl = 0.1
    snappy_dict['Title refinement cylinder tip'] = '// refinement cylinder - Tip refinement'
    snappy_dict['xuCyl'] = xuCyl
    snappy_dict['xlCyl'] = xlCyl
    snappy_dict['zCyl'] = zCyl
    snappy_dict['rCyl'] = rCyl
    # refinement cylinder 2 - Vortex refinement
    xuCyl2 = -0.3
    xlCyl2 = -3.0
    yuCyl2 = -0.0181
    ylCyl2 = -0.0153
    zuCyl2 = 0.9671
    zlCyl2 = 0.9268
    rCyl2 = 0.1
    snappy_dict['Title refinement cylinder 2'] = '// refinement cylinder 2 - Vortex refinement'
    snappy_dict['xuCyl2'] = xuCyl2
    snappy_dict['xlCyl2'] = xlCyl2
    snappy_dict['yuCyl2'] = yuCyl2
    snappy_dict['ylCyl2'] = ylCyl2
    snappy_dict['zuCyl2'] = zuCyl2
    snappy_dict['zlCyl2'] = zlCyl2
    snappy_dict['rCyl2'] = rCyl2
    # refinement cylinder 3 - TE refinement
    xCyl3 =  -0.3335
    zlCyl3 = 0
    zuCyl3 = 1.01
    rCyl3 = 0.1
    snappy_dict['Title refinement cylinder 3'] = '// refinement cylinder 3 - TE refinement'
    snappy_dict['xCyl3'] = xCyl3
    snappy_dict['zlCyl3'] = zlCyl3
    snappy_dict['zuCyl3'] = zuCyl3
    snappy_dict['rCyl3'] = rCyl3
    # refinement box LE
    xuBoxLE = 0.5
    xlBoxLE = 0
    yuBoxLE = 0.1
    ylBoxLE = -0.07
    zuBoxLE = 1.05
    zlBoxLE = 0
    snappy_dict['Title refinement box LE'] = '// refinement box LE'
    snappy_dict['xuBoxLE'] = xuBoxLE
    snappy_dict['xlBoxLE'] = xlBoxLE
    snappy_dict['yuBoxLE'] = yuBoxLE
    snappy_dict['ylBoxLE'] = ylBoxLE
    snappy_dict['zuBoxLE'] = zuBoxLE
    snappy_dict['zlBoxLE'] = zlBoxLE
    
    # refinement box TE
    xuBoxTE = 0
    xlBoxTE = -1
    yuBoxTE = 0.1
    ylBoxTE = -0.1
    zuBoxTE =  1.0
    zlBoxTE =  0
    snappy_dict['Title refinement box TE'] = '// refinement box TE'
    snappy_dict['xuBoxTE'] = xuBoxTE
    snappy_dict['xlBoxTE'] = xlBoxTE
    snappy_dict['yuBoxTE'] = yuBoxTE
    snappy_dict['ylBoxTE'] = ylBoxTE
    snappy_dict['zuBoxTE'] = zuBoxTE
    snappy_dict['zlBoxTE'] = zlBoxTE
    # refinement levels
    seisLvlMin =	1
    seisLvlMax =	2
    nSeisLayers = 1#//30;
    snappy_dict['Title refinement levels'] = '// refinement levels'
    snappy_dict['seisLvlMin'] = seisLvlMin
    snappy_dict['seisLvlMax'] = seisLvlMax
    snappy_dict['nSeisLayers'] = nSeisLayers
    
    seisMaxLvlRefDist = 0.01
    seisMinLvlRefDist = 0.02
    snappy_dict['seisMaxLvlRefDist'] = seisMaxLvlRefDist
    snappy_dict['seisMinLvlRefDist'] = seisMinLvlRefDist
    refCylLvl = 1 #6;
    refCylLvl2 = 1 #5;
    refCylLvl3 = 5 #8;
    refBoxLELvl = 1 #5;
    refBoxTELvl =  4
    snappy_dict['refCylLvl'] = refCylLvl
    snappy_dict['refCylLvl2'] = refCylLvl2
    snappy_dict['refCylLvl3'] = refCylLvl3
    snappy_dict['refBoxLELvl'] = refBoxLELvl
    snappy_dict['refBoxTELvl'] = refBoxTELvl
    
    
    
    
    value_dict = {}
    
    value_dict['rho'] = rho
    value_dict['visc'] = visc
    value_dict['flowVelocity'] = flowVelocity
    value_dict['magSpeed'] = magSpeed
    value_dict['pressure'] = pressure
    value_dict['turbulentKE'] = turbulentKE
    value_dict['turbulentOmega'] = turbulentOmega
    value_dict['turbulentOmegaWall'] = turbulentOmegaWall
    value_dict['refArea'] = refArea
    value_dict['refL'] = refL
    
    block_mesh_dict = {}
    
    block_mesh_dict['nx'] = nx
    block_mesh_dict['ny'] = ny
    block_mesh_dict['nz'] = nz
    block_mesh_dict['xu'] = xu
    block_mesh_dict['xl'] = xl
    block_mesh_dict['yu'] = yu
    block_mesh_dict['yl'] = yl
    block_mesh_dict['zl'] = zl
    block_mesh_dict['zu'] = zu
    array_checker(value_dict)
    array_checker(block_mesh_dict)
    array_checker(snappy_dict)
    file = make_file(points_mat,'runConditions_test',value_dict,block_mesh_dict,snappy_dict)
    file.make_header()
    file.make_run_conditions()
    file.end()