#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug  1 18:04:24 2017

@author: nicholas
"""
import numpy as np
from matplotlib import pyplot as plt
#N_elements = 100
#q = np.zeros((N_elements,6))
#E = 2e5
#nu = 0.3
#L = 100
#kappa = 5/6
#h = 10 * np.ones(N_elements)
#b = h
#q[:,1] = 500
def timo_3D_explicit(N_elements,q,E,nu,kappa,h,b,L):
    '''
    Computes the displacements and rotations of a cantilever beam with 6 degrees
    of freedom.
    
    Note degrees of freedom:
        
    z
    X  -------->  x
    |
    |
    |
    y
    
    q[e,:] = [qx,qy,qz,mx,my,mz]
    
    where:
    q are pressures in the x, y and z-directions
    m are moments about the axis. mx is moment about x, mz is moment about y, 
    my is moment about z. See Figure X in report for more details.
    '''
    
    #Location matrix
    LM = np.zeros((12, N_elements), dtype=np.int)
    for e in range(N_elements):
        LM[0, e] = (e*6)
        LM[1, e] = (e*6)+1
        LM[2, e] = (e*6)+2
    
        LM[3, e] = (e*6)+3
        LM[4, e] = (e*6)+4
        LM[5, e] = (e*6)+5
        
        LM[6, e] = (e*6)+6
        LM[7, e] = (e*6)+7
        LM[8, e] = (e*6)+8
    
        LM[9, e] = (e*6)+9        
        LM[10, e] = (e*6)+10
        LM[11, e] = (e*6)+11
    G =  E/(2*(1+nu))
    
    nodes, dx = np.linspace(0, L, N_elements+1, retstep=True)
    
    I = [(t**3)*b[i]/12 for i,t in enumerate(h)]
    J = b[0]*h[0]**3*((16/3)-3.36*(h[0]/b[0])*(1-(h[0]**4/(12*b[0]**4))))
    I_zz = [(t)*(b[i]**3)/12 for i,t in enumerate(h)]
    
    I_yy = [(t)*(L**3)/12 for i,t in enumerate(h)]
    
    A = [(t)*(b[i]) for i,t in enumerate(h)] 
    phi_z = (12*E*I_zz[0])/(kappa*G*A[0]*L**2)
    phi_y = (12*E*I_yy[0])/(kappa*G*A[0]*L**2)
    
    phi_y_bar = 1/(1+ phi_y)
    phi_z_bar = 1/(1 + phi_z)
    # explcitly define each element stiffness matrix term
    k_z_const = (E*I_zz[0])/((1+phi_z)*dx**2)
    k_y_const = (E*I_yy[0])/((1+phi_z)*dx**2)
    
    k_11_z = k_z_const*12
    k_12_z = k_z_const*6*dx
    k_13_z = k_z_const*-12
    k_14_z = k_z_const*6*dx
    k_22_z = k_z_const*(4+phi_y)*dx**2
    k_23_z = k_z_const*-6*dx
    k_24_z = k_z_const*(2-phi_y)*dx**2
    k_33_z = k_z_const*12
    k_34_z = k_z_const*-6*dx
    k_44_z = k_z_const*(4+phi_y)*dx**2
    
    
    k_11_y = k_y_const*12
    k_12_y = k_y_const*-6*dx
    k_13_y = k_y_const*-12
    k_14_y = k_y_const*-6*dx
    k_22_y = k_y_const*(4+phi_z)*dx**2
    k_23_y = k_y_const*6*dx
    k_24_y = k_y_const*(2-phi_z)*dx**2
    k_33_y = k_y_const*12
    k_34_y = k_y_const*6*dx
    k_44_y = k_y_const*(4+phi_z)*dx**2
    
    
    ke = np.zeros((12,12))
    
    ke[0,0] = (E*A[0])/dx
    ke[0,6] = -ke[0,0]
    ke[1,1] = k_11_z
    ke[1,5] = k_12_z
    ke[1,7] = k_13_z
    ke[1,11] = k_14_z
    ke[2,2] = k_11_y
    ke[2,4] = k_12_y
    ke[2,8] = k_13_y
    ke[2,10] = k_14_y
    ke[3,3] = G*J/dx
    ke[3,9] = -ke[3,3]
    ke[4,2] = k_12_y
    ke[4,4] = k_22_y
    ke[4,8] = k_23_y
    ke[4,10] = k_24_y
    ke[5,1] = k_12_z
    ke[5,5] = k_22_z
    ke[5,7] = k_23_z
    ke[5,11] = k_24_z
    
    ke[6,0] = -(E*A[0])/dx
    ke[6,6] = -ke[6,0]
    ke[7,1] = k_13_z
    ke[7,5] = k_23_z
    ke[7,7] = k_33_z
    ke[7,11] = k_34_z
    ke[8,2] = k_13_y
    ke[8,4] = k_23_y
    ke[8,8] = k_33_y
    ke[8,10] = k_34_y
    ke[9,3] = -G*J/dx
    ke[9,9] = -ke[9,3]
    ke[10,2] = k_14_y
    ke[10,4] = k_24_y
    ke[10,8] = k_34_y
    ke[10,10] = k_44_y
    ke[11,1] = k_14_z
    ke[11,5] = k_24_z
    ke[11,7] = k_34_z
    ke[11,11] = k_44_z
    fe = np.zeros(12)
    K = np.zeros((6*(N_elements+1),6*(N_elements+1)))
    F = np.zeros((6*(N_elements+1)))
    # define each force vector term
    for e in range(N_elements):
    
        fe[0] = -q[e,0]*0.5*dx
        fe[1] = -0.5*q[e,1]*dx +q[e,5]
        fe[2] = -0.5*q[e,2]*dx +q[e,4]
        
        fe[3] = -q[e,3]*0.5*dx
        fe[4] = q[e,2]/12*dx**2 - 0.5*(phi_z/(1+phi_z))*q[e,4]*dx
        fe[5] = -q[e,1]/12*dx**2 - 0.5*(phi_y/(1+phi_y))*q[e,5]*dx
        
        fe[6] = -q[e,0]*0.5*dx
        fe[7] = -0.5*q[e,1]*dx -q[e,5]
        fe[8] = -0.5*q[e,2]*dx -q[e,4]
        
        fe[9] = -q[e,3]*0.5*dx
        fe[10] = -q[e,2]/12*dx**2 - 0.5*(phi_z/(1+phi_z))*q[e,4]*dx
        fe[11] = q[e,1]/12*dx**2 - 0.5*(phi_y/(1+phi_y))*q[e,5]*dx
        
        
        # populate global matrices
        for i in range(12):
            for j in range(12):
                K[LM[i,e],LM[j,e]] += ke[i,j]
                
                
            F[LM[i,e]] += fe[i]
    # fix nodes at root to make cantilever beam
    fixedNodeU = [0]
    fixedNodeV = [1]
    fixedNodeW = [2]
    
    fixedNodePhi = [3]
    fixedNodePsi = [4]
    fixedNodeTheta = [5]
    
    # define active degrees of freedom
    free_dof = np.delete(np.arange(0,6*(N_elements+1)),[fixedNodeW,fixedNodeU,fixedNodeV,fixedNodeTheta,fixedNodePhi,fixedNodePsi])
    active_nodes_A,active_nodes_B = np.meshgrid(free_dof,free_dof)
    # solve matrices to obtain displacement
    U = np.linalg.solve(K[active_nodes_A,active_nodes_B],F[free_dof])
    '''
    plt.figure(1)
    plt.plot(U[0::6]) #u
    plt.figure(2)
    plt.plot(U[1::6]) #v
    plt.figure(3)
    plt.plot(U[2::6]) #w
    plt.figure(4)
    plt.plot(U[3::6]) # phi
    plt.figure(5)
    plt.plot(U[4::6]) # theta 
    plt.figure(6)
    plt.plot(U[5::6]) # psi
    '''
    
    if N_elements == 1:
        errors = 0
        for i in range(12):
            for j in range(12):
                if ke[i,j] != K[i,j]:
                    print('Error')
                    errors+=1
                    
        print(errors)
        
        for i in range(12):
            if fe[i] != F[i]:
                print('Error in F')
    return U


def test_torsion():
    N_elements = 100
    q = np.zeros((N_elements,6))
    E = 2e5
    nu = 0.3
    L = 100
    kappa = 5/6
    h = 10 * np.ones(N_elements)
    b = h
    q[:,3] = 500
    U = timo_3D_explicit(N_elements,q,E,nu,kappa,h,b,L)
    plt.figure(1)
    plt.plot(U[0::6]) #u
    plt.figure(2)
    plt.plot(U[1::6]) #v
    plt.figure(3)
    plt.plot(U[2::6]) #w
    plt.figure(4)
    plt.plot(U[3::6]) # phi
    plt.figure(5)
    plt.plot(U[4::6]) # theta 
    plt.figure(6)
    plt.plot(U[5::6]) # psi
    
    
if __name__ == '__main__':
    test_torsion()