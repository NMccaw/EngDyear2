#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 24 12:07:39 2017

@author: nicholas
"""
import numpy as np
from matplotlib import pyplot as plt

def shape_func(xi):
    '''
    Defines the shape functions used in FEA
    
    Parameters
    ----------
    xi: float
        local coordinate of element
        
    Returns
    -------
    N: array of float
        shape function
    '''
    N = np.zeros(2)
    N[0] = 0.5*(1-xi)
    N[1]= 0.5*(1+xi)
    return N
    
def shape_func_dev(xi):
    '''
    Defines the derivative of the shape function w.r.t the local coordinate.
    
    Parameters
    ----------
    xi: float
        local coordinate 
    
    Returns
    -------
    dNdxi: array of float
            derivative of shape function w.r.t xi
    '''
    dNdxi = np.zeros(2)
    dNdxi[0] = -0.5
    dNdxi[1] = 0.5
    return dNdxi
    
def finite_element_1d(N_elements,q,L,E,nu,kappa,h,b):
    '''
    Main FEA program. Computes deflection and rotation of timoshenko beam using
    finite elements. Computes the stiffness matrix and force vector and populates
    appropriately. Computes deflection and rotation based on cantilever beam 
    boundary conditions.
    
    Parameters
    ----------
    N_elements: int
        number of finite elements
    q: array of floats
        load acting on each element
    L: float
        length of beam
    E: float
        modulus of elasticity(Youngs Modulus)
    nu: float
        poissions ratio
    kappa: float
        curvature constant
    h: float
        beam thickness
    b: float
        beam bredth
    Returns
    -------
    U: array of floats
        deflections and rotations of beams at nodes
    free_dof: array of floats
        array of the nodes which are free to deflect and rotate.
    '''
    G =  E/(2*(1+nu))
    
    I = [(t**3)*b[i]/12 for i,t in enumerate(h)]
    
     
    EI = [E*i for i in I]

    
    # Even grid
    nodes, dx = np.linspace(0, L, N_elements+1, retstep=True)
    M = np.zeros_like(q)
    M[-1] = 1000
    # Location matrix
    LM = np.zeros((2, N_elements), dtype=np.int)
    for e in range(N_elements):
        LM[0, e] = e
        LM[1, e] = e+1
    
    # Global stiffness matrix and force vector
    K = np.zeros((2*(N_elements+1), 2*(N_elements+1)))
    
    F = np.zeros((2*(N_elements+1),))
    dxdxi = dx/2
    integration_points = [-np.sqrt(1/3),np.sqrt(1/3)]
    # Loop over elements
    for e in range(N_elements):
        for xi in integration_points:
            N = shape_func(xi)
            dNdxi = shape_func_dev(xi)
            dNdx = dNdxi * (1/dxdxi)
            B = np.zeros((2,4))
            B[0,2:] = dNdx
            A,C = np.meshgrid([LM[:,e], LM[:,e] + N_elements+1],[LM[:,e], LM[:,e] + N_elements+1])
            # deflection stiffness in z-direction
            K[A,C] += np.dot(np.transpose(B),B) * EI[e] * dxdxi
            
            # force in z-direction
            F[LM[:,e]] += N * q[e] * dxdxi
    
            
            B = np.zeros((2,4))
            B[1,:2] = dNdx
            B[1,2:] = -N
            # stiffness in rotation
            K[A,C] += np.dot(np.transpose(B),B) * kappa * h[e] * b[e] * G * dxdxi
            
            #F[LM[:,e]] += N * M[e] * dxdxi
            
    #Boundary Conditions
    fixedNodeW = [0]
    fixedNodeTheta = [0+N_elements+1]
    
    # Boundary conditions- fixed at 0.25c
    #fixedNodeW = [int(N_elements/4)]
    #fixedNodeTheta = [int(N_elements/4) + N_elements+1]
    print(K)
    free_dof = np.delete(np.arange(0,2*(N_elements+1)),[fixedNodeW,fixedNodeTheta])
    
    
    active_nodes_A,active_nodes_B = np.meshgrid(free_dof,free_dof)
    
    U = np.linalg.solve(K[active_nodes_A,active_nodes_B],F[free_dof])
    
    
    
    return U,free_dof

def test_moment():
    
    N_elements = 2
    q = np.ones(N_elements) * 0 
    
    L = 1
    E = 500e6
    kappa = 5/6
    nu = 0.3 
    h = np.ones_like(q) * 0.2
    b = np.ones_like(h)  
    
    U,free_dof = finite_element_1d(N_elements,q,L,E,nu,kappa,h,b)
    
    plt.plot(free_dof[:N_elements],U[:N_elements])
    
if __name__ == '__main__':
    test_moment()
    