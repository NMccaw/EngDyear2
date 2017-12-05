#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 25 11:32:18 2017

@author: nicholas
"""
import numpy as np
from matplotlib import pyplot as plt

def shape_func_3D(E,I_yy,G,L,I_zz,kappa,A,xi):
    N  = np.zeros((6,12))
    phi_z = (12*E*I_zz)/(kappa*G*A*L**2)
    phi_y = (12*E*I_yy)/(kappa*G*A*L**2)
    
    phi_y_bar = 1/(1- phi_y)
    phi_z_bar = 1/(1 - phi_z)
    
    N[0,0] = (1-xi) #N1    
    N[0,6] = xi #N2
    N[1,1] = phi_y_bar*(2*xi**3 -3*xi**2 + phi_y*xi+1-phi_y) # Hv1
    N[1,7] = phi_y_bar*(-2*xi**3 +3*xi**2 - phi_y*xi) # Hv2
    N[1,4] = L*phi_y_bar*(xi**3 +(0.5*phi_y-2)*xi**2 +(1-0.5*phi_y)*xi) #Ht1
    N[1,10] = L*phi_y_bar*(xi**3 -(0.5*phi_y+1)*xi**2 +(0.5*phi_y)*xi) #Ht2
    N[2,2] = phi_z_bar*(2*xi**3 -3*xi**2 + phi_z*xi+1-phi_z) # Hw1
    N[2,8] = phi_z_bar*(-2*xi**3 +3*xi**2 - phi_z*xi) # Hw2
    N[2,5] = L*phi_z_bar*(xi**3 + (0.5*phi_z -2)*xi**2 + (1- 0.5*phi_z)*xi) #Hpsi1
    N[2,11] = L*phi_z_bar*(xi**3 - (1+0.5*phi_z)*xi**2 + (0.5*phi_z)*xi) #Hpsi2
    N[3,3] = 1 - xi # N1
    N[3,9] = xi # N2
    N[4,1] = 6*phi_y_bar*(-xi +xi**2)/L #Gv1
    N[4,7] = 6*phi_y_bar*(xi-xi**2)/L# Gv2
    N[4,4] = phi_y_bar*(3*xi**2 + (phi_y-4)*xi + 1 - phi_y) #Gt1
    N[4,10] = phi_y_bar*(3*xi**2 - (phi_y +2)*xi) #Gt2
    N[5,2] = 6*phi_z_bar*(-xi +xi**2)/L #Gw1
    N[5,8] = 6*phi_z_bar*(xi-xi**2)/L #Gw2 
    N[5,5] = phi_z_bar*(3*xi**2 + (phi_z-4)*xi + 1 - phi_y) #Gpsi1
    N[5,11] = phi_z_bar*(3*xi**2 - (phi_y +2)*xi) #Gpsi2
    
    
    
    return N

def shape_func_derive_3D(E,I_yy,G,L,I_zz,kappa,A,xi):
    dN  = np.zeros((6,12))
    phi_z = (12*E*I_zz)/(kappa*G*A*(L**2))
    phi_y = (12*E*I_yy)/(kappa*G*A*(L**2))
    
    phi_y_bar = 1/(1 - phi_y)
    phi_z_bar = 1/(1 - phi_z)
    
    dN[0,0] = -1 #N1    
    dN[0,6] = 1 #N2
    dN[1,1] = phi_y_bar*(6*xi**2 -6*xi + phi_y)
    dN[1,7] = phi_y_bar*(-6*xi**2 +6*xi - phi_y)
    dN[1,4] = L*phi_y_bar*(3*xi**2 + 2*(0.5*phi_y-2)*xi +(1-0.5*phi_y))
    dN[1,10] = L*phi_y_bar*(3*xi**2 - 2*(0.5*phi_y+1)*xi +(0.5*phi_y))
    dN[2,2] = phi_z_bar*(6*xi**2 -6*xi + phi_z)
    dN[2,8] = phi_z_bar*(-6*xi**2 +6*xi - phi_z)
    dN[2,5] = L*phi_z_bar*(3*xi**2 + 2*(0.5*phi_z -2)*xi + (1- 0.5*phi_z))
    dN[2,11] = L*phi_z_bar*(3*xi**2 - 2*(1+0.5*phi_z)*xi + (0.5*phi_z))
    dN[3,3] = -1
    dN[3,9] = 1
    dN[4,1] = 6*phi_y_bar*(-1 +2*xi)/L #Gv1
    dN[4,7] = 6*phi_y_bar*(1-2*xi)/L # Gv2
    dN[4,4] = phi_y_bar*(6*xi + (phi_y-4))
    dN[4,10] = phi_y_bar*(6*xi - (phi_y +2))
    dN[5,2] = 6*phi_z_bar*(-1 + 2*xi)/L #Gw1
    dN[5,8] = 6*phi_z_bar*(1-2*xi)/L #Gw2 
    dN[5,5] = phi_z_bar*(6*xi + (phi_z-4))
    dN[5,11] = phi_z_bar*(6*xi - (phi_y +2))
    
    return dN


def shape_func_2nd_derive_3D(E,I_yy,G,L,I_zz,kappa,A,xi):
    ddN  = np.zeros((6,12))
    phi_z = (12*E*I_zz)/(kappa*G*A*(L**2))*0
    phi_y = (12*E*I_yy)/(kappa*G*A*(L**2))*0
    
    phi_y_bar = 1/(1 - phi_y)
    phi_z_bar = 1/(1 - phi_z)
    
    ddN[0,0] = 0 #N1    
    ddN[0,6] = 0 #N2
    ddN[1,1] = phi_y_bar*(12*xi -6)
    ddN[1,7] = phi_y_bar*(-12*xi +6)
    ddN[1,4] = L*phi_y_bar*(6*xi + 2*(0.5*phi_y-2))
    ddN[1,10] = L*phi_y_bar*(6*xi - 2*(0.5*phi_y+1))
    ddN[2,2] = phi_z_bar*(12*xi -6)
    ddN[2,8] = phi_z_bar*(-12*xi +6)
    ddN[2,5] = L*phi_z_bar*(6*xi + 2*(0.5*phi_z -2))
    ddN[2,11] = L*phi_z_bar*(6*xi - 2*(1+0.5*phi_z))
    ddN[3,3] = 0
    ddN[3,9] = 0
    ddN[4,1] = 6*phi_y_bar*(2)/L #Gv1
    ddN[4,7] = 6*phi_y_bar*(2)/L # Gv2
    ddN[4,4] = phi_y_bar*(6)
    ddN[4,10] = phi_y_bar*(6)
    ddN[5,2] = 6*phi_z_bar*(2)/L #Gw1
    ddN[5,8] = 6*phi_z_bar*(2)/L #Gw2 
    ddN[5,5] = phi_z_bar*(6)
    ddN[5,11] = phi_z_bar*(6)
    
    return ddN


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
    
    I_zz = [(t)*(b[i]**3)/12 for i,t in enumerate(h)]
    
    I_yy = [(t)*(L**3)/12 for i,t in enumerate(h)]
    
    A = [(t)*(b[i]) for i,t in enumerate(h)] 
    
    Ay = [(t)*(L)*0 for i,t in enumerate(h)] 
    
    Az = [(L)*(b[i])*0 for i,t in enumerate(h)] 
    
    
    J = 0.2e5 *0
    
    EI = [E*i for i in I]
    
    
    
    # Even grid
    nodes, dx = np.linspace(0, L, N_elements+1, retstep=True)
    
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
        
 # Location matrix
#    LM = np.zeros((2, N_elements), dtype=np.int)
#    for e in range(N_elements):
#        LM[0,e] = e
#        LM[1,e] = e+1
    # Global stiffness matrix and force vector
    K = np.zeros((6*(N_elements+1), 6*(N_elements+1)))
    ke = np.zeros((12, 12))
    
    F = np.zeros(6*(N_elements+1))
    fe = np.zeros((12,))
    dxdxi = dx
    integration_points = [-np.sqrt(1/3),np.sqrt(1/3)]
    
    D = np.zeros((6,6))
    D[0,0] = E*A[0]
    D[0,1] = E *Ay[0]
    D[0,2] = E*Az[0]
    D[1,0] = E*Ay[0]
    D[1,1] = E*I_yy[0]
    D[1,2] = E*J
    D[2,0] = E*Az[0]
    D[2,1] = E*J
    D[2,2] = E*I_zz[0]
    D[3,3] = kappa*G*A[0]
    D[3,5] = -kappa*G*Az[0]
    D[4,4] = kappa*G*A[0]
    D[4,5] = kappa*G*Ay[0]
    D[5,3] = -kappa*G*Az[0]
    D[5,4] = D[4,5]
    D[5,5] = kappa*G*(I_yy[0]+I_zz[0])
    # Loop over elements
    #B = np.zeros((12,12))
    for e in range(N_elements):
        for xi in integration_points:
            
            N = shape_func_3D(E,I_yy[e],G,dx,I_zz[e],kappa,A[e],xi)
            
            dNdxi = shape_func_derive_3D(E,I_yy[e],G,dx,I_zz[e],kappa,A[e],xi)
            
            ddNdxi = shape_func_2nd_derive_3D(E,I_yy[e],G,dx,I_zz[e],kappa,A[e],xi)
            
            ddNdx = ddNdxi* (1/dxdxi)
            
            dNdx = dNdxi * (1/dxdxi)
            
            
            B = dNdx
            
            A_2,C = np.meshgrid([LM[:,e]],[LM[:,e]])
            
            K[A_2,C] += np.dot(np.dot(np.transpose(B),D),B) * dxdxi/2
            #K[A_2,C] += np.dot(np.transpose(B),B)*E*I_yy[0] * dxdxi
            ke = np.dot(np.dot(np.transpose(B),D),B) * dxdxi
            F[LM[:,e]] += np.dot(q[e,:],N) * dxdxi
            
           
        '''
        #populate global stiffness matrix
            for i in range(12):
                for j in range(12):
                    K[LM[i,e],LM[j,e]] += ke[i,j]
                    
                    
                F[LM[i,e]] += fe[i]
        '''
    #print(F[1::6])
    #F[1::6] = 1000
    #Boundary Conditions
    fixedNodeU = [0]
    fixedNodeV = [1]
    fixedNodeW = [2]
    
    fixedNodePhi = [3]
    fixedNodePsi = [4]
    fixedNodeTheta = [5]
    
    # Boundary conditions- fixed at 0.25c
    #fixedNodeW = [int(N_elements/4)]
    #fixedNodeTheta = [int(N_elements/4) + N_elements+1]
    free_dof = np.delete(np.arange(0,6*(N_elements+1)),[fixedNodeW,fixedNodeU,fixedNodeV,fixedNodeTheta,fixedNodePhi,fixedNodePsi])
    active_nodes_A,active_nodes_B = np.meshgrid(free_dof,free_dof)
    print(ke[2,:6])
    print(ke[2,6:])
    ''' 
    for i in range(12):
        for j in range(12):
            print('('+str(i)+','+str(j)+')')
            print(ke[i,j])
    print(12*kappa*G*A[0]*E*I_yy[0]*((12*E*I_yy[0])+(kappa*G*A[0]*(L**2)))/(L*(((12*E*I_yy[0]) - (kappa*G*A[0]*(L**2)))**2))) # K[1,1]      
    '''  
    '''      
    print(K[1,7])
    print(12*kappa*G*A[0]*E*I_yy[0]*(12*E*I_yy[0]+(kappa*G*A[0]*(L**2)))/(L*((12*E*I_yy[0] - kappa*G*A[0]*(L**2))**2))) # K[1,1]
    
    print(6*kappa*G*A[0]*E*I_yy[0]*(12*E*I_yy[0]+kappa*G*A[0]*L**2)/((12*E*I_yy[0] - kappa*G*A[0]*L**2)**2)) # K[1,5]
    print(K[2,2])
    print(12*kappa*G*A[0]*E*I_zz[0]*(12*E*I_zz[0]+(kappa*G*A[0]*(L**2)))/(L*((12*E*I_zz[0] - kappa*G*A[0]*(L**2))**2))) # k[2,2]
    
    assert K[1,5] == 6*kappa*G*A[0]*E*I_yy[0]*(12*E*I_yy[0]+kappa*G*A[0]*L**2)/(L*(12*E*I_yy[0] - kappa*G*A[0]*L**2)**2)
    assert K[0,0] == K[6,6] == -K[0,6] 
    assert K[1,1] == -K[1,7]
    assert K[1,5] == K[1,11]
    assert K[2,2] == K[2,10]
    '''
    U = np.linalg.solve(K[active_nodes_A,active_nodes_B],F[free_dof])
    
    
    return U,free_dof
    
if __name__ == '__main__':
    
    N_elements = 1
    q = np.zeros((N_elements,6))
    E = 2e5
    nu = 0.3
    L = 100
    kappa = 5/6
    h = 10 * np.ones(N_elements)
    b = h
    q[:,1] = 500
    
    
    
    U,free_dof = finite_element_1d(N_elements,q,L,E,nu,kappa,h,b)
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
    
    print(U)