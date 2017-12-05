#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 27 16:15:28 2017

@author: Nicholas
"""
import numpy as np
from matplotlib import pyplot as plt
from scipy import linalg as LA

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
    
def finite_element_1d(N_elements,q,L,E,nu,kappa,h,b,rho):
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
    
    # Location matrix
    LM = np.zeros((2, N_elements), dtype=np.int)
    for e in range(N_elements):
        LM[0, e] = e
        LM[1, e] = e+1
    
    # Global stiffness matrix and force vector
    K = np.zeros((2*(N_elements+1), 2*(N_elements+1)))
    M = np.zeros_like(K)
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
            
            K[A,C] += np.dot(np.transpose(B),B) * EI[e] * dxdxi
            M[A,C] += np.dot(np.transpose(N),N) * rho* b[e]* h[e] * dxdxi
            
            F[LM[:,e]] += N * q[e] * dxdxi
    
            
            B = np.zeros((2,4))
            B[1,:2] = dNdx
            B[1,2:] = -N
            K[A,C] += np.dot(np.transpose(B),B) * kappa * h[e] * b[e] * G * dxdxi
    
    #Boundary Conditions
    fixedNodeW = [0]
    fixedNodeTheta = [0+N_elements+1]
    
    # Boundary conditions- fixed at 0.25c
    #fixedNodeW = [int(N_elements/4)]
    #fixedNodeTheta = [int(N_elements/4) + N_elements+1]
    
    free_dof = np.delete(np.arange(0,2*(N_elements+1)),[fixedNodeW,fixedNodeTheta])
    
    
    active_nodes_A,active_nodes_B = np.meshgrid(free_dof,free_dof)
    
    U = np.linalg.solve(K[active_nodes_A,active_nodes_B],F[free_dof])
    
    
    
    return U,free_dof,M,K
  
N_elements = 100
# q is force per unit length 
q = np.zeros(N_elements)
L = 10
elements_per_unit_length = N_elements//L
q[:] = -1000



def BEMT_FEA_mesh_interp(N_elements_FEA,N_elements_BEMT,Cl,L):
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
    for e in range(N_elements_BEMT):
        LM_BEMT[0,e] = e*L/N_elements_BEMT
        LM_BEMT[1,e] = (e+1)*L/N_elements_BEMT

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
                
                
    
def test_BEMT_FEA_mesh_interp():
    '''
    tests the interpolation between meshes occurs correctly
    '''
    N_BEMT = 4
    N_FEA = 10
    Cl = np.linspace(1,10,N_BEMT)
    L= 10
    q = BEMT_FEA_mesh_interp(N_FEA,N_BEMT,Cl,L)
    assert np.all(q == np.array([1,1,2.5,4,4,7,7,8.5,10,10]))

    
    
    
    
def test_timoshenko_uniform_load():
    '''
    verification that timoshenko beam theory code produces same results as other 
    codes
    '''
    E = 5e6
    nu = 0.3
    N_elements = 1000
    L = 10
    kappa = 5/6
    h = 2
    # q is pressure per element area not force
    q = np.zeros(N_elements)
    q[:] = -1000
    U,free_dof = finite_element_1d(N_elements,q,L,E,nu,kappa,h)
    plt.plot(free_dof[:N_elements],U[:N_elements])
    plt.savefig('timo_test.pdf')
    assert np.allclose(U[N_elements-1],-0.3906)
    
def test_timoshenko_uniform_load_experiment():
    '''
    test the code produces the same results as experimental data.
    '''
    E = 194.3e9
    nu = 0.3
    N_elements = 100
    L = 0.4
    kappa = 5/6
    h = 0.0004 *np.ones(N_elements)
    b = 0.025 * np.ones(N_elements)
    # q is pressure per element area not force
    q = np.zeros(N_elements)
    q[:] = -0.758
    U,free_dof = finite_element_1d(N_elements,q,L,E,nu,kappa,h,b)
    plt.plot(free_dof[:N_elements],U[:N_elements])
    print(U[N_elements-1])
    #assert np.allclose(U[N_elements-1],-0.3906)
    
def test_timoshenko_analytical(N_elements):
    '''
    verification that timoshenko beam theory code produces same results as other 
    codes
    '''
    E = 2e5
    nu = 0.3
    
    L = 100
    kappa = 5/6
    h = 10 * np.ones(N_elements)
    b = h
    # q is pressure per element area not force
    q = np.zeros(N_elements)
    q[:] = -500
    U,free_dof = finite_element_1d(N_elements,q,L,E,nu,kappa,h,b)
    plt.plot(free_dof[:N_elements],U[:N_elements])
    
    plt.savefig('timo_test.pdf')
    plt.figure(2)
    plt.plot(free_dof[:N_elements],U[N_elements:])
    return U[N_elements-1]

def vibration_modes(K,M):
    print('computing modes...')
    wn2, B = LA.eig(K,M)
    print('Modes computed!')
    wn = wn2**(0.5) #1D array of natural frequencies
    print('Computed natural frequency = {}'.format(wn[0]))
    
    plt.plot(B[0::6,0:4],'-o') #plot modal shapes
def test_mass_matrix():
    N_elements = 100
    E = 2.1e11
    G = 8.1e10 
    rho = 7860 
    L = 1 
    b = 0.02 * np.ones(N_elements)
    h = 0.08 * np.ones(N_elements) #0.03
    q = np.zeros(N_elements)
    #q[:,1] = 500
    nu = 0.3
    kappa = 5/6   
    U,free_dof,M,K = finite_element_1d(N_elements,q,L,E,nu,kappa,h,b,rho)
    print(K)
    print(M)
    vibration_modes(K,M)
    
    
    
if __name__ == '__main__':
    
    '''
    q_list = np.zeros(9)
    N_elements_list = np.zeros_like(q_list)
    error = np.zeros_like(q_list)
    for i,N_elements in enumerate(range(6,13)):
       print(i) 
       q_list[i] = test_timoshenko_analytical(2**N_elements)
       error[i] = abs(q_list[i] - q_list[i-1])
       print('completed n_elements = {}'.format(2**N_elements))
       N_elements_list[i] = 2**N_elements
    plt.figure(2)
    plt.loglog(N_elements_list[1:],error[1:])
    plt.title('Error Convergence of FEA')
    plt.xlabel('Number of elements')
    plt.ylabel('Error')
    plt.savefig('FEA_convergence.eps')
    '''
    #test_timoshenko_uniform_load_experiment()
    test_mass_matrix()