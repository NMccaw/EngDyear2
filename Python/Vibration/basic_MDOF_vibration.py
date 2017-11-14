#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 14 15:29:35 2017

@author: nicholas
"""

import numpy as np
from scipy import linalg as LA
import matplotlib.pyplot as plt
M=np.array([[1,0,0],[0,2,0],[0,0,3.5]])
# mass matrix
K=np.array([[2.,-1,0],[-1,3,-2],[0,-2,4.5]])
# stiffness matrix
wn2, B = LA.eig(K,M)

wn = wn2**(0.5) #1D array of natural frequencies
plt.plot(B,'-o') #plot modal shapes

print(B)
print(wn)