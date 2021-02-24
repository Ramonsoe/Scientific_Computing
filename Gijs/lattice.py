# -*- coding: utf-8 -*-
"""
Created on Wed Feb 24 11:44:14 2021

@author: gijsv
"""

class Lattice:
    import numpy as np
    import matplotlib.pyplot as plt
    import copy
    
    def __init__(self, N):
        self.N = N
        self.L = 1/N
        self.dx = self.dy = 1/N
        self.lattice = np.zeros((self.N, self.N))
        self.objects = np.zeros((self.N, self.N))
        self.delta = 100
        
    
    def initial_lattice(self, x, y, x_len, y_len):
        self.lattice[0, :] = 1
        self.objects[x:x+x_len, y:y+y_len] = 1
        
    def SOR(self, omega):
        m = np.copy(self.lattice)
        for i in range(1, self.N-1):
            # Left boundary:
            m[i, 0]= (omega/4) * (m[i+1, 0] + m[i-1, 0] + m[i, 1] + m[i, self.N-1]) + (1 - omega) * m[i, 0]   
        
            # Middle part:
            for j in range(1, self.N-1):
                if (self.objects[i,j] == 0):
                    m[i, j] = (omega/4) * (m[i+1, j] + m[i-1, j] + m[i, j+1] + m[i, j-1]) + (1 - omega) * m[i, j]  
            
            # Right boundary:
            m[i, self.N-1] = (omega/4) * (m[i+1, self.N-1] + m[i-1, self.N-1] + m[i, self.N-2] + m[i, 0]) + (1 - omega) * m[i, self.N-1]     
            
        self.delta = np.max(np.subtract(m, self.lattice))
        self.lattice = m
        
        
    def print_lattice(self):
        fig, ax = plt.subplots()
        
        ax.imshow(self.lattice, cmap="rainbow", interpolation='nearest', extent=[0,1,0,1], aspect='auto')
        ax.set_aspect(1)
        
        plt.xlabel('x')
        plt.ylabel('y')
        plt.show()
        
        
        