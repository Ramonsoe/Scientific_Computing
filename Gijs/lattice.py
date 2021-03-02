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
        self.walker = [0, np.random.choice(self.N)]
    
    def initial_lattice(self, y, x, y_len, x_len):
        self.lattice[0, :] = 1
        left_border = x
        right_border = x + x_len
        top_border = y
        bottom_border = y + y_len
        
        self.objects[top_border:bottom_border, left_border:right_border] = 1
        
    def SOR(self, omega):
        m = np.copy(self.lattice)
        for i in range(1, self.N-1):
            # Left boundary:
            m[i, 0]= (omega/4) * (m[i+1, 0] + m[i-1, 0] + m[i, 1] + m[i, self.N-1]) + (1 - omega) * m[i, 0]   
            # self.MC_lattice()
            # Middle part:
            for j in range(1, self.N-1):
                if (self.objects[i,j] != 1):
                    m[i, j] = (omega/4) * (m[i+1, j] + m[i-1, j] + m[i, j+1] + m[i, j-1]) + (1 - omega) * m[i, j]  
                    # self.MC_lattice()
            
            # Right boundary:
            m[i, self.N-1] = (omega/4) * (m[i+1, self.N-1] + m[i-1, self.N-1] + m[i, self.N-2] + m[i, 0]) + (1 - omega) * m[i, self.N-1]     
            # self.MC_lattice()
            
        self.delta = np.max(np.subtract(m, self.lattice))
        self.lattice = m
        
        
    def print_lattice(self, bool):
        fig, ax = plt.subplots()
        self.lattice[self.lattice == 0] = np.nan
        if (bool == True):
            ax.imshow(self.lattice, cmap="rainbow" , interpolation='nearest', extent=[0,1,0,1], aspect='auto')
        else:
            ax.imshow(self.objects, cmap="rainbow" , interpolation='nearest', extent=[0,1,0,1], aspect='auto')
        ax.set_aspect(1)
        
        plt.xlabel('x')
        plt.ylabel('y')
        plt.show()
        
    def check_neighbours(self, i,j):
        """
        This function checks the surrounding coordinates of a given coordinate set (i,j)
        If an instance of an object can be found in the lattice surrounding these
        current coordinates, True is returned
        
        For the border cases, special if-statements are implemented.
        
        If no object is found in the vicinity, False is returned. 
        """
        
        if (i == 0):
            if (self.objects[i+1, j] == 1 or self.objects[i, j+1] == 1 or self.objects[i, j-1] == 1):
                return True
        elif (i == self.N-1):
            if (self.objects[i-1, j] == 1 or self.objects[i, j+1] == 1 or self.objects[i, j-1] == 1):
                return True
        elif (j == 0):
            if (self.objects[i+1, j] == 1 or self.objects[i-1, j] == 1 or self.objects[i, j+1] == 1):
                return True
        elif (j == self.N-1):
            if (self.objects[i+1, j] == 1 or self.objects[i-1, j] == 1 or self.objects[i, j-1] == 1):
                return True
        elif (self.objects[i+1, j] == 1 or self.objects[i-1, j] == 1 or self.objects[i, j+1] == 1 or self.objects[i, j-1] == 1):
                return True
        return False
            
            
            
    
    def grow_object(self, eta):
        p = np.zeros((self.N, self.N))
        for i in range(1, self.N-1):
            for j in range(0, self.N):
                if (self.lattice[i,j] < 0):
                    self.lattice[i,j] = 0.01
                if (self.lattice[i,j] < 0):
                    print("Huh, hoe dan")
                numerator = self.lattice[i,j]**eta
                denominator = np.sum(self.lattice[1: self.N-1, :])**eta
                if (denominator == np.nan):
                    print("tering he")
                if (denominator != np.nan):
                    if (denominator > 0 and numerator >= 0):
                        p[i,j] = numerator/denominator
                    else:
                        p[i,j] = 0
                else:
                    p[i,j] = 0
                    
                c = np.random.rand(1)
                if (self.check_neighbours(i,j) == True): 
                    # print("We came across object, current coordinates = ", i, ",", j)
                    if(p[i,j] > c):
                        self.objects[i,j] = 1
                        self.lattice[i,j] = 0
                        # print("Object has been extended")
    
    

    def initialize_walker(self):
        """
        This function (re)initializes the random walker to the top of the lattice.
        """
        
        start = np.random.choice(self.N)
        self.walker = [0, start]
        
                    
    def check_boundary(self):
        """
        The random walker is walking through the lattice. At each step, some
        conditions need to be satisfied: we can't go lower than i=0, can't go
        higher than i=N-1. Furthermore, the periodic boundaries are 
        implemented: if j < 0, then j = N-2, and vice versa (j > N-1 -> j=1)
        """
        
        # Top and Bottom checks:
        if (self.walker[0] < 0 or self.walker[0] >= self.N-1):
            # The walker has exited the lattice. A new walker is created.
            self.initialize_walker()
        
        # Periodic boundaries:
        if (self.walker[1] < 0):
            self.walker[1] = self.N-2
            
        if (self.walker[1] >= self.N-1):
            self.walker[1] = 1
    
    
    
    def MC_lattice(self):
        """
        This method picks a random direction in which the random walker will
        move. After moving, check the boundary conditions and check whether
        an object is next to the random walker. 
            If the latter holds, the random walker becomes an object,
        and a new walker is created.
        """
        up_or_down = np.random.choice([-1, 1]) #, p=[0.01,0.99])
        direction = np.random.choice([0, 1]) #, p=[0.99, 0.01])
               
        self.walker[direction] = self.walker[direction] + up_or_down
        
        # Check for boundary conditions:
        self.check_boundary()
        
        i = self.walker[0]
        j = self.walker[1]
        # If we are next to the object, the random walker should be added to it
        if (self.check_neighbours(i,j) == True):
            # print("Dat noemen wij een bingo")
            self.objects[i,j] = 1
            self.lattice[i,j] = 0
            self.initialize_walker()
            
        
        
        