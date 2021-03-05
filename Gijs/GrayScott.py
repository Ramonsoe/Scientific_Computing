# -*- coding: utf-8 -*-
"""
Created on Wed Mar  3 10:25:53 2021

@author: gijsv
"""

class GrayScott:
    import numpy as np
    import matplotlib.pyplot as plt
    import copy
    
    def __init__(self, N, timesteps):
        self.N = N
        self.timesteps = timesteps
        self.Du = self.Dv = self.f = self.K = 0
        self.dt = 1
        self.dx = self.dy = 1
        self.system = []
        
        self.u = 0
        self.v = 1
    
    def initialize_system(self):
        """
        This function initializes the system with uniform u = 0.5 
        across both dimensions with a small square of v = 0.25 in
        the centre. 
        """
        sys = np.zeros((self.N, self.N), dtype = np.dtype('2f'))
        for i in range(self.N):
            for j in range(self.N):
                sys[i,j] = (0.5, 0)                     # Value for u and v
                
        # Also, start with a small square of v = 0.25 in the centre:
        i = j = np.int(0.49 * self.N)
        sys[i,j] = (0.5, 0.25)
        
        for t in range(self.timesteps):
            self.system.append(sys)
        
        print("INITIAL SYSTEM FOR U:")
        self.print_system(0, 0)
        
    def safe_diffusion(self, i,j):
        """
        This function checks boundaries of the system and returns an integer
        based on what diffusion could take place
        """
        if (i == 0):
            if (j == 0):
                # Bottom left border
                return 5
            elif (j == self.N-1):
                # Bottom right border:
                return 7
            else:
                # Bottom border:
                return 1
        elif (i == self.N-1):
            if (j == 0):
                # Top left border:
                return 6
            if (j == self.N-1):
                # Top right border:
                return 8
            else:
                # Top border
                return 2
        elif (j == 0):
            # Left border
            return 3
        elif (j == self.N-1):
            # Right border:
            return 4
        else:
            # If we get here: we are just in the middle, so plain vanilla 5 point stencil
            return 0

        
    def next_step(self, i, j, k, el):
        """
        Based on our code from assignment 1, we can solve the diffusion
        equation for either u or v based on the indicator el.
        Indicator el = 0 --> solve for u
        Indicator el = 1 --> solve for v
        
        Also, checks for the boundary conditions are implemented. 
        """
        if (el == 0):
            D = self.Du
        elif (el == 1):
            D = self.Dv
        
        
        s = self.system
        if (self.safe_diffusion(i,j) == 0):
            next_timestep = ((self.dt*D)/(self.dx**2)) * (s[k][i+1,j][el] + s[k][i-1,j][el] + s[k][i,j+1][el] + s[k][i,j-1][el] - 4*s[k][i,j][el]) # + s[k][i,j][el] 
        elif (self.safe_diffusion(i,j) == 1):
            # We are against bottom border
            next_timestep = ((self.dt*D)/(self.dx**2)) * (s[k][i+1,j][el] + s[k][i,j+1][el] + s[k][i,j-1][el] - 3*s[k][i,j][el]) # + s[k][i,j][el] 
        elif (self.safe_diffusion(i,j) == 2):
            # We are against top border
            next_timestep = ((self.dt*D)/(self.dx**2)) * (s[k][i-1,j][el] + s[k][i,j+1][el] + s[k][i,j-1][el] - 3*s[k][i,j][el]) # + s[k][i,j][el] 
        elif (self.safe_diffusion(i,j) == 3):
            # We are against left border
            next_timestep = ((self.dt*D)/(self.dx**2)) * (s[k][i+1,j][el] + s[k][i-1,j][el] + s[k][i,j+1][el] - 3*s[k][i,j][el]) # + s[k][i,j][el] 
        elif (self.safe_diffusion(i,j) == 4):
            # We are against right border
            next_timestep = ((self.dt*D)/(self.dx**2)) * (s[k][i+1,j][el] + s[k][i-1,j][el] + s[k][i,j-1][el] - 3*s[k][i,j][el]) # + s[k][i,j][el] 
        elif (self.safe_diffusion(i,j) == 5):
            # We are against bottom-left border
            next_timestep = ((self.dt*D)/(self.dx**2)) * (s[k][i+1,j][el] + s[k][i,j+1][el] - 2*s[k][i,j][el]) # + s[k][i,j][el] 
        elif (self.safe_diffusion(i,j) == 6):
            # We are against top-left border
            next_timestep = ((self.dt*D)/(self.dx**2)) * (s[k][i-1,j][el] + s[k][i,j+1][el] - 2*s[k][i,j][el]) # + s[k][i,j][el] 
        elif (self.safe_diffusion(i,j) == 7):
            # We are against bottom-right border
            next_timestep = ((self.dt*D)/(self.dx**2)) * (s[k][i+1,j][el] + s[k][i,j-1][el] - 2*s[k][i,j][el]) # + s[k][i,j][el] 
        elif (self.safe_diffusion(i,j) == 8):
            # We are against top-right border
            next_timestep = ((self.dt*D)/(self.dx**2)) * (s[k][i-1,j][el] + s[k][i,j-1][el] - 2*s[k][i,j][el]) # + s[k][i,j][el] 
        
        return next_timestep
        
        
        
    def solve_gray_scott(self, Du, Dv, f, K, timestep):
        # First, set parameters of the reaction-diffusion system
        # We put them here to also allow for different parameters over time,
        # just in case.
        self.Du = Du
        self.Dv = Dv
        self.f = f
        self.K = K                                                             # Capital K to not be confused with discrete timestepping
        t = timestep
        # Then, we can iterate over the system:
        dudt = dvdt = np.zeros((self.N, self.N))
        for i in range(self.N):
            for j in range(self.N):
                dudt[i,j] = self.next_step(i,j, t, self.u) - (self.system[t][i,j][self.u] * self.system[t][i,j][self.v]**2) + self.f*(1 - self.system[t][i,j][self.u])
                # self.system[t+1][i,j][self.u] = self.system[t][i,j][self.u] + dudt
                
                dvdt[i,j] = self.next_step(i,j, t, self.v) + (self.system[t][i,j][self.u] * self.system[t][i,j][self.v]**2) - (self.f + self.K)*self.system[t][i,j][self.v]
                # self.system[t+1][i,j][self.v] = self.system[t][i,j][self.v] + dvdt
        for i in range(self.N):
            for j in range(self.N):
                self.system[t+1][i,j][0] = dudt[i,j] + self.system[t][i,j][0]
                self.system[t+1][i,j][1] = dvdt[i,j] + self.system[t][i,j][1]
     
           
    def solve_gray_scott_noise(self, Du, Dv, f, K, timestep, rand):
        """
        Solving the Gray Scott model, except that a little bit of noise 
        has been added to the system.
        """
        # First, set parameters of the reaction-diffusion system
        # We put them here to also allow for different parameters over time,
        # just in case.
        self.Du = Du
        self.Dv = Dv
        self.f = f
        self.K = K                                                             # Capital K to not be confused with discrete timestepping
        t = timestep
        # Then, we can iterate over the system:
        dudt = dvdt = np.zeros((self.N, self.N))
        for i in range(self.N):
            for j in range(self.N):
                dudt[i,j] = self.next_step(i,j, t, self.u) - (self.system[t][i,j][self.u] * self.system[t][i,j][self.v]**2) + self.f*(1 - self.system[t][i,j][self.u])
                # self.system[t+1][i,j][self.u] = self.system[t][i,j][self.u] + dudt
                
                dvdt[i,j] = self.next_step(i,j, t, self.v) + (self.system[t][i,j][self.u] * self.system[t][i,j][self.v]**2) - (self.f + self.K)*self.system[t][i,j][self.v]
                # self.system[t+1][i,j][self.v] = self.system[t][i,j][self.v] + dvdt
        for i in range(self.N):
            for j in range(self.N):
                self.system[t+1][i,j][0] = np.random.uniform(low=1-rand, high=1+rand) * dudt[i,j] + self.system[t][i,j][0]
                self.system[t+1][i,j][1] = np.random.uniform(low=1-rand, high=1+rand) * dvdt[i,j] + self.system[t][i,j][1]
    
            
    def print_system(self, el, t):
        fig, ax = plt.subplots()
        sys = np.zeros((self.N, self.N))
        for i in range(self.N):
            for j in range(self.N):
                sys[i,j] = self.system[t][i,j][el]
        ax.imshow(sys, cmap="hot", interpolation='nearest', aspect="auto") #, extent=[0,1,0,1], aspect='auto')
        ax.set_aspect(1)
        
        plt.xlabel('x')
        plt.ylabel('y')
        plt.show()
        
        
    """
    FUNCTION WITH WRONG INDEXING:
        
    def next_step(self, i, j, k, el):
        
        Based on our code from assignment 1, we can solve the diffusion
        equation for either u or v based on the indicator el.
        Indicator el = 0 --> solve for u
        Indicator el = 1 --> solve for v
        
        Also, checks for the boundary conditions are implemented. 
        
        if (el == 0):
            D = self.Du
        elif (el == 1):
            D = self.Dv
        
        
        s = self.system
        if (self.safe_diffusion(i,j) == 0):
            next_timestep = s[k][i,j[el] + ((self.dt*D)/(self.dx**2)) * (s[k][i+1,j,k][el] + s[i-1,j,k][el] + s[i,j+1,k][el] + s[i,j-1,k][el] - 4*s[i,j,k][el])
        elif (self.safe_diffusion(i,j) == 1):
            # We are against bottom border
            next_timestep = s[i,j, k][el] + ((self.dt*D)/(self.dx**2)) * (s[i+1,j,k][el] + s[i,j+1,k][el] + s[i,j-1,k][el] - 3*s[i,j,k][el])
        elif (self.safe_diffusion(i,j) == 2):
            # We are against top border
            next_timestep = s[i,j, k][el] + ((self.dt*D)/(self.dx**2)) * (s[i-1,j,k][el] + s[i,j+1,k][el] + s[i,j-1,k][el] - 3*s[i,j,k][el])
        elif (self.safe_diffusion(i,j) == 3):
            # We are against left border
            next_timestep = s[i,j, k][el] + ((self.dt*D)/(self.dx**2)) * (s[i+1,j,k][el] + s[i-1,j,k][el] + s[i,j+1,k][el] - 3*s[i,j,k][el])
        elif (self.safe_diffusion(i,j) == 4):
            # We are against right border
            next_timestep = s[i,j, k][el] + ((self.dt*D)/(self.dx**2)) * (s[i+1,j,k][el] + s[i-1,j,k][el] + s[i,j-1,k][el] - 3*s[i,j,k][el])
        elif (self.safe_diffusion(i,j) == 5):
            # We are against bottom-left border
            next_timestep = s[i,j, k][el] + ((self.dt*D)/(self.dx**2)) * (s[i+1,j,k][el] + s[i,j+1,k][el] - 2*s[i,j,k][el]) 
        elif (self.safe_diffusion(i,j) == 6):
            # We are against top-left border
            next_timestep = s[i,j, k][el] + ((self.dt*D)/(self.dx**2)) * (s[i-1,j,k][el] + s[i,j+1,k][el] - 2*s[i,j,k][el]) 
        elif (self.safe_diffusion(i,j) == 7):
            # We are against bottom-right border
            next_timestep = s[i,j, k][el] + ((self.dt*D)/(self.dx**2)) * (s[i+1,j,k][el] + s[i,j-1,k][el] - 2*s[i,j,k][el])
        elif (self.safe_diffusion(i,j) == 8):
            # We are against top-right border
            next_timestep = s[i,j, k][el] + ((self.dt*D)/(self.dx**2)) * (s[i-1,j,k][el] + s[i,j-1,k][el] - 2*s[i,j,k][el])
        
        if (next_timestep == np.nan):
            print("oeioeioei")
        return next_timestep
     """  
        
        
        
        
        
        
        