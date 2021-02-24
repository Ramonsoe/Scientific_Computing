import numpy as np
import matplotlib.pyplot as plt
import copy






def main():
    # Initialize lattice object:
    l = Lattice(100)
    
    # Initializing ojbect in the lattice (bacteria) and empty lattice:
    x = np.int(0.99*l.N)
    y = np.int(0.49*l.N)
    
    l.initial_lattice(x,y, 1,1)             # Top row = 1, object in lattice
    
    
    epsilon = 10**-5
    num_iterations= 0
    while(l.delta > epsilon):
        l.SOR(omega = 1.5) 
        num_iterations+=1
        
    l.print_lattice()
    print("Number of iterations needed to get here = ", num_iterations)
    
    
    
if __name__ == "__main__":
    main()