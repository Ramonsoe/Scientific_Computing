import numpy as np
import matplotlib.pyplot as plt
import copy






def main():
    # Initialize lattice object:
    l = Lattice(50)
    
    # Initializing ojbect in the lattice (bacteria) and empty lattice:
    y = np.int(0.98*l.N)
    x = np.int(0.40*l.N)
    l.initial_lattice(y,x, y_len=1, x_len=20)            # Top row = 1 -> index 0 in matrix
                                                        # Bottom row is where object of length 8 is situated.
   
    
    epsilon = 10**-5
    num_iterations= 0
    while(num_iterations < 1000000): #l.delta > epsilon):
        # l.SOR(omega = 1.3) 
        l.MC_lattice()
        num_iterations+=1
        
    l.print_lattice(False)
    print("Number of iterations needed to get here = ", num_iterations)
    print("Number of items in object = ", np.count_nonzero(l.objects))
    
    # l.MC_lattice()
    
if __name__ == "__main__":
    main()