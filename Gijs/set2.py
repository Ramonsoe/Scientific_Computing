import numpy as np
import matplotlib.pyplot as plt
import copy
import time





def main():
    # PART 3: GRAY SCOTT:
    """
    N = 50
    timesteps = 100
    
    s = GrayScott(N, timesteps)
    s.initialize_system()
    
    for i in range(timesteps-1):
        s.solve_gray_scott(Du=0.16 , Dv=0.26 , f=0.035 , K=0.06 , timestep=i)
        # s.solve_gray_scott_noise(Du=0.16 , Dv=0.08 , f=0.09 , K=0.06 , timestep=i, rand=0.3)
        if (i%10 == 0 or i > 560):
            print("System at timestep ", i)
            print("u in the corner = ", s.system[i][1,1][0])
            # s.print_system(0, i)
            print("u:")
            s.print_system(0,i)
            print("v:")
            s.print_system(1, i)
    
    print("System at timestep ", i-1)
    s.print_system(0,i)
    s.print_system(1,i)
    
    # print("Hoe ziet matrixje er eigenlijk uit nu?")
    # print(s.system[1][:,:])
    
    
        
    
    # PART 1 AND 2: LATTICE:
    """
    # Initialize lattice object:
    l = Lattice(50)
    
    # Initializing ojbect in the lattice (bacteria) and empty lattice:
    y = np.int(0.98*l.N)
    x = np.int(0.47*l.N)
    l.initial_lattice(y,x, y_len=1, x_len=6)            # Top row = 1 -> index 0 in matrix
                                                        # Bottom row is where object of length 8 is situated.
   
    
    epsilon = 10**-5
    num_iterations= 0
    start_time = time.time()
    while(num_iterations < 1000000): #l.delta > epsilon):
        # l.SOR(omega = 1)
        # l.grow_object(eta=0.5)
        l.MC_lattice_probability(eta = 0.1)                                     # Random Walker with probabilityt of sticking, eta simulation
        # l.MC_lattice()
        num_iterations+=1
    print("--- %s seconds ---" % (time.time() - start_time))
    l.print_lattice(False)                                                       # If True: print entire lattice (including SOR). If False: only prints objects
    print("Number of iterations needed to get here = ", num_iterations)
    print("Number of items in object = ", np.count_nonzero(l.objects))
    
    # l.MC_lattice()
    
if __name__ == "__main__":
    main()