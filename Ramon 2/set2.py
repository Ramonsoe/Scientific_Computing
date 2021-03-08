#!/usr/local/bin/python3

import numpy as np
import matplotlib.pyplot as plt
import copy
from lattice import Lattice
import time

def main():
    # Initialize lattice object:
    l = Lattice(50)

    # Initializing ojbect in the lattice (bacteria) and empty lattice:
    x = int(0.99*l.N)
    x = l.N  - 1
    y = int(0.49*l.N)

    l.initial_lattice(int(x),y, 1, 1)             # Top row = 1, object in lattice

    epsilon = 10**-3
    num_iterations = 0
    start_time = time.time()

    # Question 2.1
    while(l.delta > epsilon):
        l.SOR(omega = 1)
        l.grow_object(eta = 1)
        num_iterations+=1

    print("--- %s seconds ---" % (time.time() - start_time))
    print(">>>>",len(np.where(l.objects == 1)[0]))
    print("Number of iterations needed to get here = ", num_iterations)
    l.print_lattice()

    # iterations = 1
    # Question 2.2
    # while (iterations < 100000):
    #     l.walker()
    #     iterations += 1

    # l.print_lattice()
if __name__ == "__main__":
    main()
