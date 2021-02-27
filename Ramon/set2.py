#!/usr/local/bin/python3

import numpy as np
import matplotlib.pyplot as plt
import copy
from lattice import Lattice
def main():
    # Initialize lattice object:
    l = Lattice(100)

    # Initializing ojbect in the lattice (bacteria) and empty lattice:
    x = int(0.99*l.N)
    y = int(0.49*l.N)

    l.initial_lattice(x,y, 1, 1)             # Top row = 1, object in lattice


    epsilon = 10**-5
    num_iterations = 0
    while(l.delta > epsilon):
    # while(num_iterations < 1):
        l.SOR(omega = 1)
        # print(l.lattice)
        l.grow_object(eta = 0.5)
        num_iterations+=1

    # print(l.lattice[8,4])
    # print(np.where(l.lattice == 0.0001283884048461914))
    print("Number of iterations needed to get here = ", num_iterations)
    l.print_lattice()

if __name__ == "__main__":
    main()
