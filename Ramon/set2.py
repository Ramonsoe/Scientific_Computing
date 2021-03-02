#!/usr/local/bin/python3

import numpy as np
import matplotlib.pyplot as plt
import copy
from lattice import Lattice
def main():
    # Initialize lattice object:
    l = Lattice(50)

    # Initializing ojbect in the lattice (bacteria) and empty lattice:
    x = int(0.99*l.N)
    x = l.N  - 1
    y = int(0.49*l.N)

    l.initial_lattice(int(x),y, 1, 10)             # Top row = 1, object in lattice

    epsilon = 10**-5
    num_iterations = 0
    # while(l.delta > epsilon):
    while(num_iterations < 1000000):
        # l.SOR(omega = 1)
        # l.grow_object(eta = 0.2)
        l.random_walk()
        # print(l.walker)
        # if num_iterations%10000 == 0:
            # print(num_iterations)
        num_iterations+=1

    print(">>>>",len(np.where(l.objects == 1)[0]))
    # print(l.walker)
    print("Number of iterations needed to get here = ", num_iterations)
    # l.print_lattice()
    l.print_object()

if __name__ == "__main__":
    main()
