
"""
Created on Wed Feb 24 11:44:14 2021

@author: ramonsoe
"""
import numpy as np
import matplotlib.pyplot as plt
import copy

class Lattice:
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

    def grow_object(self, eta):

        np.seterr(all='ignore')
        object_places = np.where(np.array(self.objects) == 1)

        concentrations_potential_object = []
        look_around = [[0, 1], [0, -1],[1, 0],[-1,0]]

        for i in range(len(object_places[0])):
            for j in range(4):
                current_i = object_places[0][i]
                current_j = object_places[1][i]

                potential_object_i = current_i - look_around[j][0]
                potential_object_j = current_j - look_around[j][1]

                if potential_object_i == 0 or potential_object_i ==self.N or potential_object_i < 0:
                    pass
                elif potential_object_j == 0 or potential_object_j == self.N or potential_object_j < 0:
                    pass
                elif self.objects[potential_object_i, potential_object_j] == 1:
                    pass
                else:
                    potential_object_concentration = self.lattice[potential_object_i, potential_object_j]
                    # print(potential_object_concentration)
                    concentrations_potential_object.append([potential_object_i, potential_object_j, potential_object_concentration, 0])

        # print(concentrations_potential_object)
        concentrations_potential_object = np.unique(concentrations_potential_object, axis=0)
        sum_total_concentrations = np.sum(concentrations_potential_object**eta, axis=0)[2]

        if sum_total_concentrations == 0:
            return False
        else:
            for i in range(len(concentrations_potential_object)):
                concentration = np.maximum(concentrations_potential_object[i][2]**eta,0)
                growth_probability = concentration / sum_total_concentrations
                concentrations_potential_object[i][3] = growth_probability

        # print(concentrations_potential_object)
        p = [object[3] for object in concentrations_potential_object]
        p /= np.sum(np.array(p))
        # print(p)
        new_object_index = np.random.choice([i for i in range(len(concentrations_potential_object))], 1, p=p)
        new_object = concentrations_potential_object[new_object_index]
        self.objects[int(new_object[0][0]), int(new_object[0][1])] = 1
        self.lattice[int(new_object[0][0]), int(new_object[0][1])] = 0

    def print_lattice(self):
        fig, ax = plt.subplots()
        self.lattice[self.lattice==0] = np.nan
        ax.imshow(self.lattice, cmap="rainbow", interpolation='nearest', extent=[0,1,0,1], aspect='auto')
        ax.set_aspect(1)

        plt.xlabel('x')
        plt.ylabel('y')
        plt.show()
