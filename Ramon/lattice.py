
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
        self.walker = np.zeros((self.N, self.N))
        self.delta = 100


    def initial_lattice(self, x, y, x_len, y_len):
        self.lattice[0, :] = 1
        self.objects[x:x+x_len, y:y+y_len] = 1

    def SOR(self, omega):
        m = np.copy(self.lattice)
        for i in range(1, self.N-1):
            # Left boundary:
            # m[i, 0] = (omega/4) * (m[i+1, 0] + m[i-1, 0] + m[i, 1] + m[i, self.N-1]) + (1 - omega) * m[i, 0]

            # Periodic Boundaries:
            j = 0
            m[i, j] = m[i, self.N-1] = omega/4 * (m[i+1, j] + m[i-1, j] + m[i, j+1] + m[i, self.N-2]) + (1 - omega) * m[i, j]

            # Middle part:
            for j in range(1, self.N-1):
                if (self.objects[i,j] == 0):
                    m[i, j] = (omega/4) * (m[i+1, j] + m[i-1, j] + m[i, j+1] + m[i, j-1]) + (1 - omega) * m[i, j]

            # Right boundary:
            # m[i, self.N-1] = (omega/4) * (m[i+1, self.N-1] + m[i-1, self.N-1] + m[i, self.N-2] + m[i, 0]) + (1 - omega) * m[i, self.N-1]

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
                elif potential_object_j == self.N or potential_object_j < 0:
                    pass
                elif self.objects[potential_object_i, potential_object_j] == 1:
                    pass
                else:
                    potential_object_concentration = self.lattice[potential_object_i, potential_object_j]
                    # print("dit punt:", potential_object_concentration, potential_object_i, potential_object_j)
                    # print(self.lattice[potential_object_i, potential_object_j])
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

    def random_walk(self):
        """
        Function that performs a random step.
        """
        # check if there is a walking object or there has to be initialized one
        if len(np.where(np.array(self.walker) == 1)[0]) == 0:
            self.initialize_walker()
        else:
            direction = self.direction()
            current_i, current_j = self.current_position()
            new_i, new_j = self.next_point(direction, current_i, current_j)
            if self.check_sinks(new_i) == True:
                self.initialize_walker()
                return False
            elif self.check_periodic(new_j) == True:
                new_j = self.periodic_j(new_j)
            if self.check_object(new_i, new_j) == True:
                self.stick(new_i, new_j)
                print("jaja bingo")
                self.walker[current_i, current_j] = 0
                self.lattice[new_i, new_j] = 0
                return False
            else:
                self.do_step(new_i, new_j)
                self.walker[current_i, current_j] = 0
                return False


    def initialize_walker(self):
        """
        Initialize a random walker object at the top boundary.
        """

        self.walker = np.zeros((self.N, self.N))

        p = [1/self.N for i in range(self.N)]

        new_object_index = np.random.choice([i for i in range(self.N)], 1, p=p)

        self.walker[0, new_object_index] = 1

    def direction(self):
        """
        Determine the direction of the step.
        """

        # lef, right, up down
        directions = ['left', 'right', 'up', 'down']

        p = [1/4 for i in range(len(directions))]
        # p = [0,0,1,0]
        step_direction = np.random.choice(directions, 1, p=p)

        return step_direction

    def current_position(self):
        """
        Det. current position.
        """

        current_position = np.where(self.walker == 1)

        return current_position[0], current_position[1]

    def next_point(self, direction, current_i, current_j):
        """
        Determine next grid point to walk to.
        """

        if direction == 'left':
            new_i = current_i
            new_j = current_j - 1

        elif direction == 'right':
            new_i = current_i
            new_j = current_j + 1

        elif direction == 'up':
            new_i = current_i + 1
            new_j = current_j

        elif direction == 'down':
            new_i = current_i - 1
            new_j = current_j

        return new_i, new_j

    def check_sinks(self, new_i):
        """
        Check if new position if boundary
        """

        if new_i == self.N-1 or new_i == 0 or new_i == -1:
            return True

        else:
            return False

    def check_periodic(self, new_j):
        """
        Check if new point is at the left or right boundary.
        """

        if new_j == -1 or new_j == self.N:
            return True
        else:
            return False

    def periodic_j(self, new_j):
        """
        Change new j with periodic conditions.
        """

        if new_j == -1:
            new_j = self.N - 1
        elif new_j == self.N:
            new_j = 0

        # print(new_j)
        return new_j

    def check_object(self, new_i, new_j):
        """
        Check if walker neighbours object.
        """

        if new_j + 1 == self.N:
            check_j = 1
            if self.objects[new_i+1,new_j] == 1 or self.objects[new_i-1,new_j] == 1 or self.objects[new_i, check_j] == 1 or self.objects[new_i, new_j-1] == 1:
                return True
        elif new_j - 1 == -1:
            check_j = 1
            if self.objects[new_i+1,new_j] == 1 or self.objects[new_i-1,new_j] == 1 or self.objects[new_i, new_j+1] == 1 or self.objects[new_i, check_j] == 1:
                return True

        else:
            if self.objects[new_i+1,new_j] == 1 or self.objects[new_i-1,new_j] == 1 or self.objects[new_i, new_j+1] == 1 or self.objects[new_i, new_j-1] == 1:
                return True

    def stick(self, new_i, new_j):
        """
        Stick walker to object.
        """

        self.objects[new_i, new_j] = 1

    def do_step(self, new_i, new_j):
        """
        Perform step.
        """
        self.walker[new_i, new_j] = 1

    def print_lattice(self):
        fig, ax = plt.subplots()
        self.lattice[self.lattice==0] = np.nan
        ax.imshow(self.lattice, cmap="rainbow", interpolation='nearest', extent=[0,1,0,1], aspect='auto')
        ax.set_aspect(1)

        plt.xlabel('x')
        plt.ylabel('y')
        plt.grid()
        plt.show()

    def print_object(self):
        fig, ax = plt.subplots()
        self.objects[self.objects==0] = np.nan
        ax.imshow(self.objects, cmap="rainbow", interpolation='nearest', extent=[0,1,0,1], aspect='auto')
        ax.set_aspect(1)

        plt.xlabel('x')
        plt.ylabel('y')
        plt.grid()
        plt.show()
