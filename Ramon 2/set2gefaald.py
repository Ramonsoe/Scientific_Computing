import numpy as np
import matplotlib.pyplot as plt
import copy
import time

class World(object):
    def __init__(self, N, number_of_objects, epsilon):
        self.N = N
        self.dx = self.dy = 1/N
        self.y = []
        self.y.append(0)

        for i in range(1, N - 1):
            self.y.append(self.y[i-1] + self.dy)

        self.y.append(1)

        self.top_row = 1
        self.bottom_row = 0
        self.L = 1/N
        self.epsilon = epsilon
        self.size_solution_matrix = (N, N)
        self.objects = number_of_objects

        self.solution_matrix_current = np.zeros(self.size_solution_matrix)
        self.solution_matrix_next  = np.zeros(self.size_solution_matrix)
        self.object_matrix = np.zeros(self.size_solution_matrix)


        self.solution_matrix_current[N-1, :] = self.top_row
        self.solution_matrix_next[N-1, :] = self.top_row

        temp = np.ones(self.objects)
        middle = round(N/2)
        start = middle - round(0.5*len(temp))

        self.object_matrix[0, start:(start+len(temp))] = temp

def SOR_solver(world, omega):
    """
    Successive Over Relaxation function
    """

    N = world.N

    for i in range(1, N-1):
        # Left boundary:
        world.solution_matrix_next[i, 0] = (omega/4) * (world.solution_matrix_current[i+1, 0] + world.solution_matrix_next[i-1, 0] + world.solution_matrix_current[i, 1] + world.solution_matrix_current[i, N-1]) + (1 - omega) * world.solution_matrix_current[i, 0]

        # Middle part:
        for j in range(1, N-1):
            if (world.object_matrix[i,j] == 0):
                world.solution_matrix_next[i, j] = (omega/4)*(world.solution_matrix_current[i+1, j] + world.solution_matrix_next[i-1, j] + world.solution_matrix_current[i, j+1] + world.solution_matrix_next[i, j-1]) + (1 - omega) * world.solution_matrix_current[i, j]

        # Right boundary:
        world.solution_matrix_next[i, N-1] = (omega/4) * (world.solution_matrix_current[i+1, N-1] + world.solution_matrix_next[i-1, N-1] + world.solution_matrix_next[i, N-2] + world.solution_matrix_next[i, 0]) + (1 - omega) * world.solution_matrix_current[i, N-1]

def solve(world, omega):

    delta = 100
    convergence = True
    count = 0
    while(convergence):
        count +=1
        SOR_solver(world, omega)
        delta = np.max(np.subtract(world.solution_matrix_next, world.solution_matrix_current))

        if (delta < world.epsilon):
            convergence = False
        world.solution_matrix_current = copy.deepcopy(world.solution_matrix_next)

def plot_lattice(world):
    fig, ax = plt.subplots()

    ax.imshow(world.solution_matrix_current, cmap='rainbow', interpolation='nearest', extent=[0,1,0,1], aspect='auto')
    ax.set_aspect(1) # you may also use am.imshow(..., aspect="auto") to restore the aspect ratio
    plt.xlabel('x')
    plt.ylabel('y')
    plt.show()


if __name__ == '__main__':
    epsilon = 10**-5
    world = World(100, 10, epsilon)

    omega = 1.8
    solve(world, omega)

    plot_lattice(world)
