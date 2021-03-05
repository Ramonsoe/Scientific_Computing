import numpy as np
import matplotlib.pyplot as plt
import copy
from grayscott import Grayscott
def main():
    # Initialize lattice object:
    t = 3
    dt = 1
    dx = 1
    du = 0.16
    dv = 0.08
    f = 0.035
    k = 0.060
    g = Grayscott(5, t, dt, dx, du, dv, f, k)

    x = int(0.49*g.N)
    # x = g.N
    y = int(0.49*g.N)

    g.initialize_lattice(x, y, 1,1)
    #
    # for i in range(1, t):
    g.do_grayscott()

    # print(g.lattice_u[:,:,2])
    # print(g.lattice_v[:,:,2])
    # for i in range(t):
    #     print(g.lattice[1,1,i][0])

    # fig, ax = plt.subplots()
    #
    # ax.imshow(g.lattice[:, :, 9][0], cmap='hot', interpolation='nearest', extent=[0,1,0,1], aspect='auto')
    # ax.set_aspect(1) # you may also use am.imshow(..., aspect="auto") to restore the aspect ratio
    # plt.xlabel('x')
    # plt.ylabel('y')
    # plt.show()
if __name__ == "__main__":
    main()
