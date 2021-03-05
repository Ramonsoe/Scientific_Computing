import numpy as np
import matplotlib.pyplot as plt
import copy

class Grayscott:
    def __init__(self, N, t, dt, dx, du, dv, f, k):
        self.N = N
        self.t = t
        self.dt = dt
        self.dx = dx
        self.du = du
        self.dv = dv
        self.f = f
        self.k = k
        self.L = 1/N
        self.dx = self.dy = 1/N
        size = (self.N, self.N, t)
        self.lattice_u = np.zeros([self.N, self.N, self.t])
        self.lattice_v = np.zeros([self.N, self.N, self.t])
        self.delta = 100

    def initialize_lattice(self, x, y, xlen, ylen):

        self.lattice_u[:,:,0] = 0.5
        self.lattice_v[:,:,0] = 0

        for i in range(x, x+xlen):
            for j in range(y, y+ylen):

                self.lattice_v[i, j, 0] = 0.25


    def do_grayscott(self):
        for i in range(1, self.t-1):
            for j in range(1, self.N-1):
                for k in range(1, self.N-1):
                    if j == 0 or k == 0 or j == self.N-1 or k == self.N-1:
                        next_u, next_v = self.next_step_boundary(i, j, k)
                        self.lattice_u[k, j, i + 1] = next_u
                        self.lattice_v[k, j, i + 1] = next_v
                    else:
                        next_u, next_v = self.next_step(i, j, k)
                        # print(next_u, next_v)
                        self.lattice_u[k, j, i + 1] = next_u
                        print(self.lattice_u[k,j,i+1])
                        self.lattice_v[k, j, i + 1] = next_v

    # General equation
    def next_step(self, i, j, k):

        u = self.lattice_u
        v = self.lattice_v
        next_timestep_u_diffusion = u[j, k, i] + ((self.dt * self.du)/(self.dx**2)) * (u[j+1, k, i] + u[j-1, k, i] + u[j, k+1, i] + u[j, k-1, i] - 4*u[j, k, i])
        next_timestep_u_reaction = - u[j, k, i] * (v[j, k, i])**2 + self.f * (1 - u[j, k, i])
        next_timestep_u = next_timestep_u_diffusion + next_timestep_u_reaction

        next_timestep_v_diffusion = v[j, k, i] + (self.dt * self.dv)/(self.dx**2) * (v[j+1, k, i] + v[j-1, k, i] + v[j, k+1, i] + v[j, k-1, i] - 4*v[j, k, i])
        next_timestep_v_reaction = u[j, k, i] * (v[j, k, i])**2 - (self.f + self.k) * v[j, k, i]
        next_timestep_v = next_timestep_v_diffusion + next_timestep_v_reaction

        # print(next_timestep_u, next_timestep_v)
        return next_timestep_u, next_timestep_v

    # Boundary equation
    def next_step_boundary(self, i, j, k):

        border = self.check_boundary(j, k)
        u = self.lattice_u
        v = self.lattice_v
        if border == 'bottom_row':
            next_timestep_u_diffusion = u[j, k, i] + (self.dt * self.du)/(self.dx**2) * (u[j+1, k, i] + u[j, k+1, i] + u[j, k-1, i] - 3*u[j, k, i])
            next_timestep_u_reaction = - u[j, k, i] * (v[j, k, i])**2 + self.f * (1 - u[j, k, i])
            next_timestep_u = next_timestep_u_diffusion + next_timestep_u_reaction

            next_timestep_v_diffusion = v[j, k, i] + (self.dt * self.dv)/(self.dx**2) * (v[j+1, k, i] + v[j, k+1, i] + v[j, k-1, i] - 3*v[j, k, i])
            next_timestep_v_reaction = u[j, k, i] * (v[j, k, i])**2 - (self.f + self.k) * v[j, k, i]
            next_timestep_v = next_timestep_v_diffusion + next_timestep_v_reaction

        elif border == 'top_row':
            next_timestep_u_diffusion = u[j, k, i] + (self.dt * self.du)/(self.dx**2) * (u[j-1, k, i] + u[j, k+1, i] + u[j, k-1, i] - 3*u[j, k, i])
            next_timestep_u_reaction = - u[j, k, i] * (v[j, k, i])**2 + self.f * (1 - u[j, k, i])
            next_timestep_u = next_timestep_u_diffusion + next_timestep_u_reaction

            next_timestep_v_diffusion = v[j, k, i] + (self.dt * self.dv)/(self.dx**2) * (v[j-1, k, i] + v[j, k+1, i] + v[j, k-1, i] - 3*v[j, k, i])
            next_timestep_v_reaction = u[j, k, i] * (v[j, k, i])**2 - (self.f + self.k) * v[j, k, i]
            next_timestep_v = next_timestep_v_diffusion + next_timestep_v_reaction

        elif border == 'left_column':
            next_timestep_u_diffusion = u[j, k, i] + (self.dt * self.du)/(self.dx**2) * (u[j+1, k, i] + u[j-1, k, i] + u[j, k+1, i] - 3*u[j, k, i])
            next_timestep_u_reaction = - u[j, k, i] * (v[j, k, i])**2 + self.f * (1 - u[j, k, i])
            next_timestep_u = next_timestep_u_diffusion + next_timestep_u_reaction

            next_timestep_v_diffusion = v[j, k, i] + (self.dt * self.dv)/(self.dx**2) * (v[j+1, k, i] + v[j-1, k, i] + v[j, k+1, i] - 3*v[j, k, i])
            next_timestep_v_reaction = u[j, k, i] * (v[j, k, i])**2 - (self.f + self.k) * v[j, k, i]
            next_timestep_v = next_timestep_v_diffusion + next_timestep_v_reaction

        elif border == 'right_column':
            next_timestep_u_diffusion = u[j, k, i] + (self.dt * self.du)/(self.dx**2) * (u[j+1, k, i] + u[j-1, k, i] + u[j, k-1, i] - 3*u[j, k, i])
            next_timestep_u_reaction = - u[j, k, i] * (v[j, k, i])**2 + self.f * (1 - u[j, k, i])
            next_timestep_u = next_timestep_u_diffusion + next_timestep_u_reaction

            next_timestep_v_diffusion = v[j, k, i] + (self.dt * self.dv)/(self.dx**2) * (v[j+1, k, i] + v[j-1, k, i] + v[j, k-1, i] - 3*v[j, k, i])
            next_timestep_v_reaction = u[j, k, i] * (v[j, k, i])**2 - (self.f + self.k) * v[j, k, i]
            next_timestep_v = next_timestep_v_diffusion + next_timestep_v_reaction

        elif border == 'left_bottom':
            next_timestep_u_diffusion = u[j, k, i] + (self.dt * self.du)/(self.dx**2) * (u[j+1, k, i] + u[j, k+1, i] - 2*u[j, k, i])
            next_timestep_u_reaction = - v[j, k, i] * (v[j, k, i])**2 + self.f * (1 - u[j, k, i])
            next_timestep_u = next_timestep_u_diffusion + next_timestep_u_reaction

            next_timestep_v_diffusion = v[j, k, i] + (self.dt * self.dv)/(self.dx**2) * (v[j+1, k, i] + v[j, k+1, i] - 2*v[j, k, i])
            next_timestep_v_reaction = u[j, k, i] * (v[j, k, i])**2 - (self.f + self.k) * v[j, k, i]
            next_timestep_v = next_timestep_v_diffusion + next_timestep_v_reaction

        elif border == 'right_bottom':
            next_timestep_u_diffusion = u[j, k, i] + (self.dt * self.du)/(self.dx**2) * (u[j+1, k, i] + u[j, k-1, i] - 2*u[j, k, i])
            next_timestep_u_reaction = - u[j, k, i] * (v[j, k, i])**2 + self.f * (1 - u[j, k, i])
            next_timestep_u = next_timestep_u_diffusion + next_timestep_u_reaction

            next_timestep_v_diffusion = v[j, k, i] + (self.dt * self.dv)/(self.dx**2) * (v[j+1, k, i] + v[j, k-1, i] - 2*v[j, k, i])
            next_timestep_v_reaction = u[j, k, i] * (v[j, k, i])**2 - (self.f + self.k) * v[j, k, i]
            next_timestep_v = next_timestep_v_diffusion + next_timestep_v_reaction

        elif border == 'left_top':
            next_timestep_u_diffusion = u[j, k, i] + (self.dt * self.du)/(self.dx**2) * (u[j-1, k, i] + u[j, k+1, i] - 2*u[j, k, i])
            next_timestep_u_reaction = - u[j, k, i] * (v[j, k, i])**2 + self.f * (1 - u[j, k, i])
            next_timestep_u = next_timestep_u_diffusion + next_timestep_u_reaction

            next_timestep_v_diffusion = v[j, k, i] + (self.dt * self.dv)/(self.dx**2) * (v[j-1, k, i] + v[j, k+1, i] - 2*v[j, k, i])
            next_timestep_v_reaction = u[j, k, i] * (v[j, k, i])**2 - (self.f + self.k) * v[j, k, i]
            next_timestep_v = next_timestep_v_diffusion + next_timestep_v_reaction

        elif border == 'right_top':
            next_timestep_u_diffusion = u[j, k, i] + (self.dt * self.du)/(self.dx**2) * (u[j-1, k, i] + u[j, k-1, i] - 2*u[j, k, i])
            next_timestep_u_reaction = - u[j, k, i] * (v[j, k, i])**2 + self.f * (1 - u[j, k, i])
            next_timestep_u = next_timestep_u_diffusion + next_timestep_u_reaction

            next_timestep_v_diffusion = v[j, k, i] + (self.dt * self.dv)/(self.dx**2) * (v[j-1, k, i] + v[j, k-1, i] - 2*v[j, k, i])
            next_timestep_v_reaction = v[j, k, i] * (v[j, k, i])**2 - (self.f + self.k) * v[j, k, i]
            next_timestep_v = next_timestep_v_diffusion + next_timestep_v_reaction

        else:
            print('huh hoe kan dit')
        # next_timestep_u_diffusion = self.lattice[j, k, i][0] + (self.dt * self.du)/(self.dx**2) * (self.lattice[j+1, k, i][0] + self.lattice[j-1, k, i][0] + self.lattice[j, k+1, i][0] + self.lattice[j, k-1, i][0] - 4*self.lattice[j, k, i][0])
        # next_timestep_u_reaction = - self.lattice[j, k, i][0] * (self.lattice[j, k, i][1])**2 + self.f * (1 - self.lattice[j, k, i][0])
        # next_timestep_u = next_timestep_u_diffusion + next_timestep_u_reaction
        #
        # next_timestep_v_diffusion = self.lattice[j, k, i][1] + (self.dt * self.dv)/(slef.dx**2) * (self.lattice[j+1, k, i][1] + self.lattice[j-1, k, i][1] + self.lattice[j, k+1, i][1] + self.lattice[j, k-1, i][1] - 4*self.lattice[j, k, i][1])
        # next_timestep_v_reaction = self.lattice[j, k, i][0] * (self.lattice[j, k, i][1])**2 - (self.f + self.k) * self.lattice[j, k, i][1]
        # next_timestep_v = next_timestep_v_diffusion + next_timestep_v_reaction
        return next_timestep_u, next_timestep_v

    def check_boundary(self, j, k):

        if j == 0:
            boundary = 'bottom_row'

        if j == self.N - 1:
            boundary = 'top_row'

        if k == 0:
            boundary = 'left_column'

        if k == self.N -1:
            boundary = 'right_column'

        if j == 0 and k == 0:
            boundary = 'left_bottom'

        if j == 0 and k == self.N - 1:
            boundary = 'right_bottom'

        if j == self.N - 1 and k == 0:
            boundary = 'left_top'

        if j == self.N - 1 and k == self.N -1:
            boundary = 'right_top'

        # else:
        #     print('dit kan niet kloppen')

        return boundary

    # def print_lattice(self):
