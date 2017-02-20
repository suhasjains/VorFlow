from definitions import *
from mesh import *
from solver import *
from plotter import *


# Generate Mesh

mesh = Mesh(N,L_x,L_y,BCs)

# Define Initial Conditions



# Solve

Dx, Dy, L, Gx, Gy = time_step(mesh)

data = solve(data, Dx, Dy, L, Gx, Gy, Re, dt)
#data = solve(N,dt,T_end,nu,mesh,BCs,BC_vals)

