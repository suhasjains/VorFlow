from definitions import *
from mesh import *
from solver import *
from plotter import *


# Generate Mesh

mesh = Mesh(N,L_x,L_y,BCs);

# Define Initial Conditions



# Solve

solve(N,dt,T_end,nu,mesh,BCs,BC_vals);

