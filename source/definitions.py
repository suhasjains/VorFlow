import numpy as np


# Parameter definitions

N = 100;
L_x = 1.;
L_y = 1.;
T_end = 1.;
dt = 0.01;
nu = 0.01; # Nondimensionalise with Re
BCs = np.array([0,0,0,0]);
BC_vals = np.array([-1,-1,-1,-1]);




