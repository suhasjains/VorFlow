import numpy as np
import sys
sys.path.append('../')
from mesh import *
from plotting import *
from solver import *


N = 75**2
L_x = 2.*np.pi
L_y = 2.*np.pi
dt = 0.1
Tend = 1.
nu = 1.
rho = 1.

mesh = Mesh(N,L_x,L_y,np.zeros(4),'random')

data = Data(N);

for i in range(N):
		data.u_vel[i] = 0.5 + 0.5*np.sin(mesh.centroid[i][1])
		data.v_vel[i] = 1.

t = 0.
while t < Tend:
		data = time_step(mesh,data,dt,nu)
		mesh.update_mesh(data, dt)
		t += dt
