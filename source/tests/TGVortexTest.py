import numpy as np
import sys
sys.path.append('../')
from mesh import *
from plotting import *
from solver import *


N = 20**2
L_x = 2.*np.pi
L_y = 2.*np.pi
dt = 0.1
Tend = 1.
nu = 0.01
rho = 1.

mesh = Mesh(N,L_x,L_y,np.zeros(4),'random')

data = Data(N);

for i in range(N):
		data.u_vel[i] = np.sin(mesh.centroid[i][0])*np.cos(mesh.centroid[i][1])
		data.v_vel[i] = -np.cos(mesh.centroid[i][0])*np.sin(mesh.centroid[i][1])

t = 0.
plt.ion()
make_frame(mesh,data.u_vel,'u')
plt.pause(0.005)
while t < Tend:
		data = time_step(mesh,data,dt,nu)
		make_frame(mesh,data.u_vel,'u')
		mesh.update_mesh(data, dt)
		plt.pause(0.005)
		t += dt


