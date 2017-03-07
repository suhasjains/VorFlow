import numpy as np
import sys
sys.path.append('../')
from mesh import *
from plotting import *
from solver import *


N = 160**2
L_x = 1.
L_y = 1.
dt = 0.01
dTPlot = 0.1
Tend = 1.
nu = 1e-4
rho = 1.

mesh = Mesh(N,L_x,L_y,np.zeros(4),'random')
data = Data(N);

width = 1./3.;

for i in range(N):
		y = mesh.centroid[i][1];
		if np.abs(y-0.5*L_y) < 0.5*width:
				data.u_vel[i] = 1.
		else:
				data.u_vel[i] = -1.
		
		data.v_vel[i] = 0.


t = 0.
tprint = 0.
plt.ion()
ax = plt.gca()
make_frame(mesh,data.u_vel,'u',ax,False)
while t < Tend:
		data = time_step(mesh,data,dt,nu)
		if t > tprint:
				make_frame(mesh,data.u_vel,'u',ax,False)
				tprint += dTPlot;
		
		mesh.update_mesh(data, dt)
		t += dt
