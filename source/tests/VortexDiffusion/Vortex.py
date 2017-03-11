import numpy as np
import sys
sys.path.append('../')
from mesh import *
from plotting import *
from solver import *


Nx = 25;

N = Nx**2
L_x = 1.
L_y = 1.
dt = 1./Nx;
dTPlot = 0.1
Tend = 1.
nu = 1e-4
rho = 1.

mesh = Mesh(N,L_x,L_y,np.zeros(4),'random')
data = Data(N);

Gamma = 1.;
centre = 0.5;
R = 1./4.;
smoothing = 10.;

for i in range(N):
		y = mesh.centroid[i][1] - centre;
		x = mesh.centroid[i][0] - centre;
                r = np.sqrt(x**2 + y**2);
		data.u_vel[i] = -y * (1.+np.tanh(smoothing*(R-r)));
		data.v_vel[i] = x * (1.+np.tanh(smoothing*(R-r)));


t = 0.
tprint = 0.
plt.ion()
ax = plt.gca()
while t < Tend:
		data = time_step(mesh,data,dt,nu)
	        make_frame(mesh,data.u_vel**2 + data.v_vel**2,'Energy',ax,False)
		mesh.update_mesh(data, dt)
		t += dt
