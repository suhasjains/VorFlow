import numpy as np
import sys
sys.path.append('../../')
from mesh import *
from plotting import *
from solver import *


Nx = 25;

N = Nx**2
L_x = 1.
L_y = 1.
dt = 0.005;
dTPlot = 0.005
Tend = 10.
nu = 1e-4
rho = 1.

mesh = Mesh(N,L_x,L_y,np.zeros(4),'cartesian')
data = Data(N);

centre = 0.5;
R = 1./3.;
smoothing = 20.;

for i in range(N):
		y = mesh.centroid[i][1] - centre;
		x = mesh.centroid[i][0] - centre;
		r = np.sqrt(x**2 + y**2);
		
		if r < R/2.: # Gresho
			uPhi = r/(R/2.);
		elif r < R:
			uPhi = 2. - r/(R/2.);
		else:
			uPhi = 0.;
		
		#data.u_vel[i] = -y * (1.+np.tanh(smoothing*(R-r)));
		#data.v_vel[i] = x * (1.+np.tanh(smoothing*(R-r)));
		data.u_vel[i] = uPhi * -y/r;
		data.v_vel[i] = uPhi * x/r;

t = 0.
tprint = 0.
#plt.ion()
ax = plt.gca()
i = 0;
while t < Tend:
		data = time_step(mesh,data,dt,nu)
		#make_frame(mesh,data.u_vel**2 + data.v_vel**2,'Energy',ax,False)
		save_frame(mesh,data.u_vel**2 + data.v_vel**2,'Energy',t,ax,'energy'+'{:04d}'.format(i)+'.png',True)
		save_frame(mesh,data.press,'Pressure',t,ax,'pressure'+'{:04d}'.format(i)+'.png',True)
		mesh.update_mesh(data, dt)
		t += dt
		i += 1;
