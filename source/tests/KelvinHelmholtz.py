import numpy as np
import sys
sys.path.append('../')
from mesh import *
from plotting import *
from solver import *


N = 200**2
L_x = 1.
L_y = 1.
dt = 0.01
dTPlot = 0.1
Tend = 1.
nu = 1e-4
rho = 1.

mesh = Mesh(N,L_x,L_y,np.zeros(4),'random')
data = Data(N);

A = 1.;
k_x = 2.;
width = 1./3.;
centre = 0.5;
smoothing = 10.;

for i in range(N):
		y = mesh.centroid[i][1] - centre;
		x = mesh.centroid[i][0];
		data.u_vel[i] = 0.5*A*(np.tanh(smoothing*(y+0.5*width)) - np.tanh(smoothing*(y-0.5*width))) - 0.5*A; 
		
		data.v_vel[i] = A/10. * np.sin(2.*np.pi*k_x * x/L_x)


t = 0.
tprint = 0.
#plt.ion()
#ax = plt.gca()
#make_frame(mesh,data.u_vel,'u',ax,False)
while t < Tend:
		data = time_step(mesh,data,dt,nu)
		if t > tprint:
				#make_frame(mesh,data.u_vel,'u',ax,False)
				tprint += dTPlot;
		
		mesh.update_mesh(data, dt)
		t += dt
