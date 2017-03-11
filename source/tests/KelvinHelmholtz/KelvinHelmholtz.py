import numpy as np
import sys
sys.path.append('../../')
from mesh import *
from plotting import *
from solver import *

Nx = 120;

N = Nx**2
L_x = 1.
L_y = 1.
dt = 1./Nx;
dTPlot = 0.05
Tend = 10.
nu = 1e-4
rho = 1.

mesh = Mesh(N,L_x,L_y,np.zeros(4),'random')
data = Data(N);
tracer = np.zeros(N);

A = 1.;
uPrime = A/10.;
k_x = 2.;
width = 1./3.;
centre = 0.5;
smoothing = 100.;

for i in range(N):
		y = mesh.centroid[i][1] - centre;
		x = mesh.centroid[i][0];
		data.u_vel[i] = 0.5*A*(np.tanh(smoothing*(y+0.5*width)) - np.tanh(smoothing*(y-0.5*width))) - 0.5*A; 
		data.v_vel[i] = uPrime * np.sin(2.*np.pi*k_x * x/L_x)
		tracer[i] = data.u_vel[i]/A + 0.5


t = 0.
tprint = dTPlot
#plt.ion()
ax = plt.gca()
#make_frame(mesh,data.u_vel,'u',ax,False)
i = 0;
while t < Tend:
		data = time_step(mesh,data,dt,nu)
	        #make_frame(mesh,data.u_vel,'u',ax,False)
		if t > tprint:
			tprint += dTPlot;
			save_frame(mesh,tracer,'Tracer',t,ax,'output/tracer'+'{:04d}'.format(i)+'.png',False)
			save_frame(mesh,0.5*(data.u_vel**2+data.v_vel**2),'Energy',t,ax,'output/energy'+'{:04d}'.format(i)+'.png',False)
			save_frame(mesh,data.press,'Pressure',t,ax,'output/pressure'+'{:04d}'.format(i)+'.png',False)
			i += 1
		
		mesh.update_mesh(data, dt)
		t += dt
