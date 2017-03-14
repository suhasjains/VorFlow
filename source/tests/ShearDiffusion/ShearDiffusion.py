import numpy as np
import sys
sys.path.append('../../')
from mesh import *
from plotting import *
from solver import *
import pickle

Nx = 150;

N = Nx**2
L_x = 1.
L_y = 1.
dt = 1./Nx;
dTPlot = 0.05
Tend = 10.
nu = 1 #Re=1
rho = 1.

mesh = Mesh(N,L_x,L_y,np.zeros(4),'random')
data = Data(N);
data.tracer = np.zeros(N);

Nbin = 20
velavg = np.zeros((N,1))
ybin = np.linspace(0,1,Nbin+1)

A = 1.;
uPrime = 0.; #No perturbation
k_x = 2.;
width = 1./3.;
centre = 0.5;
smoothing = 200.;

for i in range(N):
		y = mesh.centroid[i][1] - centre;
		x = mesh.centroid[i][0];
		data.u_vel[i] = 0.5*A*(np.tanh(smoothing*(y+0.5*width)) - np.tanh(smoothing*(y-0.5*width))) - 0.5*A; 
		data.tracer[i] = data.u_vel[i]/A + 0.5


t = 0.
#data = time_step(mesh,data,0.,nu) # Project Pressure...
tprint = 0.
#ax = plt.gca()
i = 0;
while t < Tend:
		if t >= tprint:
			tprint += dTPlot;
			#save_frame(mesh,data.tracer,'Tracer',t,ax,'output/tracer'+'{:04d}'.format(i)+'.png',False)
			#save_frame(mesh,0.5*(data.u_vel**2+data.v_vel**2),'Energy',t,ax,'output/energy'+'{:04d}'.format(i)+'.png',False)
			#save_frame(mesh,data.press,'Pressure',t,ax,'output/pressure'+'{:04d}'.format(i)+'.png',False)
			filename = 'output/output'+'{:04d}'.format(i)+'.p'
			fileObject = open(filename,'wb')
			pickle.dump([mesh,data,t,velavg],fileObject)
			fileObject.close()
			i += 1
		
		data = time_step(mesh,data,dt,nu)	
		mesh.update_mesh(data, dt)
		t += dt

		# bin data
		velavg = np.zeros((N,1))
		velcount = np.zeros((N,1))
		for i in range(N):
			ypt = mesh.centroid[i][1]
			for j in range(Nbin):
				if ybin[j]<=ypt and ypt<ybin[j+1]:
					velavg += data.u_vel[i]
					velcount += 1
		velavg = np.divide(velavg,velcount)	
