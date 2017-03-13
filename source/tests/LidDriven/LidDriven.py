import numpy as np
import sys
sys.path.append('../../')
from mesh import *
from plotting import *
from solver import *
import pickle

Nx = 10;

N = Nx**2
L_x = 1.
L_y = 1.
dt = 1./Nx
dTPlot = 5*dt
Tend = 10.
nu = 1e-6
rho = 1.

BCu = [1,1,1,1]
BCuvals = [0.,0.,0.,1.]
BCv = [1,1,1,1]
BCvvals = [0.,0.,0.,0.]

mesh = Mesh(N,L_x,L_y,BCu+BCv,'random')
data = Data(N);
data.tracer = np.zeros(N);

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
			pickle.dump([mesh,data,t],fileObject)
			fileObject.close()
			i += 1
		
		data = time_step(mesh,data,dt,nu,BCu,BCuvals,BCv,BCvvals)	
		mesh.update_mesh(data, dt)
		t += dt
