import numpy as np
import sys
sys.path.append('../')
from mesh import *
from plotting import *
from solver import *


N = 40**2
L_x = 2.*np.pi
L_y = 2.*np.pi
dt = 0.01
Tend = 0.1
nu = 0.001
rho = 1.

mesh = Mesh(N,L_x,L_y,np.zeros(4),'random')

data = Data(N);

for i in range(N):
		data.u_vel[i] = np.sin(mesh.centroid[i][0])*np.cos(mesh.centroid[i][1])
		data.v_vel[i] = -np.cos(mesh.centroid[i][0])*np.sin(mesh.centroid[i][1])

t = 0.
#plt.ion()
#make_frame(mesh,data.u_vel,'u')
#plt.pause(0.005)
while t < Tend:
		data = time_step(mesh,data,dt,nu)
		#make_frame(mesh,data.u_vel,'u')
		#plt.pause(0.005)
		mesh.update_mesh(data, dt)
		t += dt


u_exact = np.zeros(N);
for i in range(N):
	u_exact[i] = np.sin(mesh.centroid[i][0])*np.cos(mesh.centroid[i][1])

print 'Error = '+str(np.sqrt(np.sum(np.dot(np.square(data.u_vel - u_exact),mesh.area))))
#make_frame(mesh,data.u_vel,'u_final')
#make_frame(mesh,u_exact,'u_exact')
make_frame(mesh,data.u_vel - u_exact,'Error')
