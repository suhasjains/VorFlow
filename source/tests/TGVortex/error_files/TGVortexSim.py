import numpy as np
import sys
sys.path.append('../')
from mesh import *
from plotting import *
from solver import *


#N = 80**2
L_x = 2.*np.pi
L_y = 2.*np.pi
dt = 0.005
Tend = 5.
nu = 0.1
rho = 1.

NN = 150
N = NN**2
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
		
	
	u_exact = np.zeros(N);
	for i in range(N):
		u_exact[i] = np.sin(mesh.centroid[i][0])*np.cos(mesh.centroid[i][1])* np.exp(-2.*nu*t)
	error = np.sqrt(np.sum(np.dot(np.square(data.u_vel - u_exact),mesh.area))) / np.sqrt(np.sum(np.dot(np.square(u_exact),mesh.area)))

	# Write error to file
	filename = 'TGVortexErrorSim50N.txt'
	target = open(filename, 'a')
	line = np.array([t, error])
	target.write("\t".join(str(elem) for elem in line))
	target.write("\n")
	target.close()

	# Update time
	print line
	t += dt
		

