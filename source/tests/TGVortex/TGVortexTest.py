import numpy as np
import sys
sys.path.append('../')
from mesh import *
from plotting import *
from solver import *


#N = 80**2
L_x = 2.*np.pi
L_y = 2.*np.pi
dt = 0.01
Tend = 0.01
nu = 0.1
rho = 1.

N0 = 10
Nend = 1000
Nstep = 10
for NN in range(N0, Nend, Nstep):
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
			t += dt


	u_exact = np.zeros(N);
	for i in range(N):
		u_exact[i] = np.sin(mesh.centroid[i][0])*np.cos(mesh.centroid[i][1])* np.exp(-2.*nu*Tend)

	error = np.sqrt(np.sum(np.dot(np.square(data.u_vel - u_exact),mesh.area))) / np.sqrt(np.sum(np.dot(np.square(u_exact),mesh.area)))
	#print 'Error = '+str(error)

	# Plot solutions and errors
	#ax1 = plt.figure(1)
	#make_frame(mesh,data.u_vel,'u_final', ax1)
	#ax2 = plt.figure(2)
	#make_frame(mesh,u_exact,'u_exact',ax2)
	#ax3 = plt.figure(3)
	#make_frame(mesh,data.u_vel - u_exact,'Error',ax3)


	# Write error to file
	filename = 'TGVortexErrorOneStep.txt'
	target = open(filename, 'a')
	line = np.array([N, error])
	target.write("\t".join(str(elem) for elem in line))
	target.write("\n")
	target.close()

	print line


