import numpy as np
import sys
sys.path.append('../')
from mesh import *
from plotting import *
from solver import *
import matplotlib as mpl


N = 100
L_x = 2.*np.pi
L_y = 2.*np.pi
dt = 0.1
Tend = 1.
nu = 0.01
rho = 1.

mesh = Mesh(N,L_x,L_y,np.zeros(4),'random')

data = Data(N);

for i in range(N):
		data.u_vel[i] = 0.5
		data.v_vel[i] = 0.5
		data.press[i] = 0.5

t = 0.
plt.ion()
ax = plt.gca()
centroid0 = mesh.centroid.copy()
make_frame(mesh,data.u_vel**2 + data.v_vel**2,'Energy',ax,True)
while t < Tend:
		data = time_step(mesh,data,dt,nu)
		make_frame(mesh,data.u_vel**2 + data.v_vel**2,'Energy',ax,True)
		plt.pause(0.005)
		t += dt
		dx = np.zeros((N,2))
		dx[:,0] = data.u_vel*t
		dx[:,1] = data.v_vel*t
		mesh.update_mesh(data, dt)
		print 'error:',np.linalg.norm((mesh.centroid-(centroid0+dx))/N)


