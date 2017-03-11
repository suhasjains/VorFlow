import numpy as np
import sys
sys.path.append('../../')
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
		data.u_vel[i] = 1.0
		data.v_vel[i] = 1.0

t = 0.
plt.ion()
ax = plt.gca()
site0 = mesh.site.copy()
plot_mesh(mesh,ax)
while t < Tend:
		data = time_step(mesh,data,dt,nu)
		plot_mesh(mesh,ax)
		t += dt
		x_exact = np.zeros((N,2))
		x_exact[:,0] = (data.u_vel*t + site0[:,0])%L_x
		x_exact[:,1] = (data.v_vel*t + site0[:,1])%L_y
		mesh.update_mesh(data, dt)
		print 'error:',np.linalg.norm((mesh.site-x_exact)/N)



