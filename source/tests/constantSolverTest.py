import numpy as np
import sys
sys.path.append('../')
from mesh import *
from plotting import *
from solver import *
import matplotlib as mpl


N = 100
L_x = 1.
L_y = 1.
dt = 0.01
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
centroid0 = mesh.centroid
make_frame(mesh,data.u_vel**2 + data.v_vel**2,'Energy',ax)
while t < Tend:
		Dx, Dy, L, Gx, Gy = time_step(mesh)
		data = solve(data, Dx, Dy, L, Gx, Gy, dt, nu)
		mesh.update_mesh(data, dt)
		make_frame(mesh,data.u_vel**2 + data.v_vel**2,'Energy',ax)
		plt.pause(0.005)
		t += dt
		dx = np.zeros((N,2))
		dx[:,0] = data.u_vel*t
		dx[:,1] = data.v_vel*t
		print np.linalg.norm((mesh.centroid-(centroid0+dx))/N)

