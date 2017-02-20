import numpy as np
import sys
sys.path.append('../')
from mesh import *
from plotting import *
from solver import *


N = 100;
L_x = 1.;
L_y = 1.;
dt = 0.1;
Tend = 1.;

mesh = Mesh(N,L_x,L_y,np.zeros(4));

data = Data(N);

for i in range(N):
		data.u_vel[i] = 1.;
		data.v_vel[i] = 1.;
		data.press[i] = mesh.area[i]; #np.exp(-((mesh.site[i,0]-L_x/2)**2 + (mesh.site[i,1]-L_y/2)**2));


t = 0.;
plt.ion()
fig = plt.figure()
while t < Tend:
		mesh.update_mesh(data,dt);
		#plot_mesh(mesh);
		make_frame(mesh,data.press,'Area')
		plt.pause(0.01);
		t += dt;







