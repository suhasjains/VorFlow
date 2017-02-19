import numpy as np
import sys
sys.path.append('../')
from mesh import *
from plotting import *
from solver import *


N = 20;
L_x = 1.;
L_y = 1.;
dt = 0.01;
Tend = 0.5;

mesh = Mesh(N,L_x,L_y,np.zeros(4));

data = Data(N);

t = 0.;
plt.ion()
fig = plt.figure()
plot_mesh(mesh);
while t < Tend:
	mesh.update_mesh(data,dt);
	plot_mesh(mesh);
	#make_frame(mesh,data.press,'Pressure')
	plt.pause(0.01);
	t += dt;







