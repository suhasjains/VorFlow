import numpy as np
import sys
sys.path.append('../')
from mesh import *
from plotting import *
from solver import *


N = 100;
L_x = 1.;
L_y = 1.;
dt = 0.01;
Tend = 0.5;

mesh = Mesh(N,L_x,L_y,np.zeros(4));

print len(mesh.voronoi.points);

plt.ion()
fig = plt.figure()
plot_mesh(mesh);
plt.pause(15);






