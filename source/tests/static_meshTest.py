import numpy as np
import sys
sys.path.append('../')
from mesh import *
from plotting import *
from solver import *


N = 16;
L = 1.;

mesh = Mesh(N,L,L,np.zeros(4),False); # Make Cartesian grid (random=False)

test = ['FAIL!', 'PASS!'];

print '\nCartesian Tests: \n'

# N_neighbours
print 'N_neighbours: '+test[np.all(mesh.N_neighbor == 4)];

# length
print 'length: '+test[np.all(np.abs(np.sum(mesh.length[1],axis=1)) == L / np.sqrt(N))];

# face
print 'face: '+test[np.all(mesh.face[0] == L / np.sqrt(N))];

# face_center
print 'face_center: '+test[np.sqrt(mesh.face_center[0][0,0]**2 + mesh.face_center[0][0,1]**2) == L / (2.*np.sqrt(N)) ];

# area
print 'area: '+test[np.all(mesh.area == L**2 / N)];

# centroid
print 'centroid: '+test[np.all(mesh.centroid == mesh.site)];

plot_mesh(mesh);







