import numpy as np
import sys
sys.path.append('../../')
from mesh import *
from plotting import *
from solver import *


N = 100;
L = 1.;



mesh = Mesh(N,L,L,np.ones(4),'cartesian'); #Use either "random", "cartesian" or "nonuniform" for mesh_type

test = ['FAIL!', 'PASS!'];

if mesh.mesh_type == 'cartesian':
                print '\nCartesian Tests: \n'
                
                # N_neighbours
                print 'N_neighbours: '+test[np.all(mesh.N_neighbor == 4)];
                
                # length
                print 'length: '+test[np.all(np.abs(np.sum(np.asarray(mesh.length),axis=2)) == L / np.sqrt(N))];
                
                # face
                print 'face: '+test[np.all(np.asarray(mesh.face) == L / np.sqrt(N))];
                
                # face_center
                print 'face_center: '+test[np.sqrt(mesh.face_center[0][0,0]**2 + mesh.face_center[0][0,1]**2) == L / (2.*np.sqrt(N)) ];
                
                # area
                print 'area: '+test[np.all(mesh.area == L**2 / N)];
                
                # centroid
                print 'centroid: '+test[np.all(mesh.centroid == mesh.site)];

elif mesh.mesh_type == 'nonuniform':
                print '\nNonuniform Tests: \n'
                
                # N_neighbours
                print 'N_neighbours: '+test[np.all(mesh.N_neighbor == 4)];
               
                # site
                print 'Site: '+test[mesh.site[0,0]==0.125**2];

                # face_center
                print 'face_center: '+test[np.sqrt(mesh.face_center[0][3,0]**2 + mesh.face_center[0][3,1]**2) == (0.375**2 - 0.125**2)/2. ];

plt.ion();
ax = plt.gca();

plot_mesh(mesh,ax);

