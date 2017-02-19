import numpy as np
from scipy.spatial import *
import matplotlib.pyplot as plt

N = 100;
Lx = 1.;
Ly = 1.;
#points = np.array([[0, 0], [0, 1], [0, 2], [1, 0], [1, 1], [1, 2], [2, 0], [2, 1], [2, 2]])
points = np.random.rand(N,2);
pp = np.zeros([9*N,2]);
pp[:,0] = Lx*np.concatenate((points[:,0]-1,points[:,0],points[:,0]+1,points[:,0]-1,points[:,0],points[:,0]+1,points[:,0]-1,points[:,0],points[:,0]+1))
pp[:,1] = Ly*np.concatenate((points[:,1]-1,points[:,1]-1,points[:,1]-1,points[:,1],points[:,1],points[:,1],points[:,1]+1,points[:,1]+1,points[:,1]+1))
#points[:,1] = Ly*points[:,1];

# This is a little odd for aspectratio ~= 1...

vor = Voronoi(pp,qhull_options='Qbb Qc')


voronoi_plot_2d(vor)
plt.axis('scaled')
plt.xlim([0,Lx]);
plt.ylim([0,Ly]);
plt.show()
