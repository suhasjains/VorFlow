import numpy as np
from scipy.spatial import *
import matplotlib.pyplot as plt

N = 100;
Lx = 2.;
Ly = 2.;
#points = np.array([[0, 0], [0, 1], [0, 2], [1, 0], [1, 1], [1, 2], [2, 0], [2, 1], [2, 2]])
points = np.random.rand(N,2);
points[:,0] = Lx*points[:,0];
points[:,1] = Ly*points[:,1];
# This is a little odd for aspectratio ~= 1...

vor = Voronoi(points)


voronoi_plot_2d(vor)
plt.xlim([0,Lx]);
plt.ylim([0,Ly]);
plt.show()
