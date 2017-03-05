import numpy as np
import sys
sys.path.append('../')
from mesh import *
from plotting import *
from solver import *


N = 16384;
L = 1.;

tic = timeit.default_timer()
mesh = Mesh(N,L,L,np.zeros(4),'random'); #Use either "random", "cartesian" or "nonuniform" for mesh_type
toc = timeit.default_timer()
print 'Single meshing time: '+'{:.2e}'.format((toc-tic))+' s'
