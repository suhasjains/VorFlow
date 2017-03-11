import numpy as np
import sys
sys.path.append('../')
from mesh import *
from plotting import *
from solver import *
 
import matplotlib.pyplot as plt 

n = 21;

t = range(n);
t1 = range(n);
N = range(n);

for i in range(n):
    N[i] = 2**i;
    L = 1.;

    tic = timeit.default_timer()
    mesh = Mesh(N[i],L,L,np.zeros(4),'random'); #Use either "random", "cartesian" or "nonuniform" for mesh_type
    toc = timeit.default_timer()
    #print 'Single meshing time for ' +str(2**i)+' N: {:.2e}'.format((toc-tic))+' s'
    print 'Meshing for N = ' +str(N[i])+' N: {:.2e}'.format((toc-tic))+' s' 
    t[i] = toc-tic;
    t1[i] = t[0]*2**(i-1);
   
plt.hold(True);   
actual, = plt.loglog(N,t,'*-',basex=2,basey=2, label='actual');
ideal, = plt.loglog(N,t1,'--',basex=2,basey=2, label='ideal O(N)');
plt.legend(handles=[actual,ideal])
plt.xlabel('Grid size')
plt.ylabel('Time taken in seconds')
plt.show();

