import numpy as np
import sys
import os
#sys.path.append($VORFLOW_DIR)
#from mesh import *
from plotting import *
import pickle


pwd = os.getcwd()

N_file = int(sys.argv[1]); # Get file number


v_rms = np.zeros(N_file);
v_exact = np.zeros(N_file);
T = np.zeros(N_file);

k = 2.*2.*np.pi;
sigma = k/2.;

for i in range(N_file):
	mesh, data, t = pickle.load( open(pwd+'/output'+'{:04d}'.format(i)+'.p', "rb" ) )
	T[i] = t;
	v_rms[i] = np.sqrt(np.mean(np.square(data.v_vel[np.abs(mesh.centroid[:,1].flatten() - 1./3.) < 0.1])));
	#v_exact[i] = v_rms[3] * np.exp(t*sigma)

plt.figure()
plt.semilogy(T,v_rms)
plt.semilogy(T,v_rms[3]*np.exp(T*sigma),'r--')
plt.semilogy(T,v_rms[3]*np.exp(T),'r--')
plt.ylim(ymax=0.3)
plt.xlabel(r'$t$')
plt.ylabel(r'$v_\mathrm{rms}$')
plt.savefig('growth.eps')



