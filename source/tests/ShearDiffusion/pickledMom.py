import numpy as np
import sys
import os
sys.path.append('../../../')
from plotting import *
import pickle

pwd = os.getcwd()

Nt = 50
momtotu = np.zeros(Nt)
momtotv = np.zeros(Nt)
enetot = np.zeros(Nt)
ttot = np.zeros(Nt)

for i in range(Nt):

	print i
	mesh, data, ttot[i], _ = pickle.load( open(pwd+'/output'+'{:04d}'.format(i)+'.p', "rb" ) )
	N = mesh.N
	for j in range(N):
		momtotu[i] += mesh.area[j]*data.u_vel[j]
		momtotv[i] += mesh.area[j]*data.v_vel[j]
		enetot[i] += mesh.area[j]*data.u_vel[j]*data.u_vel[j] + mesh.area[j]*data.v_vel[j]*data.v_vel[j]
		
print momtotu
print momtotv
print enetot

plt.plot(ttot,momtotu-momtotu[0],'-o',ms=10)
plt.plot(ttot,momtotv-momtotv[0],'-s',ms=10)
plt.xlabel('t')
plt.ylabel('Change in total momentum')
plt.ticklabel_format(style='sci',axis='y')
ax = plt.gca()
ax.yaxis.major.formatter.set_powerlimits((0,0))
plt.legend(['x-momentum','y-momentum'],loc=2)
plt.title('Diffusion of shear layer, Re = 10')
plt.savefig('plot3'+'{:04d}'.format(Nt)+'.png')

plt.figure()
plt.plot(ttot,enetot-enetot[0],'-o',ms=10)
plt.xlabel('t')
plt.ylabel('Change in total energy')
plt.ticklabel_format(style='sci',axis='y')
plt.title('Diffusion of shear layer, Re = 10')
plt.savefig('plot4'+'{:04d}'.format(Nt)+'.png')
