import numpy as np
import sys
sys.path.append('../')
from mesh import *
#from plotting import *
from solver import *


#N = 100
N = 9
L_x = 3.
L_y = 3.
dt = 0.1
Tend = 1.
Re = 1000.

#mesh = Mesh(N,L_x,L_y,np.zeros(4),'cartesian')
mesh = Mesh(N,L_x,L_y,np.ones(4),'cartesian')

data = Data(N);

for i in range(N):
		data.u_vel[i] = 1.
		data.v_vel[i] = 1.
		#data.press[i] = np.exp(-((mesh.site[i,0]-L_x/2)**2 + (mesh.site[i,1]-L_y/2)**2))
		data.press[i] = 0.5
#D = np.zeros((N,N))
#Dx = np.zeros((N,N))
#Dy = np.zeros((N,N))
#L = np.zeros((N,N))
#Gx = np.zeros((N,N))
#Gy = np.zeros((N,N))
#G = np.zeros((N,N))

Dx, Dy, L, Gx, Gy = time_step(mesh)

#print Dx
#print Dy
print L
#print Gx
#print Gy

data = solve(data, Dx, Dy, L, Gx, Gy, Re, dt)

print(data.u_vel)
print(data.v_vel)
print(data.press)
