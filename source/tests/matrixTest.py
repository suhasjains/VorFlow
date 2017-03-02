import numpy as np
import sys
sys.path.append('../')
import timeit
from mesh import *
from solver import *

N1 = 300
N = N1**2 
L_x = 3.
L_y = 3.
dt = 0.1
Tend = 1.
Re = 1000.

tic=timeit.default_timer()
mesh = Mesh(N,L_x,L_y,np.zeros(4),'cartesian')
#mesh = Mesh(N,L_x,L_y,np.ones(4),'cartesian')
toc=timeit.default_timer()
print 'meshing done in',toc-tic,'seconds'

data = Data(N);
kx = 1.
ky = 1.

tic=timeit.default_timer()
for i in range(N):
		#data.u_vel[i] = 1. #Dx test1
		#data.u_vel[i] = mesh.site[i,0] #Dx test2
		#data.v_vel[i] = 1. #Dy test1
		#data.v_vel[i] = mesh.site[i,1]*3 #Dy test2
		#data.press[i] = 1. #L test1
		#data.press[i] = np.cos(kx*2.0*np.pi*mesh.site[i,0]/L_x) #L test2
		data.press[i] = np.cos(kx*2.0*np.pi*mesh.site[i,0]/L_x)*np.cos(ky*2.0*np.pi*mesh.site[i,1]/L_y) #L test3
		#data.press[i] = 1. #G test1
		#data.press[i] = mesh.site[i,0] #G test2
#D = np.zeros((N,N))
#Dx = np.zeros((N,N))
#Dy = np.zeros((N,N))
#L = np.zeros((N,N))
#Gx = np.zeros((N,N))
#Gy = np.zeros((N,N))
#G = np.zeros((N,N))
Dx, Dy, L, Gx, Gy = time_step(mesh)
toc=timeit.default_timer()
print 'matrices built in',toc-tic,'seconds'

#print mesh.face
#print mesh.length
#print mesh.area
#print mesh.site[:,0]
#print mesh.site[:,1]
#print data.u_vel
#print data.v_vel
#print np.dot(Dx,data.u_vel)
#print np.dot(Dy,data.v_vel)
#print L
#print data.press
#print np.dot(L,data.press)
#print np.dot(L,data.press)+4.0*np.pi*np.pi/L_x/L_x*data.press
#print np.linalg.norm(np.dot(L,data.press)+4.0*np.pi*np.pi*kx*kx/L_x/L_x*data.press)
print np.linalg.norm(np.dot(L,data.press)+(4.0*np.pi*np.pi*kx*kx/L_x/L_x+4.0*np.pi*np.pi*ky*ky/L_y/L_y)*data.press)
#print np.dot(Gx,data.press)
