import numpy as np
import sys
sys.path.append('../')
from mesh import *
from plotting import *
from solver import *

N = 49
L_x = 1.
L_y = 1.
dt = 0.02
Tend = 0.5

mesh = Mesh(N,L_x,L_y,np.zeros(4),'cartesian')

data = Data(N)

t = 0.
plt.ion()
ax = plt.gca()
#fig = plt.figure()
while t < Tend:
                for i in range(N):
                        #data.u_vel[i] = -np.sin(2.*np.pi/L_y * (mesh.centroid[i][1] - L_y/2.))
                        #data.v_vel[i] = np.sin(2.*np.pi/L_x * (mesh.centroid[i][0] - L_x/2.))
		        data.u_vel[i] = -(mesh.centroid[i][1] - L_y/2.)/L_y
                        data.v_vel[i] = (mesh.centroid[i][0] - L_x/2.)/L_x
                        data.press[i] = mesh.area[i]/(2.*(L_x*L_y)/N)

                mesh.update_mesh(data,dt)
		#plot_mesh(mesh);
		make_frame(mesh,data.press,'Area',ax)
		plt.pause(0.005)
		t += dt







