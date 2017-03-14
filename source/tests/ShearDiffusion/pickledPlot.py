import numpy as np
import sys
import os
sys.path.append('../../../')
from plotting import *
import pickle

pwd = os.getcwd()

i = int(sys.argv[1]); # Get file number

mesh, data, t, velavg = pickle.load( open(pwd+'/output'+'{:04d}'.format(i)+'.p', "rb" ) )

filename = pwd+'/velavg'+'{:04d}'.format(i)+'.dat'
target = open(filename, 'a')
target.write('Time = '+str(t))
for elem in velavg:
	target.write(str(elem))
	target.write("\n")
target.close()

ax = plt.gca()
save_frame(mesh,data.tracer,'Tracer',t,ax,pwd+'/tracer'+'{:04d}'.format(i)+'.png',False)
save_frame(mesh,0.5*(data.u_vel**2+data.v_vel**2),'Energy',t,ax,pwd+'/energy'+'{:04d}'.format(i)+'.png',False)
save_frame(mesh,data.press,'Pressure',t,ax,pwd+'/pressure'+'{:04d}'.format(i)+'.png',False)


