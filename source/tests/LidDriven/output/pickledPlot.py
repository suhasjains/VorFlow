import numpy as np
import sys
import os
#sys.path.append($VORFLOW_DIR)
#from mesh import *
sys.path.append('../../../')
from plotting import *
import pickle

pwd = os.getcwd()

i = int(sys.argv[1]); # Get file number

mesh, data, t = pickle.load( open(pwd+'/output'+'{:04d}'.format(i)+'.p', "rb" ) )

ax = plt.gca()
save_frame(mesh,data.tracer,'Tracer',t,ax,pwd+'/tracer'+'{:04d}'.format(i)+'.png',False)
save_frame(mesh,0.5*(data.u_vel**2+data.v_vel**2),'Energy',t,ax,pwd+'/energy'+'{:04d}'.format(i)+'.png',False)
save_frame(mesh,data.press,'Pressure',t,ax,pwd+'/pressure'+'{:04d}'.format(i)+'.png',False)




