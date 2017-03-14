import numpy as np
import sys
import os
sys.path.append('../../../')
from plotting import *
import pickle

pwd = os.getcwd()

Nbin = 100
Nt = 10
veltot = np.zeros((Nbin,Nt))
ttot = np.zeros(Nt)
ybin = np.linspace(0,1,Nbin+1)
ybinint = 0.5 * (ybin[1:]+ybin[:-1])
ymidtot = np.zeros(Nt)

for i in range(Nt):

	print i
	_, _, t, velavg = pickle.load( open(pwd+'/output'+'{:04d}'.format(i)+'.p', "rb" ) )
	velbin = velavg.flatten()
	ttot[i] = t
	veltot[:,i] = velbin
	if not i == 0:
		velext = 0.05*(np.max(velbin)-np.min(velbin))+np.min(velbin)
		yext = np.interp(velext,np.flipud(velbin[Nbin/2:]),np.flipud(ybinint[Nbin/2:]))
		#velmid = 0.5*(np.max(velbin)-np.min(velbin))+np.min(velbin)
		#ymid = np.interp(velmid,np.flipud(velbin[Nbin/2:]),np.flipud(ybinint[Nbin/2:]))
		#print yext,ymid
		ymidtot[i] = yext - 0.5

print ymidtot
print np.sqrt(ttot)
