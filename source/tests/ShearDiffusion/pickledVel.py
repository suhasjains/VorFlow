import numpy as np
import sys
import os
sys.path.append('../../../')
from plotting import *
import pickle
import matplotlib.pyplot as plt

pwd = os.getcwd()

Nbin = 100
Nt = 50
veltot = np.zeros((Nbin,Nt))
ttot = np.zeros(Nt)
ybin = np.linspace(0,2,Nbin+1)
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
		ymidtot[i] = (yext - 1.0 - 1.0/12.0)/2.0 #first deduction to midline; second deduction to edge of shear layer from midline

print ymidtot
print np.sqrt(ttot)

#ttot=[0,0.1*0.1,0.2*0.2]
#ymidtot=[0,0.07848,0.166835]

plt.plot(np.sqrt(ttot),ymidtot,'-o',ms=10)
plt.xlabel(r'$\sqrt{t} [\sqrt{s}]$')
plt.ylabel('Nondimensional half-width of central region')
plt.title('Diffusion of shear layer, Re = 10')
plt.savefig('plot1'+'{:04d}'.format(Nt)+'.png')

plt.figure()
plt.plot(np.log(ttot[1:]),np.log(ymidtot[1:]),'-o',ms=10)
plt.xlabel(r'$\log{}(t[s])$')
plt.ylabel(r'$\log{}($Half-width$)$')
plt.title('Diffusion of shear layer, Re = 10')
plt.savefig('plot2'+'{:04d}'.format(Nt)+'.png')
