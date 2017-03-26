import matplotlib.pyplot as plt
import numpy as np

data1 = np.loadtxt("v_rms32");
data2 = np.loadtxt("v_rms_bulk32");
data3 = np.loadtxt("v_rms64");
data4 = np.loadtxt("v_rms_bulk64");
data5 = np.loadtxt("v_rms64");

data6 = np.linspace(0,3,100);


fig, ax = plt.subplots()

ax.semilogy(data1[:,0],data1[:,1], label='$u_b = 0\ on\ 32 \\times32$')
ax.semilogy(data2[:,0],data2[:,1], label='$u_b=100\ on\ 32 \\times32$')
ax.semilogy(data3[:,0],data3[:,1], label='$u_b = 0\ on\ 64 \\times64$')
ax.semilogy(data4[:,0],data4[:,1], label='$u_b=100\ on\ 64 \\times64$')
ax.semilogy(data5[:,0],data5[:,1], label='$u_b=0\ on\ 128 \\times128$')
ax.semilogy(data6,np.exp(2*np.pi*data6), label='$ideal$')
legend  = ax.legend(loc = 'lower right')
plt.ylabel('growth rate')
plt.xlabel('time')
#plt.ylim((0.007,0.3))
plt.savefig('growth.png')
plt.show(); 
