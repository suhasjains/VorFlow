import matplotlib.pyplot as plt
import numpy as np

data = np.loadtxt("v_rms");

plt.hold(True)
plt.plot(data[:,0],data[:,1])
plt.ylabel('growth rate')
plt.xlabel('time')
plt.savefig('growth.png')
plt.show(); 

