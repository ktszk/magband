import numpy as np
import matplotlib.pyplot as plt
nk=33
nw=151
data=np.loadtxt('fort.50')
x=data[:,0].reshape(nk,nw)
y=data[:,1].reshape(nk,nw)
z=data[:,2].reshape(nk,nw)
plt.contourf(x,y,z,100)
plt.jet()
plt.colorbar()
plt.show()
