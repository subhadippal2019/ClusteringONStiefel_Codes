import matplotlib.pyplot as plt
import numpy as np

data = np.genfromtxt('density1.txt',delimiter='\t')

x=data[:,0]
y=data[:,1]
z=data[:,2]

## Equivalently, we could do that all in one line with:
# x,y,z = np.genfromtxt('eye_.txt', delimiter=',', usecols=(0,1,2))

shape = np.unique(x).shape[0],np.unique(y).shape[0]
x_arr = x.reshape(shape)
y_arr = y.reshape(shape)
z_arr = z.reshape(shape)
plt.pcolormesh(x_arr,y_arr,z_arr)

plt.show()
