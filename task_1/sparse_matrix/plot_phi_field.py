import numpy as np
import matplotlib.pyplot as plt

with open("phi field.dat") as f:
    Nx, Ny = map(int, f.readline().split())

data = np.loadtxt("phi solved.dat", skiprows=1)

x   = data[:,0].reshape(Ny, Nx)
y   = data[:,1].reshape(Ny, Nx)
phi = data[:,2].reshape(Ny, Nx)

plt.figure(figsize=(8,6))

plt.contourf(x, y, phi, levels=50)
plt.colorbar(label='Phi')

plt.xlabel('x')
plt.ylabel('y')
plt.axis('equal')

plt.show()