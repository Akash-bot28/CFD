import numpy as np
import matplotlib.pyplot as plt

file = np.loadtxt("phi solved.dat", skiprows=1)

x   = file[:,0]
y   = file[:,1]
phi = file[:,2]

Nx = 80
Ny = 60

X   = x.reshape(Ny, Nx)
Y   = y.reshape(Ny, Nx)
PHI = phi.reshape(Ny, Nx)

plt.figure(figsize=(8,5))

plt.contour(X, Y, PHI, levels=30)
plt.xlabel("x")
plt.ylabel("y")
plt.axis("equal")

plt.savefig(f"iso_phi.png", dpi=900, bbox_inches='tight')

plt.show()