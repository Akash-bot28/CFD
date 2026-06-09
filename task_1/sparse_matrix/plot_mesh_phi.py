import numpy as np
import matplotlib.pyplot as plt

print("plotting mesh and iso phi contours")
# ==================================================
# Read mesh
# ==================================================
with open("mesh.dat", "r") as f:
    Nx_mesh, Ny_mesh = map(int, f.readline().split())

mesh_data = np.loadtxt("mesh.dat", skiprows=1)

x_mesh = mesh_data[:,0]
y_mesh = mesh_data[:,1]

x_unique = np.unique(x_mesh)
y_unique = np.unique(y_mesh)

# ==================================================
# Read phi field
# ==================================================
with open("phi_solved.dat") as f:
    Nx, Ny = map(int, f.readline().split())

data = np.loadtxt("phi_solved.dat", skiprows=1)

x   = data[:,0].reshape(Ny, Nx)
y   = data[:,1].reshape(Ny, Nx)
phi = data[:,2].reshape(Ny, Nx)


# ==================================================
# Plot
# ==================================================
plt.figure(figsize=(20,15))

# ---- Mesh ----
for xi in x_unique:
    plt.plot([xi, xi],[y_unique[0], y_unique[-1]],'k',linewidth=0.3,alpha=0.5)

for yj in y_unique:
    plt.plot([x_unique[0], x_unique[-1]],[yj, yj],'k',linewidth=0.3,alpha=0.5)

# ---- Iso-phi lines ----
contours = plt.contour(x, y, phi,levels=30,linewidths=1.0)

#plt.clabel(contours, inline=True, fontsize=8)

plt.xlabel("x")
plt.ylabel("y")
plt.axis("equal")

plt.savefig("mesh_with_iso_phi.png",dpi=600,bbox_inches='tight')

plt.show()