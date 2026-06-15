import numpy as np
import matplotlib.pyplot as plt

print("plotting mesh and phi and psi contours")

# Read mesh
with open("mesh.dat", "r") as f:
    Nx_mesh, Ny_mesh = map(int, f.readline().split())

mesh_data = np.loadtxt("mesh.dat", skiprows=1)

x_mesh = mesh_data[:,0]
y_mesh = mesh_data[:,1]

x_unique = np.unique(x_mesh)
y_unique = np.unique(y_mesh)

# Read phi field
# with open("phi_solved.dat") as f:
#     Nx, Ny = map(int, f.readline().split())

# data_phi = np.loadtxt("phi_solved.dat", skiprows=1)

# x   = data_phi[:,0].reshape(Ny, Nx)
# y   = data_phi[:,1].reshape(Ny, Nx)
# phi = data_phi[:,2].reshape(Ny, Nx)

#Read psi field
with open("psi_SOR_solved.dat") as f:
    Nx, Ny = map(int, f.readline().split())

data_psi = np.loadtxt("psi_SOR_solved.dat", skiprows=1)
x   = data_psi[:,0].reshape(Ny, Nx)
y   = data_psi[:,1].reshape(Ny, Nx)
psi = data_psi[:,2].reshape(Ny, Nx)


# Plot

plt.figure(figsize=(20,15))

# ---- Mesh ----
# for xi in x_unique:
#     plt.plot([xi, xi],[y_unique[0], y_unique[-1]],'k',linewidth=0.3,alpha=0.5)

# for yj in y_unique:
#     plt.plot([x_unique[0], x_unique[-1]],[yj, yj],'k',linewidth=0.3,alpha=0.5)

# ---- Iso-phi lines ----
plt.rcParams['contour.negative_linestyle'] = 'solid'

#phi_contours = plt.contour(x,y,phi,levels=100,colors="red",linewidths=0.7)
psi_contours = plt.contour(x,y,psi,levels=100,colors="blue",linewidths=0.7)

plt.contour(x, y, psi,levels=[-0.001],colors='blue',linewidths=0.7)

#plt.clabel(phi_contours,inline=True,fontsize=7)
#plt.clabel(psi_contours,inline=True,fontsize=7)

plt.xlabel("x")
plt.ylabel("y")
plt.axis("equal")

plt.savefig("mesh_with_phi_psi.png",dpi=600,bbox_inches='tight')

plt.show()