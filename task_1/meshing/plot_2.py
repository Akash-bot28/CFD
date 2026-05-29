import numpy as np
import matplotlib.pyplot as plt


data=np.loadtxt("mesh_7.dat")
x=data[:,0]
y=data[:,1]
s=data[:,2].astype(bool)

plt.figure(figsize=(15,10))

# # verical lines
# for i in range(Nx):

#     for j in range(Ny-1):

#         p1 = j*Nx + i
#         p2 = (j+1)*Nx + i

#         if not s[p1] and not s[p2]:

#             plt.plot([x[p1], x[p2]], [y[p1], y[p2]],'k',linewidth=0.8)

# #horizontal lines
# for j in range(Ny):

#     for i in range(Nx-1):

#         p1 = j*Nx + i
#         p2 = j*Nx + (i+1)

#         if not s[p1] and not s[p2]:

#             plt.plot(
#                 [x[p1], x[p2]],[y[p1], y[p2]],'k',linewidth=0.8)

# verical lines
for xi in x:
    plt.plot([xi, xi], [y[0],y[-1]], 'k', linewidth=0.8)

#horizontal lines
for yj in y:
    plt.plot([x[0],x[-1]], [yj, yj], 'k', linewidth=0.8)    

plt.axis('equal')
plt.xlabel('x')
plt.ylabel('y')
plt.grid(False)

plt.savefig("mesh_7.png", dpi=300, bbox_inches='tight')

#plt.show()