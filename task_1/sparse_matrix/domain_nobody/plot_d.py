import numpy as np
import matplotlib.pyplot as plt


print(f"Plotting grid spacing for mesh ...")

data_dx = np.loadtxt(f"dx.dat")
data_dy = np.loadtxt(f"dy.dat")

x  = data_dx[:,0]
dx = data_dx[:,1]

y  = data_dy[:,0]
dy = data_dy[:,1]

fig, axes = plt.subplots(2,1,figsize=(14,10))


axes[0].scatter(x,dx,color='black',s=12)

axes[0].set_xlabel('x')

axes[0].set_ylabel('dx')

axes[0].set_title('dx vs x')

axes[0].grid(True)


axes[1].scatter(y,dy,color='black',s=12)

axes[1].set_xlabel('y')

axes[1].set_ylabel('dy')

axes[1].set_title('dy vs y')

axes[1].grid(True)


plt.tight_layout()

plt.savefig(f"spacing.png",dpi=500,bbox_inches='tight')

plt.show()