import sys
import numpy as np
import matplotlib.pyplot as plt

mesh_id = sys.argv[1]

xdata = np.loadtxt(f"X_optimization_{mesh_id}.dat")
ydata = np.loadtxt(f"Y_optimization_{mesh_id}.dat")

iter_x   = xdata[:,0]
beta1_x  = xdata[:,1]
betas_x  = xdata[:,2]
beta2_x  = xdata[:,3]
error1_x = xdata[:,4]
error2_x = xdata[:,5]

iter_y   = ydata[:,0]
beta1_y  = ydata[:,1]
betas_y  = ydata[:,2]
beta2_y  = ydata[:,3]
error1_y = ydata[:,4]
error2_y = ydata[:,5]

fig, ax = plt.subplots(2, 1, figsize=(12,10))

#ax[0].scatter(iter_x,beta1_x,label="beta1",s=12)
#ax[0].scatter(iter_x,betas_x,label="betas",s=12)
#ax[0].scatter(iter_x,beta2_x,label="beta2",s=12)
ax[0].plot(iter_x,error1_x,label="error1",linewidth=3)
ax[0].plot(iter_x,error2_x,label="error2",linewidth=3)

ax[0].set_title("X-axis beta optimization")
ax[0].set_xlabel("Iteration")
#ax[0].set_ylabel("Beta")
ax[0].grid(True)
ax[0].legend()


#ax[1].scatter(iter_y,beta1_y,label="beta1",s=12)
#ax[1].scatter(iter_y,betas_y,label="betas",s=12)
#ax[1].scatter(iter_y,beta2_y,label="beta2",s=12)
ax[1].plot(iter_y,error1_y,label="error1",linewidth=3)
ax[1].plot(iter_y,error2_y,label="error2",linewidth=3)

ax[1].set_title("Y-axis beta optimization")
ax[1].set_xlabel("Iteration")
#ax[1].set_ylabel("Beta")
ax[1].grid(True)
ax[1].legend()



plt.tight_layout()
plt.savefig("beta_optimization_{mesh_id}.png", dpi=1000, bbox_inches='tight')

plt.show()