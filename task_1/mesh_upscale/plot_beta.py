
import numpy as np
import matplotlib.pyplot as plt

print(f"Beta optimizing for mesh...")

xdata = np.loadtxt(f"X_optimization.dat")
ydata = np.loadtxt(f"Y_optimization.dat")

iter_x  = xdata[:,0]
error_x = xdata[:,1]

iter_y  = ydata[:,0]
error_y = ydata[:,1]

fig, ax = plt.subplots(2, 1, figsize=(12,10))

ax[0].plot(iter_x,error_x,label="Maximum Interface Error",linewidth=3)


ax[0].set_title("X-axis beta optimization")
ax[0].set_xlabel("Iteration")
ax[0].set_ylabel("Max_beta_error")
ax[0].grid(True)
ax[0].legend()



ax[1].plot(iter_y,error_y,label="Maximum Interface Error",linewidth=3)


ax[1].set_title("Y-axis beta optimization")
ax[1].set_xlabel("Iteration")
ax[1].set_ylabel("Max_beta_error")
ax[1].grid(True)
ax[1].legend()



plt.tight_layout()
plt.savefig("beta_optimization.png", dpi=1000, bbox_inches='tight')

plt.show()