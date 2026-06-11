import numpy as np
import matplotlib.pyplot as plt

print("plotting convergence")

data = np.loadtxt("phi_error_converge.dat")

iter=data[:,0]
error=data[:,1]

plt.figure(figsize=(20,10))

plt.scatter(iter,error,s=5)

plt.axis('equal')
plt.title("solver convergence")
plt.xlabel("iteration")
plt.ylabel("rms error")
plt.grid(True)

plt.savefig("solver_convergence.png", dpi=800, bbox_inches='tight')

plt.show()