import numpy as np
import matplotlib.pyplot as plt

print("plotting convergence")

plt.figure(figsize=(20,10))

#data_GS = np.loadtxt("GS_phi_error_converge.dat")
#iter_GS=data_GS[:,0]
#error_GS=data_GS[:,1]
#plt.loglog(iter_GS,error_GS,'-',linewidth=2,label='Gauss-Seidel')

data_SOR = np.loadtxt("phi_error_converge.dat")
iter_SOR=data_SOR[:,0]
error_SOR=data_SOR[:,1]
plt.loglog(iter_SOR,error_SOR,'-',linewidth=2,label='SOR (ω=1.5)')

#data_BiCG= np.loadtxt("BiCG_phi_error_converge.dat")
#iter_BICG=data_BiCG[:,0]
#error_BiCG=data_BiCG[:,1]
#plt.loglog(iter_BICG,error_BiCG,'-',linewidth=2,label='BiCG')

plt.legend()
plt.title("solver convergence")
plt.xlabel("log(iteration)")
plt.ylabel("log(rms error)")
plt.grid(True)

plt.savefig("solver_convergence.png", dpi=800, bbox_inches='tight')

plt.show()