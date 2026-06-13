import numpy as np
import matplotlib.pyplot as plt

print("plotting convergence")

plt.figure(figsize=(20,10))

data_J = np.loadtxt("phi_J_converge.dat")
iter_J=data_J[:,0]
error_J=data_J[:,1]
iter_J_final = int(iter_J[-1])
plt.loglog(iter_J,error_J,'-',linewidth=1.5,label=f'Jacobi : {iter_J_final}')

data_GS = np.loadtxt("phi_GS_converge.dat")
iter_GS=data_GS[:,0]
error_GS=data_GS[:,1]
iter_GS_final = int(iter_GS[-1])
plt.loglog(iter_GS,error_GS,'-',linewidth=1.5,label=f'Gauss-Seidel : {iter_GS_final}')

data_SOR = np.loadtxt("phi_SOR_converge.dat")
iter_SOR=data_SOR[:,0]
error_SOR=data_SOR[:,1]
iter_SOR_final = int(iter_SOR[-1])
plt.loglog(iter_SOR,error_SOR,'-',linewidth=1.5,label=f'SOR (ω=1.5) : {iter_SOR_final}')

#data_BiCG= np.loadtxt("BiCG_phi_error_converge.dat")
#iter_BICG=data_BiCG[:,0]
#error_BiCG=data_BiCG[:,1]
#plt.loglog(iter_BICG,error_BiCG,'-',linewidth=2,label='BiCG')

plt.legend()
plt.title("PHI_solver convergence")
plt.xlabel("log(iteration)")
plt.ylabel("log(rms error)")
plt.grid(True)

plt.savefig("PHI_solver_convergence.png", dpi=800, bbox_inches='tight')

plt.show()