import numpy as np
import matplotlib.pyplot as plt

data=np.loadtxt("mesh_10.dat")

x=data[:,0]
y=data[:,1]

x_unique = np.unique(x)
y_unique = np.unique(y)

XC = 0.0
YC = 0.0

a  = 1.0

# interface locations

XSL = XC - a/2
XSR = XC + a/2

YSB = YC - a/2
YST = YC + a/2

plt.figure(figsize=(20,15))

# verical lines
for xi in x_unique:
    plt.plot([xi, xi], [y_unique[0],y_unique[-1]], 'k', linewidth=0.5, alpha=0.7)

#horizontal lines
for yj in y_unique:
    plt.plot([x_unique[0],x_unique[-1]], [yj, yj], 'k', linewidth=0.5, alpha=0.7)

plt.plot([XSL, XSL], [y_unique[0],y_unique[-1]], color='red', linewidth=0.9, alpha=0.7)
plt.plot([XSR, XSR], [y_unique[0],y_unique[-1]], color='red', linewidth=0.9, alpha=0.7)
plt.plot([x_unique[0],x_unique[-1]], [YSB, YSB], color='blue', linewidth=0.9, alpha=0.7)
plt.plot([x_unique[0],x_unique[-1]], [YST, YST], color='blue', linewidth=0.9, alpha=0.7)


plt.axis('equal')
plt.xlabel('x')
plt.ylabel('y')
plt.grid(False)

plt.savefig("mesh_10.png", dpi=1000, bbox_inches='tight')

plt.show()