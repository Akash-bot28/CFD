import numpy as np
import matplotlib.pyplot as plt

print(f"Plotting mesh...")

data=np.loadtxt(f"mesh.dat")

x=data[:,0]
y=data[:,1]

x_unique = np.unique(x)
y_unique = np.unique(y)

plt.figure(figsize=(20,15))

# verical lines grids
for xi in x_unique:
    plt.plot([xi, xi], [y_unique[0],y_unique[-1]], 'k', linewidth=0.5, alpha=0.7)

#horizontal lines grids
for yj in y_unique:
    plt.plot([x_unique[0],x_unique[-1]], [yj, yj], 'k', linewidth=0.5, alpha=0.7)

#block interface plotting

x_interfaces = []
y_interfaces = []

section = None

with open("interfaces.dat") as f:
    for line in f:
        line = line.strip()
        if not line: continue
        if line == "# X":
            section = "X"
            continue
        if line == "# Y":
            section = "Y"
            continue
        if section == "X":
            x_interfaces.append(float(line))
        elif section == "Y":
            y_interfaces.append(float(line))

for xi in x_interfaces:
    plt.plot([xi, xi],[y_unique[0], y_unique[-1]],color='red',linewidth=0.9,alpha=0.7)

for yi in y_interfaces:
    plt.plot([x_unique[0], x_unique[-1]],[yi, yi],color='blue',linewidth=0.9,alpha=0.7)

plt.axis('equal')
plt.xlabel('x')
plt.ylabel('y')
plt.grid(False)

plt.savefig(f"mesh.png", dpi=1000, bbox_inches='tight')

plt.show()