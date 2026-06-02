import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

print("Plotting geometry preview...")

with open("preview.dat", "r") as file:

    # Domain
    xmin, xmax, ymin, ymax = map(float,file.readline().split())
    # Number of bodies
    nBodies = int(file.readline())
    bodies = []
    for _ in range(nBodies):
        bxmin, bxmax, bymin, bymax = map(float,file.readline().split())

        bodies.append((bxmin, bxmax, bymin, bymax))

fig, ax = plt.subplots(figsize=(10, 6))


# Block interfaces
x_interfaces = []
y_interfaces = []

for bxmin, bxmax, bymin, bymax in bodies:
    x_interfaces.append(bxmin)
    x_interfaces.append(bxmax)
    y_interfaces.append(bymin)
    y_interfaces.append(bymax)

x_interfaces = sorted(set(x_interfaces))
y_interfaces = sorted(set(y_interfaces))

# Domain rectangle
domain = Rectangle((xmin, ymin),xmax - xmin,ymax - ymin,fill=False,linewidth=1.5)
ax.add_patch(domain)

# Interface lines
#for x in x_interfaces:
#    ax.axvline(x,color='red',linewidth=2)

for x in x_interfaces:
    ax.plot([x, x],[ymin, ymax],color='red',linewidth=0.9)

#for y in y_interfaces:
#    ax.axhline(y,color='blue',linewidth=2)

for y in y_interfaces:
    ax.plot([xmin, xmax],[y, y],color='blue',linewidth=0.9)

# Bodies

for bxmin, bxmax, bymin, bymax in bodies:
    body = Rectangle((bxmin, bymin),bxmax - bxmin,bymax - bymin,fill=True,alpha=0.4)
    ax.add_patch(body)

ax.set_title("Geometry Preview")
ax.set_xlabel("x")
ax.set_ylabel("y")
ax.set_aspect("equal")

margin_x = 0.05 * (xmax - xmin)
margin_y = 0.05 * (ymax - ymin)

ax.set_xlim(xmin - margin_x,xmax + margin_x)
ax.set_ylim(ymin - margin_y,ymax + margin_y)
ax.grid(True)
plt.tight_layout()
plt.savefig("preview.png",dpi=600,bbox_inches="tight")

plt.show()