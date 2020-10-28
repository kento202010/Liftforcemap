#! /usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ptick

FILENAME1 = "IntpltdFr.dat"
FILENAME2 = "IntpltdFt.dat"
FILENAME3 = "Stable_EquilibriumPoints.dat"
FILENAME4 = "Saddle_EquilibriumPoints.dat"
FILENAME5 = "data.dat"
FILENAME6 = "Center.dat"
N = 150
NODES = N+1

data = np.loadtxt(FILENAME1, unpack=True)
Y = data[0].reshape(NODES, NODES)
Z = data[1].reshape(NODES, NODES)
Fr = data[2].reshape(NODES, NODES)
data = np.loadtxt(FILENAME2, unpack=True)
Ft = data[2].reshape(NODES, NODES)
data = np.loadtxt(FILENAME3, unpack=True)
Ye = data[0]
Ze = data[1]
data = np.loadtxt(FILENAME4, unpack=True)
Ye1 = data[0]
Ze1 = data[1]
data = np.loadtxt(FILENAME5, unpack=True)
y = data[0]
z = data[1]
Fy = data[2]
Fz = data[3]
data = np.loadtxt(FILENAME6, unpack=True)
Ye2 = data[0]
Ze2 = data[1]

fig = plt.figure(figsize=(8, 7))

### axis
ax = fig.add_subplot(1, 1, 1)
ax.spines['top'].set_linewidth(3)
ax.spines['bottom'].set_color('grey')
ax.spines['left'].set_color('grey')
ax.spines['right'].set_linewidth(3)

### range
plt.xlim(0, 0.5)
plt.ylim(0, 0.5)
plt.xticks( [ -0.0625, 0, 0.125, 0.25, 0.375,0.5] )
plt.yticks( [ -0.0625, 0, 0.125, 0.25, 0.375,0.5] )
label_l = [' ' , '0',  '', '0.25', '','0.5']
ax.set_xticklabels(label_l, rotation=0, fontname='Arial', fontsize=15)
ax.set_yticklabels(label_l, rotation=90, fontname='Arial', fontsize=15)

### contour
cont = ax.contour(Y, Z, Fr, levels=[0.0], colors="#000000", linestyles="dashed",zorder=1)
cont = ax.contour(Y, Z, Ft, levels=[0.0], colors="#000000",zorder=1)

### points
plt.scatter(Ye, Ze, marker="o", s=200, c="red",zorder=2,linewidths="2",edgecolor="black")
ax.set_aspect("equal")
plt.scatter(Ye1, Ze1, marker="o", s=200, c="white",zorder=2,linewidths="2",edgecolor="black")
ax.set_aspect("equal")
plt.scatter(Ye2, Ze2, marker="s", s=200, c="white",zorder=2,linewidths="2",edgecolor="black")
ax.set_aspect("equal")

### marker
ax.scatter(y, z,marker = ".", s=50, c= "black")

### bector
R = np.sqrt(Fy ** 2 + Fz ** 2)
im=plt.quiver(y, z, Fy/R, Fz/R, R, cmap="jet", width = 0.006)
plt.colorbar(im)

### border
xmin, xmax = -0.425, 0.425
ymin, ymax = -0.425, 0.425
plt.plot([xmin, xmax],[0.425, 0.425], "lightgrey", linestyle='solid')
plt.plot([0.425, 0.425],[ymin, ymax], "lightgrey", linestyle='solid')

### output
plt.savefig("out/nullcline.eps")
plt.savefig("out/nullcline.pdf")
plt.show()