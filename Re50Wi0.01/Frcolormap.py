#! /usr/bin/env python3
import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ptick

FILENAME = "IntpltdFr.dat"
N = 150
NODES = N+1

data = np.loadtxt(FILENAME, unpack=True)
Y = data[0].reshape(NODES, NODES)
Z = data[1].reshape(NODES, NODES)
F = data[2].reshape(NODES, NODES)

fig = plt.figure(figsize=(8, 7))

### ax
ax = fig.add_subplot(1, 1, 1)
ax.spines['top'].set_linewidth(3)
ax.spines['bottom'].set_linewidth(3)
ax.spines['left'].set_linewidth(3)
ax.spines['right'].set_linewidth(3)

### range
plt.xlim(-0.5, 0.5)
plt.ylim(-0.5, 0.5)
plt.xticks( [-0.5, -0.4375, -0.375, -0.3125,  -0.25, -0.1875, -0.125, -0.0625, 0, 0.0625, 0.125, 0.1875, 0.25, 0.3125, 0.375, 0.4375,0.5] )
plt.yticks( [-0.5, -0.4375, -0.375, -0.3125,  -0.25, -0.1875, -0.125, -0.0625, 0, 0.0625, 0.125, 0.1875, 0.25, 0.3125, 0.375, 0.4375,0.5] )
label_l = ['-0.5', '', '', '', '-0.25', '', '', '','0', '', '', '', '0.25', '', '', '','0.5']
ax.set_xticklabels(label_l, rotation=0, fontname='Arial', fontsize=15)
ax.set_yticklabels(label_l, rotation=90, fontname='Arial', fontsize=15)

### contour
cont = ax.contour(Y, Z, F, levels=[0.0], colors="#000000", linewidths=1.5)

### color map
val1 = np.amin(F)
val2 = np.amax(F)
c = ax.pcolormesh(Y, Z, F, cmap="jet", vmin=val1, vmax=val2)
plt.colorbar(c, format="%2.2lf")
ax.set_aspect("equal")

### output
plt.savefig("out/Frcolormap.eps")
plt.savefig("out/Frcolormap.pdf")
plt.show()