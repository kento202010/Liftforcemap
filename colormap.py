#! /usr/bin/env python3
import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ptick
from matplotlib.colors import ListedColormap
from mpl_toolkits.axes_grid1 import make_axes_locatable

FILENAME = sys.argv[1]
N = 160
NODES = N+1

data = np.loadtxt(FILENAME, unpack=True)
Y = data[0].reshape(NODES, NODES)
Z = data[1].reshape(NODES, NODES)
F = data[2].reshape(NODES, NODES)

###axis
fig = plt.figure(figsize=(8, 7))
ax = fig.add_subplot(1, 1, 1)
ax.spines['top'].set_linewidth(3)
ax.spines['bottom'].set_color('grey')
ax.spines['left'].set_color('grey')
ax.spines['right'].set_linewidth(3)
###range
plt.xlim(-0.0525, 0.5)
plt.ylim(-0.0525, 0.5)
plt.xticks( [  0, 0.125, 0.25, 0.375,0.5] )
plt.yticks( [  0, 0.125, 0.25, 0.375,0.5] )
label_l = [ '$0$',  '$y$', '$0.25$', '','$0.5$']
label_m = [ '$0$',  '$z$', '$0.25$', '','$0.5$']
plt.rcParams["mathtext.fontset"] = 'cm'
plt.rcParams['mathtext.default'] = 'it'
ax.set_xticklabels(label_l, rotation=0,fontsize=15)
ax.set_yticklabels(label_m, rotation=90,fontsize=15)
### contour
cont = ax.contour(Y, Z, F, levels=[0.0], colors="#000000", linewidths=0.4)
cont.clabel(fmt="%1.1f", fontsize=9)
### color map
val1 = np.amin(F)
val2 = np.amax(F)
c = ax.pcolormesh(Y, Z, F, cmap="jet", vmin=val1, vmax=val2)
ax.set_aspect("equal")

### color map
val1 = np.amin(F)
val2 = np.amax(F)
c = ax.contourf(Y, Z, F,50, cmap="jet", vmin=val1, vmax=val2)
divider = make_axes_locatable(ax)
cax = divider.append_axes("top", size="2.7%", pad=0.3)
cbar = plt.colorbar(c, format="$%2.4le$",ticks=[val1, val2], orientation = "horizontal", cax=cax)
cax.xaxis.set_ticks_position('top')
ax.set_aspect("equal")
plt.tight_layout()

plt.savefig("output.ps")
#   plt.show()
