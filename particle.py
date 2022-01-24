#! /usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ptick
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

FILENAME1 = "IntpltdFr.dat"
FILENAME2 = "IntpltdFt.dat"
FILENAME3 = "Center.dat"
FILENAME4 = "Stable_EquilibriumPoints.dat"
FILENAME5 = "Saddle_EquilibriumPoints.dat"
FILENAME6 = "particle.dat"

#メモリの長さ
length = 15
#メモリの太さ
width = 4

N = 160
NODES = N+1

data = np.loadtxt(FILENAME1, unpack=True)
Y = data[0].reshape(NODES, NODES)
Z = data[1].reshape(NODES, NODES)
Fr = data[2].reshape(NODES, NODES)
data = np.loadtxt(FILENAME2, unpack=True)
Ft = data[2].reshape(NODES, NODES)
data = np.loadtxt(FILENAME3, unpack=True)
Ye0 = data[0]
Ze0 = data[1]
data = np.loadtxt(FILENAME4, unpack=True)
Ye1 = data[0]
Ze1 = data[1]
data = np.loadtxt(FILENAME5, unpack=True)
Ye2 = data[0]
Ze2 = data[1]
data = np.loadtxt(FILENAME6, unpack=True)
Ye3 = data[0]
Ze3 = data[1]

fig = plt.figure(figsize=(8, 7))

###axis
ax = fig.add_subplot(1, 1, 1)
ax.spines['top'].set_linewidth(4)
ax.spines['bottom'].set_linewidth(4)
ax.spines['left'].set_linewidth(4)
ax.spines['right'].set_linewidth(4)
###range
plt.xlim(-0.5, 0.5)
plt.ylim(-0.5, 0.5)
plt.xticks( [-0.5,  -0.25,  0, 0.25, 0.5] )
plt.yticks( [-0.5,  -0.25,  0, 0.25, 0.5] )
plt.rcParams["mathtext.fontset"] = 'cm'
plt.rcParams['mathtext.default'] = 'it'

label_l = [' ', ' ',' ', ' ',' ']
label_m = ['$-D/2$', ' ','$0$', ' ','$D/2$']
ax.set_xticklabels(label_l, rotation=0,fontsize=20)
ax.set_yticklabels(label_m, rotation=0,fontsize=20)
ax.tick_params(direction="in", length=length, width=width)
#ax.xaxis.tick_top()
#ax.yaxis.tick_right()
ax.xaxis.set_ticks_position('both')
ax.yaxis.set_ticks_position('both')
#ax.xaxis.set_label_position("top")
#ax.yaxis.set_label_position("left")


### contour
cont = ax.contour(Y, Z, Fr, levels=[0.0], colors="black", linewidths=2,zorder=1)
cont = ax.contour(Y, Z, Ft, levels=[0.0], colors="black",linestyles="dashed",linewidths=2,zorder=1)
### points

plt.scatter(Ye3/2, Ze3/2, marker="o", s=10, c="red",zorder=0, linewidths=1)
ax.set_aspect("equal")

plt.savefig("output.ps")
#plt.show()