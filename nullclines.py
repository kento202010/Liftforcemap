#! /usr/bin/env python3
import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ptick

FILENAME1 = "IntpltdFr.dat"
FILENAME2 = "IntpltdFt.dat"
FILENAME3 = "Center.dat"
FILENAME4 = "Stable_EquilibriumPoints.dat"
FILENAME5 = "Saddle_EquilibriumPoints.dat"
FILENAME6 = "Lift_Direction.dat"

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
y = data[0]
z = data[1]
r = data[2]
Fy = data[3]
Fz = data[4]
Fr_moto = data[5]
Fz_moto = data[6]


fig = plt.figure(figsize=(8, 7))

###axis
ax = fig.add_subplot(1, 1, 1)
ax.spines['top'].set_linewidth(4)
ax.spines['bottom'].set_color('white')
ax.spines['left'].set_color('white')
ax.spines['right'].set_linewidth(4)
###range
plt.xlim(-0.0525, 0.5)
plt.ylim(-0.0525, 0.5)
plt.xticks( [  0, 0.125, 0.25, 0.375,0.5] )
plt.yticks( [  0, 0.125, 0.25, 0.375,0.5] )
label_l = [ '$0$',  '$y$', '$0.25$', '','$0.5$']
label_m = [ '$0$',  '$z$', '$0.25$', '','$0.5$']
plt.rcParams["mathtext.fontset"] = 'cm'
plt.rcParams['mathtext.default'] = 'it'

ax.set_xticklabels(label_l, rotation=0,fontsize=20)
ax.set_yticklabels(label_m, rotation=90,fontsize=20)
ax.tick_params(direction="in", length=length, width=width)
ax.xaxis.tick_top()
ax.yaxis.tick_right()
ax.xaxis.set_label_position("top")
ax.yaxis.set_label_position("right")

### contour
cont = ax.contour(Y, Z, Fr, levels=[0.0], colors="black", linewidths=5,zorder=1)
cont = ax.contour(Y, Z, Ft, levels=[0.0], colors="black",linestyles="dashed",linewidths=3,zorder=1)
### points
plt.scatter(Ye0, Ze0, marker="s", s=200, c="white",zorder=2,edgecolor="black", linewidths=3)
ax.set_aspect("equal")
plt.scatter(Ye1, Ze1, marker="o", s=200, c="red",zorder=2,edgecolor="black", linewidths=3)
ax.set_aspect("equal")
plt.scatter(Ye2, Ze2, marker="o", s=200, c="white",zorder=2,edgecolor="black", linewidths=3)
ax.set_aspect("equal")

###bector
R = np.sqrt(Fy ** 2 + Fz ** 2)
plt.quiver(y, z, Fy/R, Fz/R, R, cmap="jet")
plt.quiver(z, y, Fz/R, Fy/R, R, cmap="jet")
plt.quiver(y, z, Fy/R, Fz/R, edgecolor="black",facecolor="None",linewidth=0.5)
plt.quiver(z, y, Fz/R, Fy/R, edgecolor="black",facecolor="None",linewidth=0.5)

plt.savefig("output.ps")
#   plt.show()

