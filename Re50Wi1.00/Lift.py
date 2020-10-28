#! /usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ptick

FILENAME5 = "data.dat"

data = np.loadtxt(FILENAME5, unpack=True)
y = data[0]
z = data[1]
Fy = data[2]
Fz = data[3]

fig = plt.figure(figsize=(8, 7))

###axis
ax = fig.add_subplot(1, 1, 1)
ax.spines['top'].set_linewidth(3)
ax.spines['bottom'].set_color('grey')
ax.spines['left'].set_color('grey')
ax.spines['right'].set_linewidth(3)

#ax.axvspan(-0.425, 0.425 ,color = "lightgrey",zorder=1)

###range
plt.xlim(0, 0.5)
plt.ylim(0, 0.5)
plt.xticks( [ -0.0625, 0, 0.125, 0.25, 0.375,0.5] )
plt.yticks( [ -0.0625, 0, 0.125, 0.25, 0.375,0.5] )
label_l = [' ' , '0',  '', '0.25', '','0.5']
ax.set_xticklabels(label_l, rotation=0, fontname='Arial', fontsize=15)
ax.set_yticklabels(label_l, rotation=90, fontname='Arial', fontsize=15)

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


plt.savefig("out/Lift.pdf")
plt.savefig("out/Lift.eps")
plt.show()