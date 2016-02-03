from flows import *
import pdb
import numpy as np;
from matplotlib import pyplot as plt
from matplotlib import animation

# tv = np.linspace(0, 10, 100);

xv = np.linspace(0, 1, 100);
yv = np.linspace(-1, 1, 10);

xm, ym = np.meshgrid(xv, yv);

xy = [xm, ym];

freq = 1.0;
amp = 1.0;

args1 = [amp, freq];

# for k in tv:
uv = PoiseuilleFlow(0, xy, args1);

u = uv[0];
v = uv[1];

fig, ax = plt.subplots(1,1);
Q = ax.quiver(xm, ym, u, v);

#     plt.quiver(xm, ym, u, v);
#     plt.show();
#



def update_quiver(t, Q, xy, args1):
    uv = PoiseuilleFlow(t, xy, args1);
    u = uv[0];
    v = uv[1];
    Q.set_UVC(u, v);
    return Q,
    
anim = animation.FuncAnimation(fig, update_quiver, fargs = (Q, xy, args1), interval = 10, blit = False)

plt.show()
