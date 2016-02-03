import numpy as np
from matplotlib import pyplot as plt
from matplotlib import animation
from flows import *
import pdb

def update_quiver(k, Q, X, Y, ax, ipts):
    """updates the horizontal and vertical vector components by a
    fixed increment on each frame
    """
    xy = [X, Y];
    uv = HamaVelocity(t[k], xy, input_args);
    U = uv[0];
    V = uv[1];

    # Update the plot
    Q.set_UVC(U,V)
    
    # Set some limits
    ax.set_xlim(-1, 9);
    ax.set_ylim(-0.15, 0.15);

    # Return the plot
    return Q,

# Make some coordinate vectors
xv = np.linspace(0, 8, 30);
yv = np.linspace(-0.1, 0.1, 20);

# Form the coordinates into a mesh
X, Y = np.meshgrid(xv, yv);

# Form the coordinates into a list
xy = [X, Y];

# Initial and final times
t0 = 0;
tf = 100;

# Number of time steps
n_steps = 500;

# Time vector
t = np.linspace(t0,tf, num = n_steps)

# Inputs to the hama flow are:
# amplitude a
# wave number alpha
# wave speed c
a = 0.2;
alpha = 1.0;
c = 1.0;

# Input arguments to the hama flow.
input_args = [a, alpha, c];

# Calculate the velocity
uv = HamaVelocity(t0, xy, input_args);
U = uv[0];
V = uv[1];

# Make a figure
fig, ax = plt.subplots(1,1)

# Make a quiver plot
Q = ax.quiver(X, Y, U, V, pivot='mid', color='k', units='inches', scale = 5)
# Q = ax.streamplot(X, Y, U, V, density=[0.5, 1])



# you need to set blit=False, or the first set of arrows never gets
# cleared on subsequent frames
anim = animation.FuncAnimation(fig, update_quiver,
                             frames = range(0, t.size),
                             fargs=(Q, X, Y, ax, input_args),
                               interval=5, blit=False)

plt.show(anim)


