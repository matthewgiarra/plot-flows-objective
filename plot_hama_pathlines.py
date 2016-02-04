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
# Number of points
num_points = 10;
xv = np.linspace(0, 0, num_points);
yv = np.linspace(-1, 1, num_points);

# Form the coordinates into a list
y0 = np.append(xv, yv);
#
# y0 = np.meshgrid(xv, yv)

# Initial and final times
t0 = 0;
tf = 10;
dt = 0.1;

# Number of time steps
n_steps = 100;

# Time vector
t = np.linspace(t0, tf, num = n_steps)

# Inputs to the hama flow are:
# amplitude a
# wave number alpha
# wave speed c
a = 0.05;
alpha = 1.0;
c = 1;

# Input arguments to the hama flow.
input_args_hama = [a, alpha, c];

vel = 1;

input_args_uniform = vel;

m = Particle();
m.Update("hama", t0 = t0, dt = dt, input_args = input_args_hama);

pdb.set_trace();

print("Hellow!")

# uv = UniformVelocity(xy = y0, t = t[0], my_args = vel);


# pdb.set_trace();

# Plot the pathlines
# plot_pathlines(flow_type = "uniform", t = t, y0 = y0, my_args = input_args_uniform);
# plot_pathlines(flow_type = "hama", t = t, y0 = y0, my_args = input_args_hama);
# uv = HamaVelocity(xy = y0, t = t[0], my_args = [a, alpha, c]);










