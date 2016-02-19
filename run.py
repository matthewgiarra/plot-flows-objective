import flows
import numpy as np
import pdb

num_particles_y = 7
num_particles_x = 1

# Coordinates
xv = np.linspace(0, 0, num_particles_x)
yv = np.linspace(-1, 1, num_particles_y)

xm, ym = np.meshgrid(xv, yv)

x = xm.flatten()
y = ym.flatten()

# Z values
z = np.zeros(len(x))

# Make them a vector
y0 = [x, y, z]

a = 0.05
alpha = 1
c = 1

# Domain
xd = (-1.1, 24.0)
yd = (-5, 5)
zd = (-1, 1)

# Hama flow constants
extra_args = [a, alpha, c]

# Plot type
plot_type = "streak"

# Flow type
flow_type = "hama"

# New particle distance
NewParticleDistance = 0.1

# Time stuff
sim_time = [0, 50, 0.1]

sim = flows.Simulation(x_domain=xd,
                       y_domain=yd,
                       z_domain=zd,
                       y0=y0,
                       flow_type=flow_type,
                       plot_type=plot_type,
                       time=sim_time,
                       NewParticleDistance=NewParticleDistance,
                       extra_args=extra_args)

# Run the simulation                                
sim.Run()
