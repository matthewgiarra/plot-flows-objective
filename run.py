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

# Flow type
flow_type = "womersley"

if flow_type == "hama":
	# Domain
	xd = (-1.1, 24)
	yd = (-3, 3)
	zd = (-1, 1)

    # Hama flow constants
	a = 0.05
	alpha = 1
	c = 1
	extra_args = [a, alpha, c]
elif flow_type == "poiseuille":
	# Domain
	xd = (-1.1, 10)
	yd = (-3, 3)
	zd = (-1, 1)
	
	# Velocity Constant
	a = 5
	
	extra_args = [a]
elif flow_type == "uniform":
	# Domain
	xd = (-1.1, 10)
	yd = (-3, 3)
	zd = (-1, 1)
	
	# Velocity Constant
	a = 5

	extra_args = [a]
elif flow_type == "womersley":
	# Domain
	xd = (-5, 5)
	yd = (-3, 3)
	zd = (-1, 1)

	#Womersley flow constants
	
	#Womersley Number
	Wn = 0.1
	
	# Complex Number
	c = 12 + 30j
				
	#Frequency 
	f = 7
	
	extra_args = [Wn, c, f]
elif flow_type == "oscillatingplane":

    # Domain
    xd = (-15, 15)
    yd = (-3, 3)
    zd = (-1, 1)
	
	# Frequency
    f = 4 
	
    extra_args = [f]

	
# Plot type
plot_type = "streak"

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
