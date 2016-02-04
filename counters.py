import flows
import numpy as np
import pdb

num_particles = 20;

# Coordinates
x = list(np.linspace(0, 0, num_particles));
y = list(np.linspace(-1, 1, num_particles));
z = list(np.linspace(0, 0, num_particles));

# Make them a vector
y0 = [x, y, z];

a = 0.02;
alpha = 1;
c = 1;

# Domain
xd = (-1.1, 1.1);
yd = (-1.1, 1.1);
zd = (-1, 1);

# Hama flow constants
extra_args = [a, alpha, c];

# Plot type
plot_type = "streak";

# Flow type
flow_type = "hama";

sim = flows.Simulation(x_domain = xd, 
                        y_domain = yd,
                        z_domain = zd, 
                        y0 = y0, 
                        flow_type = flow_type, 
                        plot_type = plot_type, 
                        extra_args = extra_args);
  
# Run the simulation                                
sim.Run();
                                
                                