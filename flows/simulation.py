# Numpy for array operations, etc.
# import numpy as np

# Tools for python iterables
import itertools

# Python debugger
import pdb

# Particle stuff
from .particles import *

# Time
import time

# plotting
import matplotlib.pyplot as plt
from .plots import *
from mpl_toolkits.mplot3d import Axes3D

# This class represents parameters for a simulation
class Parameters:
    def __init__(self, plot_type = None,
                flow_type = None,
                extra_args = None,
                periodic_domain = True, NewParticleDistance = 0.5):
        self.PlotType = plot_type;
        self.FlowType = flow_type;
        self.ExtraArgs = extra_args;
        self.PeriodicDomain = True;
        self.RegenerateParticles = True;
        self.NewParticleDistance = NewParticleDistance;
        
        # Set options for the different types of plots
        #
        # Streakline plot
        if "streak" in plot_type.lower():
            self.PeriodicDomain = False;
            self.RegenerateParticles = False;
        
        # Pathline plot
        if "path" in plot_type.lower():
            self.PeriodicDomain = False;
            self.RegenerateParticles = False;
   
        # Vector field plot
        if "vector" in plot_type.lower():
            self.PeriodicDomain = False;
            self.RegenerateParticles = False;			
                      
# Extents of a domain
class Domain:
    def __init__(self, x = (0, 0), y = (0, 0), z = (0, 0)):
        self.X = x;
        self.Y = y;
        self.Z = z;

# Timing. Probably need to make this iterable.
class Time:
    def __init__(self, start = 0, stop = np.inf, step = 0.1):
        self.Start = start;
        self.Stop = stop;
        self.Step = step;
        self.Current = start;

# A flow simulation.
class Simulation:
    def __init__(self, x_domain = (-1,1), y_domain = (-1, 1), z_domain = (-1, 1), 
                                    y0 = [[0,], [0,], [0,]], time = [0, np.inf, 0.1], 
                                    flow_type = None, plot_type = "streak", 
                                    NewParticleDistance = 0.5, extra_args = None):
            
            # Define the extents of the domain
            domain = Domain(x_domain, y_domain, z_domain);
            
            # Assign the domain to the simulation
            self.Domain = domain;
            
            # Some simulation parameters
            self.Parameters = Parameters(plot_type = plot_type,
                                        flow_type = flow_type,
                                        extra_args = extra_args, 
                                        NewParticleDistance = NewParticleDistance);

            # Parse the particle positions
            x = y0[0];
            y = y0[1];
            z = y0[2];
            
            #pdb.set_trace();
			
            if "vector" in plot_type.lower():
                field = VelocityField(x_domain, y_domain, z_domain);
				
                self.VelocityField = field;
				
            elif "streak" in plot_type.lower() or "path" in plot_type.lower():
			    # Create a particle field
                field = ParticleField(x, y, z);
            
                # Assign the field to the simulation;
                self.ParticleField = field;
            
            # Assign the initial time and the time step
            start_time = time[0];
            stop_time  = time[1];
            step_time  = time[2];
            
            # Make the time object
            time = Time(start_time, stop_time, step_time);
            
			# Assign the time object to the simulation
            self.Time = time;
            
            # Counter to keep track of the iterations
            self.IterationNumber = itertools.count(0);
                        
    # Run the simulation
    def Run(self, make_plots = True):
        """

        :rtype: object
        """
        current_time = self.Time.Current;
        stop_time = self.Time.Stop;
        
        # Get domain
        xd = list(self.Domain.X);
        yd = list(self.Domain.Y);
        
        # if make_plots:
        #fig = plt.figure();
        #ax = fig.add_subplot(projection = '3d');
        #line1, = ax.plot([], [], [], '-k');
        #fig.show();
       
        
        # Read the plot type
        plot_type = self.Parameters.PlotType.lower();
        
        # and (self.ParticleField.Count) > 0
        # Run while the current time is
        # less than the stopping time
        try:        
            while current_time < stop_time:
                
                # Proceed with the next time step
                self.Step();
               
                # Clear axes
                ax.clear();
                plt.hold(True)
                    
                if "vector" in plot_type:
                    VectorField(ax, VelocityField = self.VelocityField); 
                # Check streaklines plot
                     
                elif "streak" in plot_type:
                    StreakPlot(ax, ParticleField = self.ParticleField);
                    #x_streak, y_streak, z_streak, d_streak = self.ParticleField.GetStreaklines();
                    
                elif "path" in plot_type:
                    PathlinePlot(ax, ParticleField = self.ParticleField);
                    
                elif "time" in plot_type:
                    TimelinePlot(ax, ParticleField = self.ParticleField);
                #pdb.set_trace() 
                ax.set_xlim(xd);
                ax.set_ylim([-2, 2]);
                fig.canvas.draw();
                time.sleep(0.05)
                current_time = self.Time.Current;
                plt.pause(0.0001)
                #pdb.set_trace() 
        except KeyboardInterrupt:
            pass
        
       
    # This iterates the simulation        
    def Step(self):
        
        plot_type = self.Parameters.PlotType.lower();
		
		
        # Advance the iteration number
        counter = next(self.IterationNumber);
        
		# Age the particles by one step
        if "vector" not in plot_type:
            self.ParticleField.Age();
        
        # Print
        #print("On iteration " + str(counter));
        
        # Starting time of the iteration
        t0 = self.Time.Current;
        dt = self.Time.Step;
        #print(t0)
        #print(dt)
        
		#Choose between plot type
        if "vector" in plot_type:
		    # Advect the field
            self.VelocityField.Advect(flow_type = self.Parameters.FlowType,
                t0 = t0, dt = dt, extra_args = self.Parameters.ExtraArgs);
				
            self.Time.Current += dt;
        else: 
			# Advect the particles
            self.ParticleField.Advect(flow_type = self.Parameters.FlowType,
                t0 = t0, dt = dt, extra_args = self.Parameters.ExtraArgs);
			
			# Remove the dead particles
			# and add some more back in
            self.CircleOfLife();
			
			# Count the number of particles
            num_alive = self.ParticleField.Count;
			
			# Print the number alive
			#print("Live particles: " + str(num_alive));
			
			# Increment the time counter
            self.Time.Current += dt;
			
			# Print a blank line
			#print();

    # This method creates and destroys particles
    # between time steps.
    def CircleOfLife(self): 
        
        # Read some parameters
        # 
        # Boolean: Regenerate dead particles?
        regenerate_particles = self.Parameters.RegenerateParticles;
        
        # Plot type
        plot_type = self.Parameters.PlotType;
        
        # Kill the particles that are outside the
        # domain and return their indices
        self.KillEscapees();
        
        # Regenerate the dead particles
        # if that option is specified.
        if regenerate_particles:
            self.RegenerateParticles();
                                
        # Put new particles at the streakline injection
        # if that option was specified.
        if "streak" in plot_type.lower():
            self.InjectStreaklines();
        
        # Remove the dead particles
        self.ParticleField.RemoveDead();            

    # This function checks if particles have left a domain
    # and returns their index within the particle field
    # structure, as well as flags specifying the direction
    # in which the particles left the domain.
    #
    # The function can also kill the particles
    # that leave the domain.
    def KillEscapees(self):
        
        # Loop over all of the particles
        # and mark as dead any that have
        # left the domain.
        for k in range(self.ParticleField.Count):
            
            # Check if the particle has left the domain for each axis
            # The function called here returns true 
            # for the axes of the domain
            # where the particle has not left.
            x_in, y_in, z_in = self.ParticleField.Particles[k].CheckInOut(self.Domain);
            
            # Mark the particle dead
            # if it left the domain.
            if not(x_in and y_in and z_in):
                self.ParticleField.Particles[k].Kill();
                
    # This method regenerates dead particles    
    def RegenerateParticles(self):
        
        # Read some options
        periodic_domain = self.Parameters.PeriodicDomain;
                
        # Allocate new coordinates
        xnew = [];
        ynew = [];
        znew = [];
        
        # Read the domain
        x_domain = self.Domain.X;
        y_domain = self.Domain.Y;
        z_domain = self.Domain.Z;
        
        # Loop over the particles that have left
        for k in range(self.ParticleField.Count):
            
            # Replace the dead particles.
            if self.ParticleField.Particles[k].Alive is False:
                
                # Read particle positions
                x = self.ParticleField.Particles[k].Position.Current.X;
                y = self.ParticleField.Particles[k].Position.Current.Y;
                z = self.ParticleField.Particles[k].Position.Current.Z;
                
                # Check if particle is outside the domain
                x_in, y_in, z_in = self.ParticleField.Particles[k].CheckInOut(self.Domain);
                
                # Periodic domain regeneration
                if periodic_domain is True:
                    # Figure out which edge of the domain 
                    # the particle exited through, and 
                    # set the new particle starting position
                    # opposite that face. 
                    #
                    # Check X faces
                    if not(x_in):
                        if x < x_domain[0]:
                            x = x_domain[1];
                        else:
                            x = x_domain[0];
            
                    # Check Y faces
                    if not(y_in):
                        if y < y_domain[0]:
                            y = y_domain[1];
                        else:
                            y = y_domain[0];
                    
                    # Check Z faces
                    if not(y_in):
                        if z < z_domain[0]:
                            z = z_domain[1];
                        else:
                            z = z_domain[0];
                else:
                    # Read the particle origin
                    particle_origin = self.ParticleField.Particles[k].Position.Origin;
                    
                    # Parse the coordinates
                    x = particle_origin[0];
                    y = particle_origin[1];
                    z = particle_origin[2];
                    
                # Append the new positions to the list      
                x_new.append(x);
                y_new.append(y);
                z_new.append(z);
                
                # Make the new particles
                self.ParticleField.CreateParticles(x_new, y_new, z_new);    
            
    # This function injects streakline particles.
    def InjectStreaklines(self):
        
        # Number of streakline starting points.
        num_streaklines = len(self.ParticleField.InitialPositions.X);
        
        # Get the coordinates of all of the
        # streakline starting points.
        x_init = self.ParticleField.InitialPositions.X;
        y_init = self.ParticleField.InitialPositions.Y;
        z_init = self.ParticleField.InitialPositions.Z;
        
        # Loop over the injection points
        for n in range(num_streaklines):
            
            # Read the startign position of
            # the n'th streak injection point
            y0 = (x_init[n], y_init[n], z_init[n]);
            
            # Extract the streakline
            xs, ys, zs, ds = self.ParticleField.GetStreakline(y0);
        
            # Components of the vector from the origin
            dx = xs[0] - x_init[n];
            dy = ys[0] - y_init[n];
            dz = zs[0] - z_init[n];
            
            # Distance of the particle from the origin
            dr = np.sqrt(dx**2 + dy**2 + dz**2);
            
            # Create a new particle if the previous particle
            # has moved far enough away from its origin.
            if dr > self.Parameters.NewParticleDistance:
                self.ParticleField.CreateParticles([x_init[n],], [y_init[n],], [z_init[n],]);
    
