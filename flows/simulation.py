
# Numpy for array operations, etc.
# import numpy as np

# Plotting
import matplotlib.pyplot as plt

# Tools for python iterables
import itertools

# Python debugger
import pdb

# Particle stuff
from .particles import *

# Velocity functions
from .velocities import *

# Time
import time

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
        if "streak" in plot_type.lower():
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

            # pdb.set_trace();
            
            # Initial positions of the particles
            self.InitialPositions = CartesianCoordinate(x, y, z);

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
        current_time = self.Time.Current;
        stop_time = self.Time.Stop;
        
        # Get domain
        xd = list(self.Domain.X);
        yd = list(self.Domain.Y);
        
        # if make_plots:
        fig = plt.figure();
        ax = fig.add_subplot(111);
        line1, = ax.plot([], [], '-k');
        fig.show();
        
        # Run while the current time is
        # less than the stopping time
        try:
            while current_time < stop_time and (self.ParticleField.Count) > 0:
                
                # Proceed with the next time step
                self.Step();
                
                # Get a streakline
                x_streak, y_streak, z_streak, d_streak = self.GetStreaklines();
                
                # Number of streaklines
                n_streaks = len(x_streak);
            
                # Make the plots
                if make_plots:
                    # xplot, yplot, zplot = self.ParticleField.GetCoordinates();
                    # line1.set_xdata(xplot);
                    # line1.set_ydata(yplot);
                    # line1.set_xdata(x_streak[0]);
                    # line1.set_ydata(y_streak[0]);
                    ax.clear();
                    for n in range(n_streaks):

                        # Number of points in this streakline
                        npoints = len(x_streak[n]);

                        # Plot this streak line
                        ax.plot(x_streak[n], y_streak[n], '.');
                    
                    # pdb.set_trace();
                    ax.set_xlim(xd);
                    ax.set_ylim([-2, 2]);
                    fig.canvas.draw();
                    time.sleep(0.01)
        except KeyboardInterrupt:
            pass
    
    # This iterates the simulation        
    def Step(self):
        
        # Advance the iteration number
        counter = next(self.IterationNumber);
        
        # Age the particles by one step
        self.ParticleField.Age();
        
        # Print
        print("On iteration " + str(counter));
        
        # Starting time of the iteration
        t0 = self.Time.Current;
        dt = self.Time.Step;
        
        # Advect the particles
        self.ParticleField.Advect(flow_type = self.Parameters.FlowType,
        t0 = t0, dt = dt, extra_args = self.Parameters.ExtraArgs);
        
        # Remove the dead particles
        # and add some more back in
        self.CircleOfLife();
        
        # Count the number of particles
        num_alive = self.ParticleField.Count;
        
        # Print the number alive
        print("Number alive: " + str(num_alive));
        
        # Increment the time counter
        self.Time.Current += dt;

    # This method creates and destroys particles
    # between time steps.
    def CircleOfLife(self):
        
        # Query the simulation iteration number
        counter = self.IterationNumber;      
        
        # Read some options
        periodic_domain = self.Parameters.PeriodicDomain;
        regenerate_particles = self.Parameters.RegenerateParticles;
        
        # Flow type
        flow_type = self.Parameters.FlowType;
        
        # Plot type
        plot_type = self.Parameters.PlotType;
        
        # Read in the domain
        x_domain = self.Domain.X;
        y_domain = self.Domain.Y;
        z_domain = self.Domain.Z;
        
        # Count the number of particles
        num_particles = self.ParticleField.Count;
        
        # Allocate space for new particle positions
        x_new = [];
        y_new = [];
        z_new = [];
        
        # Loop over the particles
        for k in range(num_particles):
            
            # Read particle positions
            x = self.ParticleField.Particles[k].Position.Current.X;
            y = self.ParticleField.Particles[k].Position.Current.Y;
            z = self.ParticleField.Particles[k].Position.Current.Z;
            
            # Determine whether the particle left the
            # domain on each face.
            x_out = not(x_domain[0] <= x <= x_domain[1]);
            y_out = not(y_domain[0] <= y <= y_domain[1]);
            z_out = not(z_domain[0] <= z <= z_domain[1]);
            
            # Kill the particle if it's outside the range
            if x_out or y_out or z_out:
                
                # This flags the particle as dead,
                # but doesn't remove it from the
                # particle field yet. That part
                # is handled by ParticleField.RemoveDead()
                self.ParticleField.Particles[k].Kill();
                
                # Check if particles should be regenerated.
                if regenerate_particles is True:
                
                    # This is for periodic domains:
                    # the new particle will be generated
                    # on the face opposite to where
                    # it exited the domain
                
                    if periodic_domain is True:
                        # Figure out which edge of the domain 
                        # the particle exited through, and 
                        # set the new particle starting position
                        # opposite that face. 
                        #
                        # Check X faces
                        if x_out:
                            if x < x_domain[0]:
                                x = x_domain[1];
                            else:
                                x = x_domain[0];
                
                        # Check Y faces
                        if y_out:
                            if y < y_domain[0]:
                                y = y_domain[1];
                            else:
                                y = y_domain[0];
                        
                        # Check Z faces
                        if z_out:
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
                    
                    # Make new particles
                    self.ParticleField.CreateParticles(x_new, y_new, z_new);       
        
        # Put new particles at the streakline injection
        # if it's a streak plot
        if "streak" in plot_type.lower():
            print("Make new particle!")
            
            # Number of streakline starting points.
            num_pos = len(self.InitialPositions.X);
            
            # Get the coordinates of all of the
            # streakline starting points.
            x_init = self.InitialPositions.X;
            y_init = self.InitialPositions.Y;
            z_init = self.InitialPositions.Z;
            
            # Loop over the injection points
            for n in range(num_pos):
                
                # Read the startign position of
                # the n'th streak injection point
                y0 = (x_init[n], y_init[n], z_init[n]);
                
                # Extract the streakline
                xs, ys, zs, ds = self.GetStreakline(y0);
            
            
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
            
            
            # Make the new particles
            # self.ParticleField.CreateParticles(x_new, y_new, z_new);
        
        # Remove the dead particles
        self.ParticleField.RemoveDead();
             

    # This function gets all of the streaklines
    def GetStreaklines(self):
        # Number of starting particles
        num_starting_points = len(self.InitialPositions.X);
        
        # Initialize vectors
        x = [];
        y = [];
        z = [];
        durations = [];
        
        # Loop over all the starting points
        for k in range(num_starting_points):
            
            # Get positions
            x0 = self.InitialPositions.X[k];
            y0 = self.InitialPositions.Y[k];
            z0 = self.InitialPositions.Z[k];
            
            # Get a single streakline
            x_streak, y_streak, z_streak, d_streak = self.GetStreakline(y0 = (x0, y0, z0));
            
            # Append the streaklines
            x.append(x_streak);
            y.append(y_streak);
            z.append(z_streak);
            durations.append(d_streak);
            
        # Return the list
        return x, y, z, durations
            
    
    # This function gets the coordinates of all of the
    # particles that started at a certain point, 
    # i.e., a streakline.  
    def GetStreakline(self, y0 = (0, 0, 0), smooth = False, num_smooth_points = 100):
        # Read the number of particles
        
        # Number of particles
        num_particles = self.ParticleField.Count;
        
        # Create vectors for the coordinates
        x = [];
        y = [];
        z = [];
        durations = [];

        # Loop over all the particles and 
        # compare the origin of each
        # with the queried starting position
        for k in range(num_particles):
            particle = self.ParticleField.Particles[k];
            
            # Particle origin
            particle_origin = particle.Position.Origin;
            
            # Check if the input location corresponds 
            # to the particle origin
            if y0 == particle_origin:
                x_new = particle.Position.Current.X;
                y_new = particle.Position.Current.Y;
                z_new = particle.Position.Current.Z;
                durations_new = particle.Duration;
                
                # Append to list
                x.append(x_new);
                y.append(y_new);
                z.append(z_new);
                durations.append(durations_new);
                
        # Now sort them by age
        xs = [k for (durations, k) in sorted(zip(durations, x))];
        ys = [k for (durations, k) in sorted(zip(durations, y))];
        zs = [k for (durations, k) in sorted(zip(durations, z))];
        ds = sorted(durations);
        
        # Smoothing
        if smooth is True:
            
            if len(xs) > 3:
                # pdb.set_trace();
                tck_x = interpolate.splrep(ds, xs, s = 0);
                tck_y = interpolate.splrep(ds, ys, s = 0);
                tck_z = interpolate.splrep(ds, zs, s = 0);
                dur_new = np.linspace(min(ds), max(ds), num_smooth_points);
            
                xs = interpolate.splev(dur_new, tck_x, der = 0);
                ys = interpolate.splev(dur_new, tck_y, der = 0);
                zs = interpolate.splev(dur_new, tck_z, der = 0)
                ds = dur_new;
        
        # Return the coordinates and durations (the ages)
        # of the points on the streakline, sorted by duration.
        return xs, ys, zs, ds
        