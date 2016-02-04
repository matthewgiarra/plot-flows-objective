import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from scipy import interpolate
import pdb
import itertools
import time

"""
flows.py

Classes and methods for plotting vectors, pathlines,
streamlines, and timelines of time-varying flows.

Object hierarchy:
Flow
    Flow.Velocity.U
    Flow.Velocity.V
    Flow.Velocity.W
    
    Flow.Domain.X
    Flow.Domain.Y
    Flow.Domain.Z
    Flow.Domain.T
    
    # Methods:
    # Plot vectors:
    Flow.PlotVelocityField(self, Flow_Type = None, time = None);
    
    # Plot stream lines
    Flow.PlotStreamLines(self, Flow_Type = None, time = None);
    
    # Plot path lines
    Flow.PlotPathLines(self, Flow_Type = None, Start_Time = 0, End_Time = None, X0 = None);
    
    # Plot streak lines
    Flow.PlotStreakLines(self, Flow_Type = None, Start_Time = 0, End_Time = None, X0 = None);
    
"""

# Define pi
pi = np.pi;

# Counter
particle_counter = itertools.count();
iter_counter = itertools.count();

class Velocity:
    def __init__(self, U = 0.0, V = 0.0, W = 0.0):
        self.U = U;
        self.V = V;
        self.W = W;
        
        
class Position:
    def __init__(self, X = 0.0, Y = 0.0, Z = 0.0):
        
        # Updateable particle coordinates
        self.X = X;
        self.Y = Y;
        self.Z = Z;
        
        # Particle origin.
        # This is where the particle was created.
        self.Origin = (X, Y, Z);

def Parse_Vector_2d(xy):
    num_points = len(xy);
    breakpoint = int(num_points / 2);
    x = xy[0 : breakpoint];
    y = xy[breakpoint : num_points];
    return x, y;
    
def Get_Pathline(xy, k = 0):
    num_points = int(xy.shape[1] / 2);
    x = xy[:, k];
    y = xy[:, k + num_points];
    
    return x, y

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

# Parameters for the simulation
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

            # Initial positions of the particles
            self.InitialPositions = Position(x, y, z);

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
            
            # Keep track of the iterations
            self.IterationNumber = 0;
            
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
        counter = next(iter_counter);
        self.IterationNumber = counter;
        
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
            x = self.ParticleField.Particles[k].Position.X;
            y = self.ParticleField.Particles[k].Position.Y;
            z = self.ParticleField.Particles[k].Position.Z;
            
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
            if particle.Position.Origin == y0:
                x_new = particle.Position.X;
                y_new = particle.Position.Y;
                z_new = particle.Position.Z;
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
        

class ParticleField:
    def __init__(self, x = [0,], y = [0,], z = [0,], durations = None):
        
        # Number of particles to create
        num_particles = np.min([len(x), len(y), len(z)]);
        
        # Default ages
        if durations is None:
            durations = list(np.zeros(num_particles));
        
        # Create the particles     
        Particles = [Particle(x[count], y[count], z[count], durations[count]) 
                        for count in range(num_particles)];
                        
        # Assign the field
        self.Particles = Particles;
        
        # Count the number of particles
        self.Count = num_particles;
        
        # Count the alive particles
        self.Alive = num_particles;
        
        # Count dead particles
        self.Dead = 0;
     
    def Advect(self, flow_type = None, t0 = 0, dt = 0, extra_args = None):

        if flow_type is not None:
    
            # Number of particles
            num_particles = len(self.Particles);
        
            # Get the positions of all of the particles
            x0, y0, z0 = self.GetCoordinates();
        
            # Time vector
            tf = t0 + dt;
        
            # Time vector
            t = np.array([t0, tf]);
        
            # Initial positions as a numpy vector
            # in the form compatible with the 
            # velocity functions
            xy0 = np.array(x0 + y0);
            # xy0 = np.array([x0, y0]);
    
            # Choose between fields
            if "hama" in flow_type.lower():
                xy = (odeint(HamaVelocity, y0 = xy0, t = t, args = (extra_args,)))[-1, :];
            
            # Extract the new positions
            x_new, y_new = Parse_Vector_2d(xy);
            
            # Set the new coordinates
            self.SetCoordinates(x = x_new, y = y_new);    
                
    # Remove dead particles
    def RemoveDead(self):
        # Update the number of particles to start
        self.Count = len(self.Particles);
        
        # Loop over particles    
        for k in range(self.Count - 1, -1, -1):
            if self.Particles[k].Alive is False:
                del self.Particles[k];
                
        # Update the number of particles
        self.Count = len(self.Particles);
    
    # Make new particles
    def CreateParticles(self, x = [0, ], y = [0, ], z = [0, ], durations = None):
        
        # New particle field
        new_field = ParticleField(x, y, z, durations);
        
        # Extend the particle list
        self.Particles.extend(new_field.Particles);
        
        # Update the count
        self.Count = len(self.Particles);
        
    # Kill particles
    def KillParticles(self, inds = None):
        if inds is not None:
            num_particles = self.Count;
            
            # Set particles to dead
            for k in inds:
                self.Particles[k].Kill;     
        
            # Clean up the list
            self.RemoveDead();
    
    # This function reads the coordinates of all of the particles.        
    def GetCoordinates(self):
        
        # Count the particles
        num_particles = self.Count;
        
        # Initialize positions
        x = [];
        y = [];
        z = [];
        
        # Loop over all the particles
        for k in range(num_particles):
            x.append(self.Particles[k].Position.X);
            y.append(self.Particles[k].Position.Y);
            z.append(self.Particles[k].Position.Z);
        
        # Return the vectors    
        return x, y, z
    
    # This function sets the coordinates of all the particles in a field
    # as the values of some vectors x_new, y_new, z_new    
    def SetCoordinates(self, x = None, y = None, z = None):
        
        # Read number of particles
        num_particles = self.Count;
        
        # Loop over all the particles
        for k in range(num_particles):
            if x is not None:
                self.Particles[k].Position.X = x[k];
            if y is not None:
                self.Particles[k].Position.Y = y[k];
            if z is not None:
                self.Particles[k].Position.Z = z[k];
            
        
    # This function increases
    # the ages of all the particles
    # by an amount
    def Age(self, amount = 1):
        
        # Number of particles
        num_particles = self.Count;
        
        # Loop over all the particles
        for k in range(num_particles):
            self.Particles[k].Age(amount);
    
class Particle:
    def __init__(self, x = 0.0, y = 0.0, z = 0.0, u = 0.0, v = 0.0, w = 0.0, duration = 0):
        
        # Velocity
        vel = Velocity(u, v, w);
        self.Velocity = vel;
        
        # Position
        pos = Position(x, y, z);
        self.Position = pos;
        
        # Alive flag
        self.Alive = True;
        
        # Set the ID
        self.ID = next(particle_counter);
        
        # Set the age
        self.Duration = duration;
    
    # Update particle position based on some velocity field    
    def Advect(self, flow_type = None, t0 = 0.0, dt = 0.0, extra_args = None):
        if flow_type is not None:
            x0 = self.Position.X;
            y0 = self.Position.Y;
            z0 = self.Position.Z;
            
            # Final time
            tf = t0 + dt;
            
            # Time vector
            t = np.array([t0, tf]);
            
            # Initial positions as a numpy vector,
            # in the form compatible with the velocity functions.
            xy0 = np.array([x0, y0]);
            
            # Choose between velocity fields
            if "hama" in flow_type.lower():
                
                # New position
                try:
                    xy = (odeint(HamaVelocity, y0 = xy0, t = t, args = (extra_args,)))[1, :];
                except:
                    print("Error integrating!");
                    
                # Extract the new positions
                self.Position.X = xy[0];
                self.Position.Y = xy[1];
                
                # New velocity
                uv = HamaVelocity(xy0, t, extra_args)

                # Extract the new velocities
                self.Velocity.U = uv[0];
                self.Velocity.V = uv[1];
    
    # Kill particles
    def Kill(self):
        self.Alive = False;
        
        
    # This function increases
    # the age of a single particle
    # by an amount
    def Age(self, amount = 1):

        # Age the particle
        self.Duration += amount;
            

# Vectorize means make this work element-wise
# @np.vectorize  
def PoiseuilleVelocity(xy, t, args1):

    # Parameters
    amp = args1[0];
    freq = args1[1];
    
    # Coordinates
    x = xy[0];
    y = xy[1];
    
    # Calculate velocities
    u_vel = amp * (1 - y**2);
    v_vel = amp * np.sin(2 * np.pi * (freq * x  + t));
    
    # Return u and v
    return [u_vel.astype(float), v_vel.astype(float)];
    
def HamaVelocity(xy, t, extra_args):

    # Amplitude (0.2 seems to work)
    a = extra_args[0];
    
    # Wave number (unity in the paper)
    alpha = extra_args[1]
    
    # Wave velocity (unity in the paper)
    c = extra_args[2];
    
    # pdb.set_trace();
    
    # Parse the input
    x, y = Parse_Vector_2d(xy);
    
    # Mean horizontal and vertical velocities
    u_o = 1 + np.tanh(y);
    v_o = 0;
    
    # Fluctuating velocity (horizontal)
    u_prime = 2 * a * 1 / np.cosh(y) * np.tanh(y) * np.sin(alpha * (x - c * t));
    
    # Fluctuating velocoity (vertical)
    v_prime = 2 * a * 1 / np.cosh(y) * np.cos(alpha * (x - c * t));
    
    # Add steady and unsteady velocities
    u_vel = u_o + u_prime;
    v_vel = v_o + v_prime;
    
    # Append into 1D
    vels = np.append(u_vel, v_vel);
    
    # Concatonate into a list
    return vels;
    
def UniformVelocity(xy, t, extra_args):
    
    # Parse the input
    x, y = Parse_Vector_2d(xy);
    
    num_points = len(x);
    
    # Velocity
    u_o = extra_args * np.ones((num_points,), dtype = np.float);
    v_o = 0 * np.ones((num_points,), dtype = np.float);
    
    u_vel = u_o;
    v_vel = v_o;
    
    vels = np.append(u_vel, v_vel);
    
    return vels;
         
def plot_streaklines(flow_type = None, t = [0], y0 = [0], my_args = None):
    if flow_type is not None:
        
        # Size of the time step
        time_step_size = 0.1;
        
        # Number of flow tracer sources
        num_sources = int(len(y0) / 2);
        
        
        
        # List of X coordinates of all the particles
        # that will be tracked
        
        
        # Allocate all of the coordinate positions
        # so we don't have to re-allocate memory all the time
                    
            



    
def plot_pathlines(flow_type = None, t = [0], y0 = [0], my_args = None):
    if flow_type is not None:
        
        fig = plt.figure();
        
        if "hama" in flow_type.lower():
            # Extract the constants related to the Hama flow
                        
            # pdb.set_trace();
                        
            # Integrate
            xy = odeint(HamaVelocity, y0 = y0, t = t, args = (my_args,));
            ax = plt.axes(xlim=(-1, 9), ylim=(-2.0, 2.0))
            
            
        elif "uniform" in flow_type.lower():
            
            xy = odeint(UniformVelocity, y0 = y0, t = t, args = (my_args,));
            ax = plt.axes(xlim=(-1, 9), ylim=(-2, 2))
            
        # Count number of path lines
        num_pathlines = int(xy.shape[1] / 2);
        
        # Do the plotting
        for k in range(num_pathlines):
            print("Figure: " + str(k));
            x, y = Get_Pathline(xy, k);
            ax.plot(x, y, '-k');
            
        plt.show(ax)
        










