import numpy as np
from .velocities import HamaVelocity, Parse_Vector_2d

# ODE integrator
from scipy.integrate import odeint

# Velocity vector class
class Velocity:
    def __init__(self, U = 0.0, V = 0.0, W = 0.0):
        self.U = U;
        self.V = V;
        self.W = W;
 
class CartesianCoordinate:
    def __init__(self, X = 0.0, Y = 0.0, Z = 0.0):
        self.X = X;
        self.Y = Y;
        self.Z = Z;
    
# Position vector class (cartesian)       
class Position:
    def __init__(self, X = 0.0, Y = 0.0, Z = 0.0):
        
        # Current position
        self.Current = CartesianCoordinate(X, Y, Z);
        
        # Position history
        self.History = CartesianCoordinate([X,], [Y,], [Z,]);
        
        # Position origin
        # I made this a single tuple
        # rather than assigning separate
        # fields for X, Y, Z mainly 
        # to make the origin immutable.
        self.Origin = (X, Y, Z);  
        
# This class represents a list of particles
# and methods that pertain to the entire field.
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
    
    # This function advects the entire field of particles
    # according to some velocity function. 
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
            x.append(self.Particles[k].Position.Current.X);
            y.append(self.Particles[k].Position.Current.Y);
            z.append(self.Particles[k].Position.Current.Z);
        
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
                self.Particles[k].Position.Current.X = x[k];
            if y is not None:
                self.Particles[k].Position.Current.Y = y[k];
            if z is not None:
                self.Particles[k].Position.Current.Z = z[k];
                  
    # This function increases
    # the ages of all the particles
    # by an amount
    def Age(self, amount = 1):
        
        # Number of particles
        num_particles = self.Count;
        
        # Loop over all the particles
        for k in range(num_particles):
            self.Particles[k].Age(amount);
                      
 # This class represents an individal particle.           
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
        
        # Set the age
        self.Duration = duration;
    
    # Update particle position based on some velocity field    
    def Advect(self, flow_type = None, t0 = 0.0, dt = 0.0, extra_args = None):
        if flow_type is not None:
            x0 = self.Position.Current.X;
            y0 = self.Position.Current.Y;
            z0 = self.Position.Current.Z;
            
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
                self.Position.Current.X = xy[0];
                self.Position.Current.Y = xy[1];
                
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
        
        
        
        