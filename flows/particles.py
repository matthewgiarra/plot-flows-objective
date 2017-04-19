import numpy as np
from .velocities import *
from scipy.integrate import odeint
import pdb


# Velocity vector class
class Velocity:
    def __init__(self, U = 0.0, V = 0.0, W = 0.0):
        self.U = U;
        self.V = V;
        self.W = W;
 
class CartesianCoordinates:
    def __init__(self, X = 0.0, Y = 0.0, Z = 0.0):
        self.X = X;
        self.Y = Y;
        self.Z = Z;
  
# Position vector class (cartesian)       
class Position:
    def __init__(self, X = 0.0, Y = 0.0, Z = 0.0):
        
        # Current position
        self.Current = CartesianCoordinates(X, Y, Z);
        
        # Position history
        self.History = CartesianCoordinates([X,], [Y,], [Z,]);
        
        # Position origin
        # I made this a single tuple
        # rather than assigning separate
        # fields for X, Y, Z mainly 
        # to make the origin immutable.
        self.Origin = (X, Y, Z);  

# This class represents the entire field with the velocities vector
class VelocityField:
    def __init__(self, x = [0,], y = [0,], z = [0,]):
		
		
		# Initialize field
        xv = np.linspace(x[0], x[1], num = 20, endpoint = True)
        yv = np.linspace(y[0], y[1], num = 20, endpoint = True)
        zv = np.linspace(z[0], z[1], num = 20, endpoint = True)
		
        #pdb.set_trace()
		
        
	    
		# Meshgrid field
        # xm, ym, zm = np.meshgrid(xv, yv, zv)
		
        # Save coordinates
        self.Coordinates = CartesianCoordinates(xv, yv, zv);
		
        #pdb.set_trace()
		# Initialize Velocity
        vel = Velocity([], [], []);
        self.Velocity = vel;

    def Advect(self, flow_type = None, t0 = 0, dt = 0, extra_args = None):
	    
        if flow_type is not None:
            # Retrieve coordinates
            x0 = self.Coordinates.X;
            y0 = self.Coordinates.Y;
            z0 = self.Coordinates.Z;
            xy0 = np.array(x0.tolist() + y0.tolist());
			
			# Get shape to resize array 
            #i = x0.shape[0];
            #j = y0.shape[0];
            #k = z0.shape[0];
			
            #pdb.set_trace();
			
            # Time vector
            tf = t0 + dt;
            
			# Choose between fields 
            if "hama" in flow_type.lower():
                uv = HamaVelocity(xy0, tf, extra_args)
            elif "poiseuille" in flow_type.lower():               
                uv = PoiseuilleVelocity(xy0, tf, extra_args)
            elif "womersley" in flow_type.lower():               
                uv = Womersley(xy0, tf, extra_args)
            elif "blasius" in flow_type.lower():              
                uv = Blasius(xy0, tf, extra_args)
            elif "oscillatingplane" in flow_type.lower():      
                uv = OscillatingPlane(xy0, tf, extra_args)
            elif "uniform" in flow_type.lower():
                uv = UniformVelocity(xy0, tf, extra_args)
				
     	    # Extract the new velocities
            u, v = Parse_Vector_2d(uv);
            #pdb.set_trace()  
			# Convert 1D to 3D
            #u_new = u.reshape((i, j, k))
            #v_new = v.reshape((i, j, k))
            #pdb.set_trace()  
			
            
			# Set new velocities (2D for now)
            self.SetVelocity(u, v);
			
			
	# This function sets the velocities of the vectorfield
    def SetVelocity(self, u_new = 0, v_new = 0):
		      
      	# Extract the new velocities
        self.Velocity.U = u_new;
        self.Velocity.V = v_new;
        self.Velocity.W = np.zeros(u_new.shape);
        #pdb.set_trace() 
	
	# This function gets the velocities of the vectorfield
    def GetVelocity(self):
        # Initialize vectors
        u = [];
        v = [];
        w = [];
        x = [];
        y = [];
        z = [];
        
				
		# Get positions within the domain
        x = self.Coordinates.X
        y = self.Coordinates.Y
        z = self.Coordinates.Z
			
		# Get velocity of one point
        u = self.Velocity.U
        v = self.Velocity.V
        w = self.Velocity.W
                       
        #pdb.set_trace()
        return x, y, z, u, v, w;
    
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
        
        # Initial positions of the particles
        self.InitialPositions = CartesianCoordinates(x, y, z);
        
		
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
            #print(xy0)
            #pdb.set_trace()
    
            # Choose between fields 
			# Time vector has incorrect shape
            if "hama" in flow_type.lower():
                xy = (odeint(HamaVelocity, y0 = xy0, t = t, args = (extra_args,)))[-1, :]
                #uv = HamaVelocity(xy0, tf, extra_args)
            elif "poiseuille" in flow_type.lower():
                xy = (odeint(PoiseuilleVelocity, y0 = xy0, t = t, args = (extra_args,)))[-1, :]
                #uv = PoiseuilleVelocity(xy0, tf, extra_args)
            elif "womersley" in flow_type.lower():
                xy = (odeint(Womersley, y0 = xy0, t = t, args = (extra_args,)))[-1, :]
                #uv = Womersley(xy0, tf, extra_args)
            elif "blasius" in flow_type.lower():
                xy = (odeint(Blasius, y0 = xy0, t = t, args = (extra_args,)))[-1, :]
                #uv = Blasius(xy0, tf, extra_args)
            elif "oscillatingplane" in flow_type.lower():
                xy = (odeint(OscillatingPlane, y0 = xy0, t = t, args = (extra_args,)))[-1, :]
                #uv = OscillatingPlane(xy0, tf, extra_args)
            elif "uniform" in flow_type.lower():
                xy = (odeint(UniformVelocity, y0 = xy0, t = t, args = (extra_args,)))[-1, :]
                #uv = UniformVelocity(xy0, tf, extra_args)
			
			
			# Extract the new positions
            x_new, y_new = Parse_Vector_2d(xy);
            
			# Extract the new velocities
            #u_new, v_new = Parse_Vector_2d(uv);
			
			# Set new velocities
            #self.SetVelocity(flow_type, tf, extra_args);
	        
            # Set the new coordinates
            self.SetCoordinates(x = x_new, y = y_new);
             
            #pdb.set_trace()

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
    def GetStreakline(self, y0 = (0, 0, 0)):
        # Read the number of particles

        # Number of particles
        num_particles = self.Count;

        # Create vectors for the coordinates
        x = [];
        y = [];
        z = [];
        durations = [];

        # Loop over all the particles and 
        # compare the origin of each
        # with the queried starting position
        for k in range(num_particles):
            particle = self.Particles[k];

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

        # Return the coordinates and durations (the ages)
        # of the points on the streakline, sorted by duration.
        return xs, ys, zs, ds
            
                
                
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
                self.Particles[k].Position.History.X.append(x[k]);
            if y is not None:
                self.Particles[k].Position.Current.Y = y[k];
                self.Particles[k].Position.History.Y.append(y[k]);
            if z is not None:
                self.Particles[k].Position.Current.Z = z[k];
                self.Particles[k].Position.History.Z.append(y[k]);
         
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
                uv = HamaVelocity(xy0, t[1], extra_args)

                # Extract the new velocities
                self.Velocity.U = uv[0];
                self.Velocity.V = uv[1];
				
            if "poiseuille" in flow_type.lower():
                
                # New position
                try:
                    xy = (odeint(PoiseuilleVelocity, y0 = xy0, t = t, args = (extra_args,)))[1, :];
                except:
                    print("Error integrating!");
                    
                # Extract the new positions
                self.Position.Current.X = xy[0];
                self.Position.Current.Y = xy[1];
                
                # New velocity
                uv = PoiseuilleVelocity(xy0, t[1], extra_args)

                # Extract the new velocities
                self.Velocity.U = uv[0];
                self.Velocity.V = uv[1];
            
            if "womersley" in flow_type.lower():
                
                # New position
                try:
                    xy = (odeint(Womersley, y0 = xy0, t = t, args = (extra_args,)))[1, :];
                except:
                    print("Error integrating!");
                    
                # Extract the new positions
                self.Position.Current.X = xy[0];
                self.Position.Current.Y = xy[1];
                
                # New velocity
                uv = Womersley(xy0, t[1], extra_args)

                # Extract the new velocities
                self.Velocity.U = uv[0];
                self.Velocity.V = uv[1];
            
            if "blasius" in flow_type.lower():
			
			    # New position
                try:
                    xy = (odeint(Blasius, y0 = xy0, t = t[1], args = (extra_args,)))[1, :];
                except:
                    print("Error integrating!");
                    
                # Extract the new positions
                self.Position.Current.X = xy[0];
                self.Position.Current.Y = xy[1];
                
                # New velocity
                uv = Blasius(xy0, t, extra_args)

                # Extract the new velocities
                self.Velocity.U = uv[0];
                self.Velocity.V = uv[1];
            if "oscillatingplane" in flow_type.lower():
                
                # New position
                try:
                    xy = (odeint(OscillatingPlane, y0 = xy0, t = t, args = (extra_args,)))[1, :];
                except:
                    print("Error integrating!");
                    
                # Extract the new positions
                self.Position.Current.X = xy[0];
                self.Position.Current.Y = xy[1];
                
                # New velocity
                uv = OscillatingPlane(xy0, t[1], extra_args)

                # Extract the new velocities
                self.Velocity.U = uv[0];
                self.Velocity.V = uv[1];

            if "uniform" in flow_type.lower():
                
                # New position
                try:
                    xy = (odeint(UniformVelocity, y0 = xy0, t = t, args = (extra_args,)))[1, :];
                except:
                    print("Error integrating!");
                    
                # Extract the new positions
                self.Position.Current.X = xy[0];
                self.Position.Current.Y = xy[1];
                
                # New velocity
                uv = UniformVelocity(xy0, t[1], extra_args)

                # Extract the new velocities
                self.Velocity.U = uv[0];
                self.Velocity.V = uv[1];

				
    # This method checks whether a particle is
    # within a domain, and returns a boolean
    # for each direction of the domain
    # that is true if the particle is within
    # the domain, and false if the particle
    # is outside of the domain in that direction
    def CheckInOut(self, Domain = None):
        if Domain is not None:
            
            # Read in the domain
            x_domain = Domain.X;
            y_domain = Domain.Y;
            z_domain = Domain.Z;
    
            # Get the current particle position
            x = self.Position.Current.X;
            y = self.Position.Current.Y;
            z = self.Position.Current.Z;
            
            # Check if the particle is within the domain
            # These statements evaluate True
            # if the particle position is within 
            # the domain, and False if it is not.
            x_in = x_domain[0] < x < x_domain[1];
            y_in = y_domain[0] < y < y_domain[1];
            z_in = z_domain[0] < z < z_domain[1];
            
            # Return the Boolean results
            return x_in, y_in, z_in;
    
    # Kill particles
    def Kill(self):
        self.Alive = False;       
        
    # This function increases
    # the age of a single particle
    # by an amount
    def Age(self, amount = 1):
        # Age the particle
        self.Duration += amount;
        
        
        
