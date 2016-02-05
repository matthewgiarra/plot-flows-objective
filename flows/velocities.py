import numpy as np

# This function parses a vector
# of dimensions [2n x 1] into 
# two vectors each of length [n x 1]
def Parse_Vector_2d(xy):
    num_points = len(xy);
    breakpoint = int(num_points / 2);
    x = xy[0 : breakpoint];
    y = xy[breakpoint : num_points];
    return x, y;

# Poiseuille velocity function
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
 
# Hama velocity function   
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