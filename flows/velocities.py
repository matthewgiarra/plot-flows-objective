import numpy as np
import pdb
import scipy
# from pycse import bvp

# This function parses a vector
# of dimensions [2n x 1] into 
# two vectors each of length [n x 1]
def Parse_Vector_2d(xy):
    num_points = len(xy);
    breakpoint = int(num_points / 2);
    x = xy[0 : breakpoint];
    y = xy[breakpoint : num_points];
    return x, y;

# This function parses a vector 
# of dimensions [3n x 1] into 
# two vectors each of length [n x 1]
'''
def Parse_Vector_3d(xyz):
	num_points = len(xyz)
	breakpoint = int(num_points / 3)
	x = xyz[0 : breakpoint]
	y = xyz[breakpoint : 2 * breakpoint]
	z = zyz[2 * breakpoint : 3 * breakpoint]
	
	return x, y, z;
'''
 
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

# UniformVelocity function
def UniformVelocity(xy, t, extra_args):
    
    # Parse the input
    x, y = Parse_Vector_2d(xy);
    
	# a represents u_o/U
    a = extra_args[0]
	
	# define U
    U = 1
	
    # Velocity
    u_o = np.full((1, x.shape[0]), a * U);
    v_o = np.zeros(y.shape[0]);
        
    vels = np.append(u_o, v_o);
    
    return vels;

# Poiseuille velocity function (Laminar)	
def PoiseuilleVelocity(xy, t, extra_args):

    # Parse the input
	x, y = Parse_Vector_2d(xy)
	
	# Velocity Constant
	a = extra_args[0]
	
	
	# Maximum Velocity
	umax = 1
	
	# Pipe radius
	R = 3
	
    # Mean horizontal velocity and vertical velocities
	u_vel = umax * a * (1 - (np.power(y, 2) / np.power(R, 2)))
	v_vel = np.zeros((y.shape[0]))
	
	#pdb.set_trace();
	vels = np.append(u_vel, v_vel)
	
	
	return vels

# Womersley velocity function
# To be revised to non-dimentional inputs
def Womersley(xy, t, extra_args):
	
	#Parse the input
	x, y = Parse_Vector_2d(xy)
	
	#Radius
	#R = extra_args[3]
	R = 1
	
	#Frequency
	Omega = extra_args[2]
	
	#Density
	#Rho = extra_args[1]
	Rho = 1
	
	#Dynamic Viscosity
	#Mu = extra_args[4]
	Mu = 1
	
	#Complex Number
	A = extra_args[1]
		
	#Kinematic Viscosity
	#Nu = Mu / Rho
	
	#pdb.set_trace();
	
	#Womersley Number
	Wm = extra_args[0]
	
	#Parts of Womersley Flow Solution
	v1 = y * np.power(1j,3/2) * Wm / R 
	v2 = np.power(1j,3/2) * Wm 
	
	# Mean horizontal velocity and vertical velocities
	u_vel = (- A * np.power(R, 2) / 1j / Mu / np.power(Wm, 2) * (1 - scipy.special.jv(0,v1)/ scipy.special.jv(0,v2)) * np.exp(1j * Omega * t)).real
	v_vel = np.zeros((y.shape[0]))
	
	#pdb.set_trace();
	
	vels = np.append(u_vel, v_vel)
	
	return vels
	
'''
# Blasius Boundary Layer
def Blasius(xy, t, extra_args):
    
    def odefun(F, x):
        f1, f2, f3 = F
        return [f2, f3, -0.5 * f1 * f3]
    def bcfun(Fa, Fb):
        return [Fa[0], Fa[1], 1.0 - Fb[1]]
    
    eta = np.linspace(0, 6, 100)
    f1init = eta
    f2init = np.exp(-eta)
    f3init = np.exp(-eta)
	
    Finit = np.vstack([f1init, f2init, f3init])
    sol = bvp(odefun, bcfun, eta, Finit)
	
	
	
    x, y = Parse_Vector_2d(xy)
	
    u_vel = U * 
    v_vel = 
	
    vels = np.append(u_vel, v_vel)
    return vels
'''


# OscillatingPlane
def OscillatingPlane(xy, t, extra_args):
    x, y = Parse_Vector_2d(xy)
	
	#Frequency
    Omega = extra_args[0]
	
	#Viscosity
    #Nu = extra_args[1]
    Nu = 1
	
	#Initial Velocity
    #u_o = extra_args[2]
    u_o = 1
	
	#Non-dimentional Convert
    Y = y / np.power((Nu/Omega),1/2)
    T = Omega * t
	
	#Exact Solution
    u_vel = u_o * np.exp(-Y/np.sqrt(2))* np.sin(T - Y / np.sqrt(2))
    v_vel = np.zeros((y.shape[0]))
	
	#Append to 1D
    vels = np.append(u_vel, v_vel)
	
    return vels








