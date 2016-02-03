import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
import pdb


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

class Velocity:
    def __init__(self, U = None, V = None, W = None):
        self.U = U;
        self.V = V;
        self.W = W;

class Domain:
    def __init__(self, X = None, Y = None, Z = None):
        if X is not None and (type(X) is list or type(X) is tuple) and len(X) == 3:
                xv = np.linspace(X[0], X[1], X[2]);
        else:
            xv = 0;
            
        if Y is not None and (type(Y) is list or type(Y) is tuple) and len(Y) == 3:
            yv = np.linspace(Y[0], Y[1], Y[2]);
        else:
            yv = 0;
            
        if Z is not None and (type(Z) is list or type(Z) is tuple) and len(Z) == 3:
            zv = np.linspace(Z[0], Z[1], Z[2]);
        else:
            zv = 0;
        
       
        # Form the matrix of coordinates
        xm, ym, zm = np.meshgrid(xv, yv, zv);
        
        # Assign outputs
        self.X = xm;
        self.Y = ym;
        self.Z = zm;
            
class Flow:
    def __init__(self, Flow_Type = None):
        
        if Flow_Type is None:
            print("Error: Specify a flow type!")
            return;
        
        # Allocate a domain        
        self.Domain = Domain();
        
        # Allocate a velocity object
        self.Velocity = Velocity();
        
        # Save the flow type
        self.FlowType = Flow_Type.lower();

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
    
def HamaVelocity(xy, t, my_args):

    # Amplitude (0.2 seems to work)
    a = my_args[0];
    
    # Wave number (unity in the paper)
    alpha = my_args[1]
    
    # Wave velocity (unity in the paper)
    c = my_args[2];
    
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
    
def UniformVelocity(xy, t, my_args):
    
    # Parse the input
    x, y = Parse_Vector_2d(xy);
    
    num_points = len(x);
    
    # Velocity
    u_o = my_args * np.ones((num_points,), dtype = np.float);
    v_o = 0 * np.ones((num_points,), dtype = np.float);
    
    u_vel = u_o;
    v_vel = v_o;
    
    vels = np.append(u_vel, v_vel);
    
    return vels;
         
    
def plot_pathlines(flow_type = None, t = [0], y0 = [0], my_args = None):
    if flow_type is not None:
        
        fig = plt.figure();
        
        if "hama" in flow_type.lower():
            # Extract the constants related to the Hama flow
                        
            # pdb.set_trace();
                        
            # Integrate
            xy = odeint(HamaVelocity, y0 = y0, t = t, args = my_args);
            ax = plt.axes(xlim=(-1, 9), ylim=(-2.0, 2.0))
            
            
        elif "uniform" in flow_type.lower():
            
            xy = odeint(UniformVelocity, y0 = y0, t = t, args = my_args);
            ax = plt.axes(xlim=(-1, 9), ylim=(-2, 2))
            
        # Count number of path lines
        num_pathlines = int(xy.shape[1] / 2);
        
        # Do the plotting
        for k in range(num_pathlines):
            print("Figure: " + str(k));
            x, y = Get_Pathline(xy, k);
            ax.plot(x, y, '-k');
            
        plt.show(ax)
            
            
            # Make a plot
            # plt.plot(xy, '.');
            # plt.show();













