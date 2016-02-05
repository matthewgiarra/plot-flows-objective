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




    


#
# def plot_streaklines(flow_type = None, t = [0], y0 = [0], my_args = None):
#     if flow_type is not None:
#
#         # Size of the time step
#         time_step_size = 0.1;
#
#         # Number of flow tracer sources
#         num_sources = int(len(y0) / 2);
#
#
#
#         # List of X coordinates of all the particles
#         # that will be tracked
#
#
#         # Allocate all of the coordinate positions
#         # so we don't have to re-allocate memory all the time
#
#
#
#
#
#
# def plot_pathlines(flow_type = None, t = [0], y0 = [0], my_args = None):
#     if flow_type is not None:
#
#         fig = plt.figure();
#
#         if "hama" in flow_type.lower():
#             # Extract the constants related to the Hama flow
#
#             # pdb.set_trace();
#
#             # Integrate
#             xy = odeint(HamaVelocity, y0 = y0, t = t, args = (my_args,));
#             ax = plt.axes(xlim=(-1, 9), ylim=(-2.0, 2.0))
#
#
#         elif "uniform" in flow_type.lower():
#
#             xy = odeint(UniformVelocity, y0 = y0, t = t, args = (my_args,));
#             ax = plt.axes(xlim=(-1, 9), ylim=(-2, 2))
#
#         # Count number of path lines
#         num_pathlines = int(xy.shape[1] / 2);
#
#         # Do the plotting
#         for k in range(num_pathlines):
#             print("Figure: " + str(k));
#             x, y = Get_Pathline(xy, k);
#             ax.plot(x, y, '-k');
#
#         plt.show(ax)
#










