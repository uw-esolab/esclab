"""
This example demonstrates how to solve the simple RLC circuit example 
using a set of linear equations with matrix inversion. The matrix
is constructed using dedicated equation collectors inside each component, 
and this means the system of equation defining flow and potential 
relationships is explicitly defined and solved at each timestep. This is 
in contrast to the first example in which the current was determined within
the voltage source component using a bracketed solver technique. 

In principle, the matrix inversion method should be more reliable and 
faster than the bracketed solver approach, but it is also more complex to
set up.

Series RLC circuit
--------------------------
|        ---->           |
|                       \\\
|                  R    ///
|                       \\\
xx +                    ///
xxx   U                  |
xx -                     |
|                       |i|
|                  L    |i|
|                       |i|
|                        |
|                       ---
|                  C     c
|                       ---
|      GGGGG             |
--------GGG---------------
"""

from esclab.simulate import *
from esclab.components.circuit_elements import *


# ############################################################################
#                   Set up the model and circuit
# ############################################################################
model = Model()

# Create components
model.vs = VoltageSource()
model.r = Resistor()
model.c = Capacitor()
model.l = Inductor()

# ----------------------------------------------------------
# Initialize
model.initialize()

# Choose a damping ratio
# https://en.wikipedia.org/wiki/RLC_circuit#/media/File:RLC_transient_plot.svg
zeta = 0.5  # damping ratio

# ----------------------------------------------------------
# Set any parameters
model.vs.V = 1
model.l.L = 1
model.c.C = 1
model.r.R = 2 * zeta * np.sqrt(model.l.L / model.c.C)  # critical damping

# ----------------------------------------------------------
# Configure simulation settings
model.settings.start_time = 0
model.time = model.settings.start_time
model.settings.stop_time = 15  #seconds
model.settings.timestep = 5e-2  # seconds
model.settings.max_iterations = 100
model.settings.tol_rel_global = 1e-6
model.settings.tol_abs_global = 1  # relax absolute tolerance - matrix solve is more stable

# ----------------------------------------------------------
# Make connections

# Voltages
model.connect( model.vs.u_out, model.r.u_in,  solve_group="voltage")
model.connect( model.r.u_out,  model.l.u_in,  solve_group="voltage")
model.connect( model.l.u_out,  model.c.u_in,  solve_group="voltage")
model.connect( model.c.u_out,  model.vs.u_in, solve_group="voltage")

# Currents
model.connect( model.vs.i_out, model.r.i_in,  solve_group="current")
model.connect( model.r.i_out,  model.l.i_in,  solve_group="current")
model.connect( model.l.i_out,  model.c.i_in,  solve_group="current")
model.connect( model.c.i_out,  model.vs.i_in, solve_group="current")

# ----------------------------------------------------------
# Set up the plotters

model.add_plotter([model.vs.u_out, model.r.u_out, model.l.u_out, model.c.u_out, model.vs.u_in],
                  y1label="Voltage (V)", 
                  update_every=4, 
                  nmax_points=1000)
model.add_plotter([model.vs.i_out, model.r.i_out, model.l.i_out, model.c.i_out, model.vs.i_in], 
                  y1label="Current (A)", 
                  update_every=4, 
                  nmax_points=1000)

model.add_network_graph()

# ----------------------------------------------------------
# Run the simulation

while model.time < model.settings.stop_time:
    model.step()

model.wait_for_plots()
