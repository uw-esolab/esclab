"""
Example of an RLC circuit with a parallel resistor branch. 
The circuit is solved using the matrix inversion method.
"""

# RLC circuit with parallel resistor branch
# ----------------TO--1-----------
# |        ---->   |             |
# |                |2           \\\
# |                |       R1   ///
# |                |            \\\
# xx +         R2 \\\           ///
# xxx   U         ///            |
# xx -            \\\            |
# |               ///           |i|
# |                |       L    |i|
# |                |            |i|
# |                |             |
# |                |            ---
# |                |       C     c
# |                |            ---
# |      GGGGG     |             |
# --------GGG-----TR--------------
from esclab.simulate import *
from esclab.components.circuit_elements import *

# ############################################################################
#                   Set up the model and circuit
# ############################################################################

# Set up the circuit
model = Model()

# Create components
model.vs = VoltageSource()
model.r1 = Resistor()
model.r1.name = "Resistor_1"
model.c = Capacitor()
model.l = Inductor()
model.r2 = Resistor()
model.r2.name = "Resistor_2"
model.tee_out = TeeOut()
model.tee_return = TeeReturn()

# ----------------------------------------------------------
# Initialize
model.initialize()

# ----------------------------------------------------------
# Set any parameters
zeta = 0.8  # damping ratio

# Set any parameters
model.vs.V = 1
model.l.L = 1
model.c.C = 1
model.r2.R = 2 * zeta * np.sqrt(model.l.L / model.c.C)  # critical damping
model.r1.R = model.r2.R/10  

# ----------------------------------------------------------
# Configure simulation settings
model.settings.start_time = 0
model.time = model.settings.start_time
model.settings.stop_time = 25  #seconds
model.settings.timestep = .05  # seconds
model.settings.max_iterations = 50
model.settings.tol_rel_global = 1e-6

# ----------------------------------------------------------
# Make connections

# Voltages
model.connect( model.vs.u_out           , model.tee_out.u_in      , solve_group="voltage")
model.connect( model.tee_out.u_out      , model.r1.u_in           , solve_group="voltage")
model.connect( model.tee_out.u_out      , model.r2.u_in           , solve_group="voltage")
model.connect( model.r1.u_out           , model.l.u_in            , solve_group="voltage")
model.connect( model.l.u_out            , model.c.u_in            , solve_group="voltage")
model.connect( model.c.u_out            , model.tee_return.u_in_1 , solve_group="voltage")
model.connect( model.r2.u_out           , model.tee_return.u_in_2 , solve_group="voltage")
model.connect( model.tee_return.u_out   , model.vs.u_in           , solve_group="voltage")

# Currents
model.connect( model.vs.i_out           , model.tee_out.i_in      , solve_group="current")
model.connect( model.tee_out.i_out_1    , model.r1.i_in           , solve_group="current")
model.connect( model.tee_out.i_out_2    , model.r2.i_in           , solve_group="current")
model.connect( model.r1.i_out           , model.l.i_in            , solve_group="current")
model.connect( model.l.i_out            , model.c.i_in            , solve_group="current")
model.connect( model.c.i_out            , model.tee_return.i_in_1 , solve_group="current")
model.connect( model.r2.i_out           , model.tee_return.i_in_2 , solve_group="current")
model.connect( model.tee_return.i_out   , model.vs.i_in           , solve_group="current")

# ----------------------------------------------------------
# Set up the plotters

model.add_plotter([model.vs.u_out, model.r1.u_out, model.r2.u_out, model.l.u_out, model.c.u_out, model.vs.u_in],
                  y1label="Voltage (V)", 
                  update_every=5, 
                  nmax_points=1000)
model.add_plotter([model.vs.i_out, model.r1.i_out, model.r2.i_out, model.l.i_out, model.c.i_out, model.vs.i_in], 
                  y1label="Current (A)", 
                  update_every=5, 
                  nmax_points=1000)

model.add_network_graph()

# ----------------------------------------------------------
# Run the simulation
while model.time < model.settings.stop_time:
    model.step()

model.wait_for_plots()
