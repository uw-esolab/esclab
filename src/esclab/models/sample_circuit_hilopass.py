"""
This example presents an RLC circuit configured as a hi/lo pass filter. This type of filter 
removes high or low frequency components from an input signal, allowing only the desired 
frequency range to appear at the resistive load. 

To call the circuit a "hi-pass" or "low-pass" filter, we can swap the positions of the 
inductor and capacitor. Use the is_hipass variable to select which configuration to simulate.

From a simulation perspective, note that the system is solved using the matrix inversion
capability, and the voltage signal in VoltageSource is time-varying and recalculated at
each timestep. The time-varying voltage is computed inside calculate() when self.context
is not None (the matrix-build phase), so the matrix always uses the current timestep value.
"""

# Hi pass or low-pass RLC circuit 
# 
# Hi pass:
#       ---------- 
# ------LLLLLLLLLL------1----------
# |     ----------    |2          |
# |         L         |           |
# |                  |C|         \\\
# xxx   U         C  |C|     R   ///
# xx +               |C|         \\\
# xx -               |C|         ///
# |      GGGGG        |           |
# --------GGG---------TR-----------
# 
# Low pass:
#       ---------- 
# ------CCCCCCCCCC------1----------
# |     ----------    |2          |
# |         C         |           |
# |                  |L|         \\\
# xxx   U         L  |L|     R   ///
# xx +               |L|         \\\
# xx -               |L|         ///
# |      GGGGG        |           |
# --------GGG-----TR---------------

from esclab.simulate import *
from esclab.components.circuit_elements import *

# ############################################################################
#                   Set up the model and circuit
# ############################################################################

# Set up the circuit
model = Model()

# Create components
model.vs = VoltageSource()
model.r = Resistor()
model.c = Capacitor()
model.l = Inductor()
model.tee_out = TeeOut()
model.tee_return = TeeReturn()

# ----------------------------------------------------------
# Initialize
model.initialize()

# Intermediate parameters
zeta = 0.8  # damping ratio
is_hipass = True
# is_hipass = False


# ----------------------------------------------------------
# Set any parameters
model.vs.V = 1
model.l.L = .5
model.c.C = 1
model.r.R = 2 * zeta * np.sqrt(model.l.L / model.c.C)  # critical damping

# filter cutoff frequency
omega_c = 1 / np.sqrt(model.l.L * model.c.C)
print("Filter cutoff frequency (rad/s):", omega_c)

model.vs.omega_hi = omega_c * 0.05  # high-pass cutoff at half the natural frequency
model.vs.omega_lo = omega_c * 2    # low-pass cutoff at twice the natural frequency

# Switch model type here...
if is_hipass:
    # Hi-pass tolopolgy. Inductor in series, capacitor in parallel
    pos_A = model.l
    pos_B = model.c
else:    
    # Low-pass tolopolgy. Capacitor in series, inductor in parallel
    pos_A = model.c
    pos_B = model.l

# ----------------------------------------------------------
# Configure simulation settings
model.settings.stop_time = 25  #seconds
model.settings.timestep = .02  # seconds
model.settings.max_iterations = 50
model.settings.tol_rel_global = 1e-6

# ----------------------------------------------------------
# Make connections

# Voltages
model.connect( model.vs.u_out           , pos_A.u_in              , solve_group="voltage")
model.connect( pos_A.u_out              , model.tee_out.u_in      , solve_group="voltage")
model.connect( model.tee_out.u_out      , model.r.u_in            , solve_group="voltage")
model.connect( model.tee_out.u_out      , pos_B.u_in              , solve_group="voltage")
model.connect( model.r.u_out            , model.tee_return.u_in_1 , solve_group="voltage")
model.connect( pos_B.u_out              , model.tee_return.u_in_2 , solve_group="voltage")
model.connect( model.tee_return.u_out   , model.vs.u_in           , solve_group="voltage")

# Currents
model.connect( model.vs.i_out           , pos_A.i_in              , solve_group="current")
model.connect( pos_A.i_out              , model.tee_out.i_in      , solve_group="current")
model.connect( model.tee_out.i_out_1    , model.r.i_in            , solve_group="current")
model.connect( model.tee_out.i_out_2    , pos_B.i_in              , solve_group="current")
model.connect( model.r.i_out            , model.tee_return.i_in_1 , solve_group="current")
model.connect( pos_B.i_out              , model.tee_return.i_in_2 , solve_group="current")
model.connect( model.tee_return.i_out   , model.vs.i_in           , solve_group="current")

# ----------------------------------------------------------
# Set up the plotters

model.add_plotter([model.vs.u_out, model.r.u_out, model.l.u_out, model.c.u_out, model.vs.u_in],
                  y1label="Voltage (V)", 
                  update_every=5,
                  nmax_points=1000,
                  )
model.add_plotter([model.vs.i_out, model.r.i_out, model.l.i_out, model.c.i_out, model.vs.i_in], 
                  y1label="Current (A)", 
                  update_every=5,
                  nmax_points=1000,
                  )

model.add_network_graph()

# ----------------------------------------------------------
# Run the simulation
while model.time < model.settings.stop_time:
    model.step()

model.wait_for_plots()
