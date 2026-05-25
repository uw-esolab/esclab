"""
This example presents an RLC circuit configured as a hi/lo pass filter. This type of filter 
removes high or low frequency components from an input signal, allowing only the desired 
frequency range to appear at the resistive load. 

To call the circuit a "hi-pass" or "low-pass" filter, we can swap the positions of the 
inductor and capacitor. Use the is_hipass variable to select which configuration to simulate.

From a simulation perspective, note that the system is again solved using the matrix 
inversion capability, and the voltage signal in VoltageSource is time-varying and recalculated
at each timestep. To get this to work properly, the calculations should be carried out inside 
the add_network_equations() method of VoltageSource rather than in calculate(), since the 
voltage value is needed during the network solve to determine the system state at each timestep.
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

class CircuitElement(Component):
    """Base class for circuit elements."""
    u_in = Component.Input()   # voltage
    i_in = Component.Input()   # current

    u_out = Component.Output()  # voltage
    i_out = Component.Output()  # current

    def add_network_equations(self, context):
        # All series elements pass current straight through: i_out = i_in
        # This must be called in child classes using super().add_network_equations(context) 
        context.add_equation({
            self.i_out: 1.0,
            context.source(self.i_in): -1.0,
        }, rhs=0.0)

class VoltageSource(CircuitElement):
    """Ideal voltage source element."""
    V = Component.Parameter()  # source voltage
    omega_hi = Component.Parameter(0.1)  # cutoff frequency for high-pass behavior
    omega_lo = Component.Parameter(2)  # cutoff frequency for low-pass behavior

    def __init__(self):
        super().__init__()
        self.V_t = self.V  # stored voltage for time-varying sources

    def add_network_equations(self, context):
        super().add_network_equations(context)  # current pass-through
        # Set the voltage rise across the source: u_out - u_in = V.
        V_t = 1
        V_t += np.sin(2 * np.pi * self.omega_hi * self.model.time) * 0.25   #Oscillates at the high-pass cutoff frequency
        V_t += np.sin(2 * np.pi * self.omega_lo * self.model.time) * 0.25  #Oscillates at the low-pass cutoff frequency
        V_t += self.V * V_t
        context.add_equation({
            self.u_out: 1.0,
            context.source(self.u_in): -1.0,
        }, rhs=V_t)
        # Ground the source input node so the circuit has a fixed 0 V reference.
        context.add_equation({
            context.source(self.u_in): 1.0,
        }, rhs=0.0)


class Resistor(CircuitElement):
    """Resistor element."""
    R = Component.Parameter()  # resistance

    def add_network_equations(self, context):
        super().add_network_equations(context)  # current pass-through
        # u_out = u_in - R*i_in    u_out - u_in + R*i_in = 0
        context.add_equation({
            self.u_out: 1.0,
            context.source(self.u_in): -1.0,
            context.source(self.i_in): self.R,
        }, rhs=0.0)

    def calculate(self):
        # The network solve already sets the resistor voltage and current.
        return 

class Capacitor(CircuitElement):
    """Capacitor element."""
    C = Component.Parameter()  # capacitance
    U_C0 = Component.Parameter(0.)  # initial voltage across the capacitor

    def presim_setup(self, **kwargs):
        self.U_C_prev = self.U_C0  # stored capacitor voltage from previous timestep

    def add_network_equations(self, context):
        super().add_network_equations(context)  # current pass-through
        # Backward Euler: u_out = u_in - U_C_prev - (dt/C)*i_in
        # u_out - u_in + (dt/C)*i_in = -U_C_prev
        dt = self.model.settings.timestep
        context.add_equation({
            self.u_out: 1.0,
            context.source(self.u_in): -1.0,
            context.source(self.i_in): dt / self.C,
        }, rhs=-self.U_C_prev)

    def calculate(self):
        if self.model.is_converged:
            # After convergence, store the capacitor voltage for the next timestep.
            U_C = self.U_C_prev + self.i_in * self.model.settings.timestep / self.C
            self.U_C_prev = U_C
        return 

class Inductor(CircuitElement):
    """Inductor element."""
    L = Component.Parameter()  # inductance
    I_L0 = Component.Parameter(0.)  # initial current through the inductor

    def presim_setup(self, **kwargs):
        self.I_L_prev = self.I_L0  # stored inductor current from previous timestep

    def add_network_equations(self, context):
        super().add_network_equations(context)  # current pass-through
        # Backward Euler: u_out = u_in - (L/dt)*i_in + (L/dt)*I_L_prev
        # u_out - u_in + (L/dt)*i_in = (L/dt)*I_L_prev
        dt = self.model.settings.timestep
        L = self.L
        context.add_equation({
            self.u_out: 1.0,
            context.source(self.u_in): -1.0,
            context.source(self.i_in): L / dt,
        }, rhs=(L / dt) * self.I_L_prev)

    def calculate(self):
        if self.model.is_converged:
            # After convergence, store the inductor current for the next timestep.
            self.I_L_prev = self.i_in
        return 


class TeeOut(Component):
    """Tee element for splitting voltage and current."""
    u_in = Component.Input()   # voltage
    i_in = Component.Input()   # current

    u_out = Component.Output()  # voltage output
    i_out_1 = Component.Output()  # current output 1
    i_out_2 = Component.Output()  # current output 2

    def calculate(self):
        return 
    
    def add_network_equations(self, context):
        # Current splits: i_in = i_out_1 + i_out_2
        context.add_equation({
            context.source(self.i_in): 1.0,
            self.i_out_1: -1.0,
            self.i_out_2: -1.0,
        }, rhs=0.0)
        # Output voltage is the same as input voltage: u_out = u_in
        context.add_equation({
            self.u_out: 1.0,
            context.source(self.u_in): -1.0,
        }, rhs=0.0)

class TeeReturn(Component):
    """Tee element for combining voltage and current."""
    u_in_1 = Component.Input()  # voltage from branch 1
    i_in_1 = Component.Input()  # current from branch 1
    u_in_2 = Component.Input()  # voltage from branch 2
    i_in_2 = Component.Input()  # current from branch 2

    u_out = Component.Output()  # voltage output
    i_out = Component.Output()  # current output

    def calculate(self):
        return
    
    def add_network_equations(self, context):
        # Current combines: i_out = i_in_1 + i_in_2
        context.add_equation({
            self.i_out: 1.0,
            context.source(self.i_in_1): -1.0,
            context.source(self.i_in_2): -1.0,
        }, rhs=0.0)
        # Output voltage is the same as input voltages: u_out = u_in_1 = u_in_2
        context.add_equation({
            self.u_out: 1.0,
            context.source(self.u_in_1): -1.0,
        }, rhs=0.0)
        context.add_equation({
            self.u_out: 1.0,
            context.source(self.u_in_2): -1.0,
        }, rhs=0.0)

# Set up the circuit
model = Model()

# Create components
model.vs = VoltageSource()
model.r = Resistor()
model.c = Capacitor()
model.l = Inductor()
model.tee_out = TeeOut()
model.tee_return = TeeReturn()

# Initialize
model.initialize()

# Intermediate parameters
zeta = 0.8  # damping ratio
is_hipass = True
# is_hipass = False


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

# Configure simulation settings
model.settings.stop_time = 25  #seconds
model.settings.timestep = .02  # seconds
model.settings.max_iterations = 50
model.settings.tol_rel_global = 1e-6

# ----------------------------------------------------------
# Make connections
# ----------------------------------------------------------

if is_hipass:
    pos_A = model.l
    pos_B = model.c
else:    
    pos_A = model.c
    pos_B = model.l

# Voltages
model.connect( model.vs.u_out           , pos_A.u_in              , solve_group="potential")
model.connect( pos_A.u_out              , model.tee_out.u_in      , solve_group="potential")
model.connect( model.tee_out.u_out      , model.r.u_in            , solve_group="potential")
model.connect( model.tee_out.u_out      , pos_B.u_in              , solve_group="potential")
model.connect( model.r.u_out            , model.tee_return.u_in_1 , solve_group="potential")
model.connect( pos_B.u_out              , model.tee_return.u_in_2 , solve_group="potential")
model.connect( model.tee_return.u_out   , model.vs.u_in           , solve_group="potential")

# Currents
model.connect( model.vs.i_out           , pos_A.i_in              , solve_group="flow")
model.connect( pos_A.i_out              , model.tee_out.i_in      , solve_group="flow")
model.connect( model.tee_out.i_out_1    , model.r.i_in            , solve_group="flow")
model.connect( model.tee_out.i_out_2    , pos_B.i_in              , solve_group="flow")
model.connect( model.r.i_out            , model.tee_return.i_in_1 , solve_group="flow")
model.connect( pos_B.i_out              , model.tee_return.i_in_2 , solve_group="flow")
model.connect( model.tee_return.i_out   , model.vs.i_in           , solve_group="flow")

# ----------------------------------------------------------

# Set up the simulation
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

while model.time < model.settings.stop_time:
    model.step()

model.wait_for_plots()
