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

    def add_network_equations(self, context):
        # Set the voltage rise across the source: u_out - u_in = V.
        context.add_equation({
            self.u_out: 1.0,
            context.source(self.u_in): -1.0,
        }, rhs=self.V)
        # Ground the source input node so the circuit has a fixed 0 V reference.
        context.add_equation({
            context.source(self.u_in): 1.0,
        }, rhs=0.0)
        super().add_network_equations(context)  # current pass-through

    def calculate(self):
        return 

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
model.r1 = Resistor()
model.r1.name = "Resistor_1"
model.c = Capacitor()
model.l = Inductor()
model.r2 = Resistor()
model.r2.name = "Resistor_2"
model.tee_out = TeeOut()
model.tee_return = TeeReturn()

# Initialize
model.initialize()

# Intermediate parameters
zeta = 0.8  # damping ratio

# Set any parameters
model.vs.V = 1
model.l.L = 1
model.c.C = 1
model.r2.R = 2 * zeta * np.sqrt(model.l.L / model.c.C)  # critical damping
model.r1.R = model.r2.R/10  

# Configure simulation settings
model.settings.start_time = 0
model.time = model.settings.start_time
model.settings.stop_time = 25  #seconds
model.settings.timestep = .05  # seconds
model.settings.max_iterations = 50
model.settings.tol_rel_global = 1e-6

# ----------------------------------------------------------
# Make connections
# ----------------------------------------------------------

# Voltages
model.connect( model.vs.u_out           , model.tee_out.u_in      , solve_group="potential")
model.connect( model.tee_out.u_out      , model.r1.u_in           , solve_group="potential")
model.connect( model.tee_out.u_out      , model.r2.u_in           , solve_group="potential")
model.connect( model.r1.u_out           , model.l.u_in            , solve_group="potential")
model.connect( model.l.u_out            , model.c.u_in            , solve_group="potential")
model.connect( model.c.u_out            , model.tee_return.u_in_1 , solve_group="potential")
model.connect( model.r2.u_out           , model.tee_return.u_in_2 , solve_group="potential")
model.connect( model.tee_return.u_out   , model.vs.u_in           , solve_group="potential")

# Currents
model.connect( model.vs.i_out           , model.tee_out.i_in      , solve_group="flow")
model.connect( model.tee_out.i_out_1    , model.r1.i_in           , solve_group="flow")
model.connect( model.tee_out.i_out_2    , model.r2.i_in           , solve_group="flow")
model.connect( model.r1.i_out           , model.l.i_in            , solve_group="flow")
model.connect( model.l.i_out            , model.c.i_in            , solve_group="flow")
model.connect( model.c.i_out            , model.tee_return.i_in_1 , solve_group="flow")
model.connect( model.r2.i_out           , model.tee_return.i_in_2 , solve_group="flow")
model.connect( model.tee_return.i_out   , model.vs.i_in           , solve_group="flow")

# ----------------------------------------------------------

# Set up the simulation
model.add_plotter([model.vs.u_out, model.r1.u_out, model.r2.u_out, model.l.u_out, model.c.u_out, model.vs.u_in],
                  y1label="Voltage (V)", 
                  update_every=5, 
                  nmax_points=1000)
model.add_plotter([model.vs.i_out, model.r1.i_out, model.r2.i_out, model.l.i_out, model.c.i_out, model.vs.i_in], 
                  y1label="Current (A)", 
                  update_every=5, 
                  nmax_points=1000)

model.add_network_graph()

while model.time < model.settings.stop_time:
    model.step()

model.wait_for_plots()
