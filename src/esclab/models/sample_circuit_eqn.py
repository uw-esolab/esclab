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

class CircuitElement(Component):
    """Base class for circuit elements."""
    u_in = Component.Input()   # voltage
    i_in = Component.Input()   # current

    u_out = Component.Output()  # voltage
    i_out = Component.Output()  # current

    def get_network_equations(self, context):
        # All series elements pass current straight through: i_out = i_in
        # This must be called in child classes using super().get_network_equations(context) 
        context.add_equation({
            self.i_out: 1.0,
            context.source(self.i_in): -1.0,
        }, rhs=0.0)

class VoltageSource(CircuitElement):
    """Ideal voltage source element."""
    V = Component.Parameter()  # source voltage

    def get_network_equations(self, context):
        # Set the voltage rise across the source: u_out - u_in = V.
        context.add_equation({
            self.u_out: 1.0,
            context.source(self.u_in): -1.0,
        }, rhs=self.V)
        # Ground the source input node so the circuit has a fixed 0 V reference.
        context.add_equation({
            context.source(self.u_in): 1.0,
        }, rhs=0.0)
        super().get_network_equations(context)  # current pass-through

    def calculate(self):
        return 

class Resistor(CircuitElement):
    """Resistor element."""
    R = Component.Parameter()  # resistance

    def get_network_equations(self, context):
        super().get_network_equations(context)  # current pass-through
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

    def get_network_equations(self, context):
        super().get_network_equations(context)  # current pass-through
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

    def get_network_equations(self, context):
        super().get_network_equations(context)  # current pass-through
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

# Set up the circuit
model = Model()

# Create components
model.vs = VoltageSource()
model.r = Resistor()
model.c = Capacitor()
model.l = Inductor()

# Initialize
model.initialize()

# Choose a damping ratio
# https://en.wikipedia.org/wiki/RLC_circuit#/media/File:RLC_transient_plot.svg
zeta = 0.5  # damping ratio

# Set component values
model.vs.V = 1
model.l.L = 1
model.c.C = 1
model.r.R = 2 * zeta * np.sqrt(model.l.L / model.c.C)  # critical damping

# Set simulation options
model.settings.start_time = 0
model.time = model.settings.start_time
model.settings.stop_time = 15  #seconds
model.settings.timestep = 5e-2  # seconds
model.settings.max_iterations = 100
model.settings.tol_rel_global = 1e-6
model.settings.tol_abs_global = 1  # relax absolute tolerance - matrix solve is more stable

# Connect the circuit

# Voltages
model.connect( model.vs.u_out, model.r.u_in,  solve_group="potential")
model.connect( model.r.u_out,  model.l.u_in,  solve_group="potential")
model.connect( model.l.u_out,  model.c.u_in,  solve_group="potential")
model.connect( model.c.u_out,  model.vs.u_in, solve_group="potential")

# Currents
model.connect( model.vs.i_out, model.r.i_in,  solve_group="flow")
model.connect( model.r.i_out,  model.l.i_in,  solve_group="flow")
model.connect( model.l.i_out,  model.c.i_in,  solve_group="flow")
model.connect( model.c.i_out,  model.vs.i_in, solve_group="flow")

# Plot the results
model.add_plotter([model.vs.u_out, model.r.u_out, model.l.u_out, model.c.u_out, model.vs.u_in],
                  y1label="Voltage (V)", 
                  update_every=4, 
                  nmax_points=1000)
model.add_plotter([model.vs.i_out, model.r.i_out, model.l.i_out, model.c.i_out, model.vs.i_in], 
                  y1label="Current (A)", 
                  update_every=4, 
                  nmax_points=1000)


while model.time < model.settings.stop_time:
    model.step()

model.wait_for_plots()
