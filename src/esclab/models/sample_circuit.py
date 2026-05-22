from esclab.simulate import *
from eeslib.functions import convert, converttemp

class CircuitElement(Component):
    """Base class for circuit elements."""
    u_in = Component.Input()   # voltage
    i_in = Component.Input()   # current

    u_out = Component.Output()  # voltage
    i_out = Component.Output()  # current

# Series RLC circuit
# --------------------------
# |        ---->           |
# |                       \\\
# |                  R    ///
# |                       \\\
# xx +                    ///
# xxx   U                  |
# xx -                     |
# |                       |i|
# |                  L    |i|
# |                       |i|
# |                        |
# |                        |
# |                       ---
# |                  C     c
# |                       ---
# |                        |
# |      GGGGG             |
# --------GGG---------------


class VoltageSource(CircuitElement):
    """Voltage source element."""
    V = Component.Parameter()  # voltage

    def __init__(self):
        super().__init__()

    def calculate(self):
        self.u_out = self.V
        # Closed-form BE update for series RLC current (within one timestep):
        # This relies on the specific configuration and is not robust to changes in the 
        # circuit topology, but it is very stable and efficient for this simple example.
        # i = (V + (L/dt)*I_prev - Uc_prev) / (R + L/dt + dt/C)
        # This removes oscillatory fixed-point behavior from the source law.
        dt = self.model.settings.timestep
        R = self.model.r.R
        L = self.model.l.L
        C = self.model.c.C
        I_prev = self.model.l.I_L_prev
        Uc_prev = self.model.c.U_C_prev

        denom = R + (L / dt) + (dt / C)
        denom = np.sign(denom) * max(abs(denom), 1e-12)
        i_target = (self.V + (L / dt) * I_prev - Uc_prev) / denom

        # Mild relaxation avoids tiny iteration jitter with aggressive connection learn rates.
        alpha = 0.8
        self.i_out = float(self.i_in) + alpha * (i_target - float(self.i_in))
        return 

class Resistor(CircuitElement):
    """Resistor element."""
    R = Component.Parameter()  # resistance

    def calculate(self):
        self.u_out = self.u_in - self.i_in * self.R
        self.i_out = self.i_in
        return 

class Capacitor(CircuitElement):
    """Capacitor element."""
    C = Component.Parameter()  # capacitance
    U_C0 = Component.Parameter(0.)  # initial voltage across the capacitor

    def presim_setup(self, **kwargs):
        self.U_C_prev = self.U_C0  # stored capacitor voltage from previous timestep

    def calculate(self):
        # Backward Euler: U_C[n] = U_C[n-1] + i_in * dt / C
        U_C = self.U_C_prev + self.i_in * self.model.settings.timestep / self.C
        self.u_out = self.u_in - U_C
        self.i_out = self.i_in  # series element: current passes through

        if self.model.is_converged:
            self.U_C_prev = U_C
        return 

class Inductor(CircuitElement):
    """Inductor element."""
    L = Component.Parameter()  # inductance
    I_L0 = Component.Parameter(0.)  # initial current through the inductor

    def presim_setup(self, **kwargs):
        self.I_L_prev = self.I_L0  # stored inductor current from previous timestep

    def calculate(self):
        # Backward Euler: V_L = L * (i_in[n] - I_L[n-1]) / dt
        didt = (self.i_in - self.I_L_prev) / self.model.settings.timestep
        self.u_out = self.u_in - self.L * didt
        self.i_out = self.i_in  # series element: current passes through

        if self.model.is_converged:
            self.I_L_prev = self.i_in
        return 

class TeeOut(Component):
    """Tee element for splitting voltage and current."""
    u_in = Component.Input()   # voltage
    i_in = Component.Input()   # current
    u_branch_1 = Component.Input()  # voltage at the terminus of branch 1
    u_branch_2 = Component.Input()  # voltage at the terminus of branch 2

    u_out = Component.Output()  # voltage output
    i_out_1 = Component.Output()  # current output 1
    i_out_2 = Component.Output()  # current output 2

    def calculate(self):
        self.u_out = self.u_in
        dU1 = max(0., float(self.u_in - self.u_branch_1))
        dU2 = max(0., float(self.u_in - self.u_branch_2))
        total_dU = dU1 + dU2
        f = dU1 / total_dU if total_dU > 1e-12 else 0.5  # equal split when voltages are indeterminate
        self.i_out_1 = self.i_in * f
        self.i_out_2 = self.i_in - self.i_out_1
        return

class TeeReturn(Component):
    """Tee element for combining voltage and current."""
    u_in_1 = Component.Input()  # voltage from branch 1
    i_in_1 = Component.Input()  # current from branch 1
    u_in_2 = Component.Input()  # voltage from branch 2
    i_in_2 = Component.Input()  # current from branch 2

    u_out = Component.Output()  # voltage output
    i_out = Component.Output()  # current output

    def calculate(self):
        self.i_out = self.i_in_1 + self.i_in_2
        i_total = abs(float(self.i_out))
        if i_total > 1e-12:
            self.u_out = (self.u_in_1 * self.i_in_1 + self.u_in_2 * self.i_in_2) / float(self.i_out)
        else:
            self.u_out = (self.u_in_1 + self.u_in_2) / 2
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

# Intermediate parameters
zeta = 0.5  # damping ratio

# Set any parameters
model.vs.V = 1
model.l.L = 1
model.c.C = 1
model.r.R = 2 * zeta * np.sqrt(model.l.L / model.c.C)  # critical damping

# Set initial values
model.r.u_in = 0
model.l.u_in = 0
model.c.u_in = 0
model.vs.u_in = 0
model.r.i_in = 0
model.l.i_in = 0
model.c.i_in = 0
model.vs.i_in = 0

# Configure simulation settings
model.settings.start_time = 0
model.time = model.settings.start_time
model.settings.stop_time = 15  #seconds
model.settings.timestep = 5e-2  # seconds
model.settings.max_iterations = 50
# model.settings.tol_rel_global = 1e-12
# model.settings.learn_rate = .98

# ----------------------------------------------------------
# Make connections
# ----------------------------------------------------------

# Voltages
model.connect( model.vs.u_out, model.r.u_in)
model.connect( model.r.u_out, model.l.u_in)
model.connect( model.l.u_out, model.c.u_in)
model.connect( model.c.u_out, model.vs.u_in)

# Currents
model.connect( model.vs.i_out, model.r.i_in)
model.connect( model.r.i_out, model.l.i_in)
model.connect( model.l.i_out, model.c.i_in)
model.connect( model.c.i_out, model.vs.i_in)

# ----------------------------------------------------------

# Set up the simulation
model.add_plotter([model.vs.u_out, model.r.u_out, model.l.u_out, model.c.u_out, model.vs.u_in],
                  [model.vs.i_out, model.r.i_out, model.l.i_out, model.c.i_out, model.vs.i_in], 
                  y1label="Voltage (V)", 
                  y2label="Current (A)", 
                  update_every=10, 
                  nmax_points=1000)


while model.time < model.settings.stop_time:
    model.step()

model.wait_for_plots()
