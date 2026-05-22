from esclab.simulate import *
from eeslib.functions import convert, converttemp

class CircuitElement(Component):
    """Base class for circuit elements."""
    u_in = Component.Input()   # voltage
    i_in = Component.Input()   # current

    u_out = Component.Output()  # voltage
    i_out = Component.Output()  # current

#   -------------TO1--1---------------------
#   |            |2                        |
#   |           \\\                       |i|
#   | +         ///  R1                L  |i|
#   xxx   U0    \\\                       |i|
#   xxxx        ///                       |i|
#   xxx          |                         |
#   | -          TO2--1---\/\/\/\/\------TR1
#   |            |2         R2             |
#   |           |c|                        |
#   |           |c|  C                     |
#   |           |c|                        |
#   |            |                         |
#   -------------TR2-----------GGGGG--------
#                               GGG            
#                                G             

class VoltageSource(CircuitElement):
    """Voltage source element."""
    V = Component.Parameter()  # voltage

    def __init__(self):
        super().__init__()

    def calculate(self):
        if self.model.is_first_step:
            self.i_last = .1 #self.i_in.v  # last current

        self.u_out = self.V
        # u_in is the return node voltage; it should converge to 0 (ground).
        # Estimate total circuit resistance from the current operating point and
        # rescale the current so the full voltage V is dropped across the circuit.
        denom = float(self.V - self.u_in)
        if abs(denom) > 1e-12:
            if denom < 0:
                denom = -denom
            self.i_out = self.i_last * (float(self.V) / denom)**.2
        else:
            self.i_out = self.i_in
        self.i_last = self.i_out.v
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
model.r1 = Resistor()
model.r2 = Resistor()
model.r2.name = 'Resistor_2'
model.c = Capacitor()
model.l = Inductor()
model.to1 = TeeOut()
model.tr1 = TeeReturn()
model.to2 = TeeOut()
model.tr2 = TeeReturn()
model.to2.name = 'TeeOut_2'
model.tr2.name = 'TeeReturn_2'

# Initialize
model.initialize()

# Set any parameters
model.vs.V = 10
model.r1.R = 100
model.r2.R = 20
model.c.C = 1e-6
model.l.L = 0.0015

# Set initial values
model.to1.u_in = 0
model.r1.u_in = 0
model.l.u_in = 0
model.to2.u_in = 0
model.c.u_in = 0
model.r2.u_in = 0
model.tr1.u_in_1 = 0
model.tr1.u_in_2 = 0
model.tr2.u_in_1 = 0
model.tr2.u_in_2 = 0
model.vs.u_in = 0
model.to1.i_in = 0
model.l.i_in = 0
model.r1.i_in = 0
model.to2.i_in = 0
model.r2.i_in = 0
model.c.i_in = 0
model.tr1.i_in_1 = 0
model.tr1.i_in_2 = 0
model.tr2.i_in_1 = 0
model.tr2.i_in_2 = 0
model.vs.i_in = 0
# ----------------------------------------------------------
# Make connections
# ----------------------------------------------------------

# Voltages
model.connect( model.vs.u_out   ,  model.to1.u_in        )
model.connect( model.to1.u_out  ,  model.r1.u_in         )
model.connect( model.to1.u_out  ,  model.l.u_in          )
model.connect( model.r1.u_out   ,  model.to2.u_in        )
model.connect( model.to2.u_out  ,  model.c.u_in          )
model.connect( model.to2.u_out  ,  model.r2.u_in         )
model.connect( model.r2.u_out   ,  model.tr1.u_in_1      )
model.connect( model.l.u_out    ,  model.tr1.u_in_2      )
model.connect( model.c.u_out    ,  model.tr2.u_in_1      )
model.connect( model.tr1.u_out  ,  model.tr2.u_in_2      )
model.connect( model.tr2.u_out  ,  model.vs.u_in         )

# Feedback voltages for the tees
model.connect( model.l.u_out    ,  model.to1.u_branch_1  )
model.connect( model.r2.u_out   ,  model.to1.u_branch_2  )
model.connect( model.tr1.u_out  ,  model.to2.u_branch_1  )
model.connect( model.c.u_out    ,  model.to2.u_branch_2  )

# Currents
model.connect( model.vs.i_out   ,  model.to1.i_in        )
model.connect( model.to1.i_out_1,  model.l.i_in          )
model.connect( model.to1.i_out_2,  model.r1.i_in         )
model.connect( model.r1.i_out   ,  model.to2.i_in        )
model.connect( model.to2.i_out_1,  model.r2.i_in         )
model.connect( model.to2.i_out_2,  model.c.i_in          )
model.connect( model.r2.i_out   ,  model.tr1.i_in_1      )
model.connect( model.l.i_out    ,  model.tr1.i_in_2      )
model.connect( model.c.i_out    ,  model.tr2.i_in_1      )
model.connect( model.tr1.i_out  ,  model.tr2.i_in_2      )
model.connect( model.tr2.i_out  ,  model.vs.i_in         )
# ----------------------------------------------------------

# Set up the simulation
model.add_plotter([model.vs.u_out, model.r1.u_out, model.r2.u_out, model.c.u_out, model.l.u_out],
                  [model.vs.i_out, model.r1.i_out, model.r2.i_out, model.c.i_out, model.l.i_out], 
                  y1label="Voltage (V)", 
                  y2label="Current (A)", 
                  update_every=1, 
                  nmax_points=1000)

model.settings.start_time = 0
model.time = model.settings.start_time
model.settings.stop_time = 1e-3  #seconds
model.settings.timestep = 1e-6  # seconds
model.settings.max_iterations = 50
model.settings.tol_rel = 1e-6
model.settings.learn_rate = .7

while model.time < model.settings.stop_time:
    model.step()

model.wait_for_plots()
