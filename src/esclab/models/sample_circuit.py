"""
RLC circuit example with a voltage source that steps from 0 to 1 at t=0. 
The current is determined by a root-solving method within the voltage source 
component using successive bracketed solver steps. 

Some care is needed to obtain convergence for this example. 

This example uses component classes that are defined locally rather than the 
circuit_elements library. This allows illustration of the differences between 
the bracketed solver approach and the matrix inversion method, which is shown in
the remaining sample circuit examples.
"""

from esclab.simulate import *

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


class VoltageSourceRLC(CircuitElement):
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

class VoltageSource(CircuitElement):
    """Voltage source element."""
    V = Component.Parameter()  # voltage

    def __init__(self):
        super().__init__()

    def presim_setup(self, **kwargs):
        # Topology-agnostic root solving method on u_in = 0.
        self.i_cmd = 0.0
        self.i_prev_step = 0.0
        self.i_obs_prev = None
        self.r_obs_prev = None
        self.i_lo = None
        self.r_lo = None
        self.i_hi = None
        self.r_hi = None
        self.probe_step = None
        self.best_i = 0.0
        self.best_abs_r = float('inf')

    def calculate(self):
        self.u_out = self.V

        # Observed point from the current network state.
        i_obs = float(self.i_in)
        r_obs = float(self.u_in)  # target is 0 at source return node

        if self.model.is_first_iteration:
            # Warm-start each timestep from the last converged current.
            self.i_cmd = float(self.i_prev_step)
            self.i_obs_prev = None
            self.r_obs_prev = None
            self.i_lo = None
            self.r_lo = None
            self.i_hi = None
            self.r_hi = None
            self.probe_step = max(0.01 * abs(self.i_cmd), 1e-4)
            self.best_i = i_obs
            self.best_abs_r = abs(r_obs)
        else:
            abs_r = abs(r_obs)
            if abs_r < self.best_abs_r:
                self.best_abs_r = abs_r
                self.best_i = i_obs

            # Update bracket from sign changes in residual.
            if r_obs <= 0.0:
                if self.r_lo is None or abs(r_obs) < abs(self.r_lo):
                    self.i_lo = i_obs
                    self.r_lo = r_obs
            if r_obs >= 0.0:
                if self.r_hi is None or abs(r_obs) < abs(self.r_hi):
                    self.i_hi = i_obs
                    self.r_hi = r_obs

            i_next = self.i_cmd

            # If bracket is valid, use a regula-falsi step.
            has_bracket = (
                self.i_lo is not None
                and self.i_hi is not None
                and self.r_lo is not None
                and self.r_hi is not None
                and (self.r_lo < 0.0 < self.r_hi)
                and abs(self.i_hi - self.i_lo) > 1e-14
            )

            if has_bracket:
                denom = self.r_hi - self.r_lo
                if abs(denom) > 1e-14:
                    i_rf = (self.i_lo * self.r_hi - self.i_hi * self.r_lo) / denom
                else:
                    i_rf = 0.5 * (self.i_lo + self.i_hi)
                i_mid = 0.5 * (self.i_lo + self.i_hi)
                f = 0.7
                i_next = f * i_rf + (1 - f) * i_mid
                i_next = float(np.clip(i_next, min(self.i_lo, self.i_hi), max(self.i_lo, self.i_hi)))
            else:
                # No bracket yet: estimate local slope, otherwise probe with increasing radius.
                used_slope = False
                if self.i_obs_prev is not None and self.r_obs_prev is not None:
                    di = i_obs - self.i_obs_prev
                    dr = r_obs - self.r_obs_prev
                    if abs(di) > 1e-14 and abs(dr) > 1e-14:
                        slope = dr / di
                        if np.isfinite(slope) and abs(slope) > 1e-10:
                            newton_step = -r_obs / slope
                            max_step = max(self.probe_step, 0.2 * max(abs(i_obs), 1e-6))
                            i_next = i_obs + float(np.clip(newton_step, -max_step, max_step))
                            used_slope = True

                if not used_slope:
                    # Typical passive networks have dr/di < 0 near the solution.
                    direction = 1.0 if r_obs > 0.0 else -1.0
                    i_next = i_obs + direction * self.probe_step
                    self.probe_step = min(self.probe_step * 1.5, max(1.0, 2.0 * abs(i_obs) + 1e-3))

            # Deadband to avoid chattering around the root.
            if abs(r_obs) < 1e-12:
                i_next = i_obs

            # Command relaxation to handle one-iteration lag in loop response.
            f = 0.5
            self.i_cmd = f * float(self.i_cmd) + (1 - f) * float(i_next)

        self.i_out = self.i_cmd
        self.i_obs_prev = i_obs
        self.r_obs_prev = r_obs

        if self.model.is_converged:
            self.i_prev_step = float(self.i_in)
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
# model.vs = VoltageSourceRLC()
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
model.settings.max_iterations = 200
model.settings.tol_rel_global = 1e-8
model.settings.tol_abs_global = 1e-6  
model.settings.learn_rate = 1.  # rely on learning rates in the component calculations for stability

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
                  update_every=4, 
                  nmax_points=1000)

model.add_network_graph()

while model.time < model.settings.stop_time:
    model.step()

model.wait_for_plots()
