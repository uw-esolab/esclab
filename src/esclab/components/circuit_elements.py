"""
Circuit component library using simple models. 
The components are set up to be compatible with the matrix inversion method.
"""

from esclab.simulate import *

class __CircuitElement(Component):
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

class VoltageSource(__CircuitElement):
    """Ideal voltage source element."""
    V = Component.Parameter()  # source voltage
    omega_hi = Component.Parameter(0)  # cutoff frequency for high-pass behavior
    omega_lo = Component.Parameter(0)  # cutoff frequency for low-pass behavior

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


class Resistor(__CircuitElement):
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

class Capacitor(__CircuitElement):
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

class Inductor(__CircuitElement):
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
