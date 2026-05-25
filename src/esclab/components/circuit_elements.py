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

    def calculate(self):
        if self.coupled_eqs is not None:
            # All series elements pass current straight through: i_out = i_in.
            # Child classes call super().calculate() to contribute this equation,
            # then add their own voltage equation(s).
            self.coupled_eqs.add_equation({
                self.i_out: 1.0,
                self.coupled_eqs.source(self.i_in): -1.0,
            }, rhs=0.0)

class VoltageSource(__CircuitElement):
    """Ideal voltage source element."""
    V = Component.Parameter()  # source voltage
    omega_hi = Component.Parameter(0)  # cutoff frequency for high-pass behavior
    omega_lo = Component.Parameter(0)  # cutoff frequency for low-pass behavior

    def __init__(self):
        super().__init__()
        self.V_t = self.V  # stored voltage for time-varying sources

    def calculate(self):
        super().calculate()  # i_out = i_in pass-through (only when coupled_eqs is set)
        if self.coupled_eqs is not None:
            V_t = 1
            V_t += np.sin(2 * np.pi * self.omega_hi * self.model.time) * 0.25
            V_t += np.sin(2 * np.pi * self.omega_lo * self.model.time) * 0.25
            V_t += self.V * V_t
            self.coupled_eqs.add_equation({
                self.u_out: 1.0,
                self.coupled_eqs.source(self.u_in): -1.0,
            }, rhs=V_t)
            self.coupled_eqs.add_equation({
                self.coupled_eqs.source(self.u_in): 1.0,
            }, rhs=0.0)


class Resistor(__CircuitElement):
    """Resistor element."""
    R = Component.Parameter()  # resistance

    def calculate(self):
        super().calculate()  # i_out = i_in pass-through
        if self.coupled_eqs is not None:
            # u_out = u_in - R*i_in    =>    u_out - u_in + R*i_in = 0
            self.coupled_eqs.add_equation({
                self.u_out: 1.0,
                self.coupled_eqs.source(self.u_in): -1.0,
                self.coupled_eqs.source(self.i_in): self.R,
            }, rhs=0.0)

class Capacitor(__CircuitElement):
    """Capacitor element."""
    C = Component.Parameter()  # capacitance
    U_C0 = Component.Parameter(0.)  # initial voltage across the capacitor

    def presim_setup(self, **kwargs):
        self.U_C_prev = self.U_C0  # stored capacitor voltage from previous timestep

    def calculate(self):
        super().calculate()  # i_out = i_in pass-through
        dt = self.model.settings.timestep
        if self.coupled_eqs is not None:
            # Backward Euler: u_out = u_in - U_C_prev - (dt/C)*i_in
            # u_out - u_in + (dt/C)*i_in = -U_C_prev
            self.coupled_eqs.add_equation({
                self.u_out: 1.0,
                self.coupled_eqs.source(self.u_in): -1.0,
                self.coupled_eqs.source(self.i_in): dt / self.C,
            }, rhs=-self.U_C_prev)
        elif self.model.is_converged:
            # After convergence, store the capacitor voltage for the next timestep.
            self.U_C_prev = self.U_C_prev + self.i_in * dt / self.C

class Inductor(__CircuitElement):
    """Inductor element."""
    L = Component.Parameter()  # inductance
    I_L0 = Component.Parameter(0.)  # initial current through the inductor

    def presim_setup(self, **kwargs):
        self.I_L_prev = self.I_L0  # stored inductor current from previous timestep

    def calculate(self):
        super().calculate()  # i_out = i_in pass-through
        dt = self.model.settings.timestep
        L = self.L
        if self.coupled_eqs is not None:
            # Backward Euler: u_out = u_in - (L/dt)*i_in + (L/dt)*I_L_prev
            # u_out - u_in + (L/dt)*i_in = (L/dt)*I_L_prev
            self.coupled_eqs.add_equation({
                self.u_out: 1.0,
                self.coupled_eqs.source(self.u_in): -1.0,
                self.coupled_eqs.source(self.i_in): L / dt,
            }, rhs=(L / dt) * self.I_L_prev)
        elif self.model.is_converged:
            # After convergence, store the inductor current for the next timestep.
            self.I_L_prev = self.i_in


class TeeOut(Component):
    """Tee element for splitting voltage and current."""
    u_in = Component.Input()   # voltage
    i_in = Component.Input()   # current

    u_out = Component.Output()  # voltage output
    i_out_1 = Component.Output()  # current output 1
    i_out_2 = Component.Output()  # current output 2

    def calculate(self):
        if self.coupled_eqs is not None:
            # Current splits: i_in = i_out_1 + i_out_2
            self.coupled_eqs.add_equation({
                self.coupled_eqs.source(self.i_in): 1.0,
                self.i_out_1: -1.0,
                self.i_out_2: -1.0,
            }, rhs=0.0)
            # Output voltage equals input voltage: u_out = u_in
            self.coupled_eqs.add_equation({
                self.u_out: 1.0,
                self.coupled_eqs.source(self.u_in): -1.0,
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
        if self.coupled_eqs is not None:
            # Current combines: i_out = i_in_1 + i_in_2
            self.coupled_eqs.add_equation({
                self.i_out: 1.0,
                self.coupled_eqs.source(self.i_in_1): -1.0,
                self.coupled_eqs.source(self.i_in_2): -1.0,
            }, rhs=0.0)
            # Output voltage equals both input voltages: u_out = u_in_1 = u_in_2
            self.coupled_eqs.add_equation({
                self.u_out: 1.0,
                self.coupled_eqs.source(self.u_in_1): -1.0,
            }, rhs=0.0)
            self.coupled_eqs.add_equation({
                self.u_out: 1.0,
                self.coupled_eqs.source(self.u_in_2): -1.0,
            }, rhs=0.0)
