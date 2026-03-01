"""Type 4008 tee-return mixer converted from Fortran."""

import numpy as np
from esclab.components.esol_properties import Incompressible as Inc

from esclab.simulate import Component


class TeeReturnSimple(Component):
    """
    TRNSYS Type 4008: ESOL4008-Tee_Return.

    Parameters
    ----------
    diameter_1, diameter_2, diameter_out : float
        Branch and outlet diameters [m].
    fluid_id : str
        Fluid identifier used for specific heat lookups.

    Inputs
    ------
    m_dot_1, t_1, p_1, m_dot_2, t_2, p_2 : float
        Branch mass flow, temperature, and pressure.
    mass_counter_1, mass_counter_2 : float
        Running mass counters.

    Outputs
    -------
    m_dot, pressure, temperature, p_1_in, p_2_in, mass_counter : float
        Mixed outlet state and diagnostics.
    """

    diameter_1 = Component.Parameter()
    diameter_2 = Component.Parameter()
    diameter_out = Component.Parameter()
    fluid_id = Component.Parameter()

    m_dot_1 = Component.Input()
    t_1 = Component.Input()
    p_1 = Component.Input()
    m_dot_2 = Component.Input()
    t_2 = Component.Input()
    p_2 = Component.Input()
    mass_counter_1 = Component.Input()
    mass_counter_2 = Component.Input()

    m_dot = Component.Output()
    pressure = Component.Output()
    temperature = Component.Output()
    p_1_in = Component.Output()
    p_2_in = Component.Output()
    mass_counter = Component.Output()
    _props = Inc()

    def calculate(self):
        fluid_name = str(self.fluid_id.v) if self.fluid_id.v == self.fluid_id.v else "Nitrate Salt"
        m_dot_tot = self.m_dot_1.v + self.m_dot_2.v
        if self.model.is_first_step:
            self.m_dot.v = m_dot_tot
            self.pressure.v = self.p_1.v
            self.temperature.v = 500.0
            self.p_1_in.v = self.p_1.v
            self.p_2_in.v = self.p_2.v
            self.mass_counter.v = self.mass_counter_1.v + self.mass_counter_2.v
            return

        if abs(m_dot_tot) < 1.0e-12:
            t_out = self.temperature.v if np.isfinite(self.temperature.v) else 0.5 * (self.t_1.v + self.t_2.v)
        else:
            cp1 = float(self._props.specheat(fluid_name, max(self.t_1.v, 273.0), self.p_1.v))
            cp2 = float(self._props.specheat(fluid_name, max(self.t_2.v, 273.0), self.p_2.v))
            cp_t = (self.m_dot_1.v * cp1 * self.t_1.v + self.m_dot_2.v * cp2 * self.t_2.v) / m_dot_tot
            t_guess = float(np.clip(self.temperature.v if np.isfinite(self.temperature.v) else 500.0, 273.0, 1000.0))
            for _ in range(50):
                cp_guess = float(self._props.specheat(fluid_name, max(t_guess, 273.0), self.p_1.v))
                err = cp_t - t_guess * cp_guess
                if abs(err) < 100.0:
                    break
                t_guess = float(np.clip(t_guess + np.sign(err), 273.0, 1000.0))
            t_out = t_guess

        if self.p_1.v * self.p_2.v < 0.0:
            p_out = max(self.p_1.v, self.p_2.v)
        elif self.p_1.v > 0.0 and self.p_2.v > 0.0:
            p_out = 0.5 * (self.p_1.v + self.p_2.v)
        else:
            p_out = 1.0

        self.m_dot.v = m_dot_tot
        self.pressure.v = p_out
        self.temperature.v = t_out
        self.p_1_in.v = self.p_1.v
        self.p_2_in.v = self.p_2.v
        self.mass_counter.v = self.mass_counter_1.v + self.mass_counter_2.v
