"""Type 4035 pipe converted from Fortran."""

import math

from esclab.components.flownetwork.SimplePipe import FricFactor_IC
from esclab.simulate import Component


class Pipe(Component):
    """
    TRNSYS Type 4035: ESOL4035-Pipe.

    Parameters
    ----------
    diameter, l_tot, n_nodes, fluid_id, roughness, n_large_elbows, n_medium_elbows,
    n_standard_elbows, n_contractions, n_expansions, n_gate_valves, init_temp, mc_mult, heat_loss : float

    Inputs
    ------
    temperature, mass_flow, t_amb, pressure, wind, mass_counter : float

    Outputs
    -------
    temperature_out, pressure_out, mass_flow_out, mass_counter_out : float
    """

    diameter = Component.Parameter()
    l_tot = Component.Parameter()
    n_nodes = Component.Parameter()
    fluid_id = Component.Parameter()
    roughness = Component.Parameter()
    n_large_elbows = Component.Parameter()
    n_medium_elbows = Component.Parameter()
    n_standard_elbows = Component.Parameter()
    n_contractions = Component.Parameter()
    n_expansions = Component.Parameter()
    n_gate_valves = Component.Parameter()
    init_temp = Component.Parameter()
    mc_mult = Component.Parameter()
    heat_loss = Component.Parameter()

    temperature = Component.Input()
    mass_flow = Component.Input()
    t_amb = Component.Input()
    pressure = Component.Input()
    wind = Component.Input()
    mass_counter = Component.Input()

    temperature_out = Component.Output()
    pressure_out = Component.Output()
    mass_flow_out = Component.Output()
    mass_counter_out = Component.Output()

    _ff_guess = 0.1

    def calculate(self):
        rho = 1000.0
        mu = 1.0e-3
        area = math.pi * max(self.diameter.v, 1.0e-6) ** 2 / 4.0
        vel = self.mass_flow.v / max(rho * area, 1.0e-9)
        re = abs(rho * vel * max(self.diameter.v, 1.0e-6) / mu)
        ff = FricFactor_IC(max(self.roughness.v, 1.0e-8) / max(self.diameter.v, 1.0e-6), max(re, 1.0), self._ff_guess)
        self._ff_guess = max(ff or 0.1, 0.05)

        k_t = (8.0 * self._ff_guess * max(self.l_tot.v, 0.0)) / (math.pi**2 * max(self.diameter.v, 1.0e-6) ** 5 * rho)
        d_p = k_t * self.mass_flow.v**2

        cp = 2200.0
        q_loss = max(self.heat_loss.v, 0.0) * max(self.l_tot.v, 0.0) * max(self.temperature.v - self.t_amb.v, 0.0)
        t_out = self.temperature.v - q_loss / max(abs(self.mass_flow.v) * cp, 1.0)

        self.mass_flow_out.v = self.mass_flow.v
        self.pressure_out.v = self.pressure.v - d_p
        self.temperature_out.v = t_out
        self.mass_counter_out.v = self.mass_counter.v * self.mc_mult.v if self.mc_mult.v == self.mc_mult.v else self.mass_counter.v