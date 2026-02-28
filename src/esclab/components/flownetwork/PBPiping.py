"""Type 6016 power-block piping converted from Fortran."""

import math

from esclab.components.flownetwork.SimplePipe import FricFactor_IC
from esclab.simulate import Component


class PBPiping(Component):
    """
    TRNSYS Type 6016: ESOL6016-PB_Piping.

    Parameters
    ----------
    pipe_id, pipe_length, roughness, elevation_change : float

    Inputs
    ------
    m_dot_in, p_in, h_in : float

    Outputs
    -------
    m_dot_out, vol_dot_out, p_out, h_out, t_out, delta_p, ff : float
    """

    pipe_id = Component.Parameter()
    pipe_length = Component.Parameter()
    roughness = Component.Parameter()
    elevation_change = Component.Parameter()

    m_dot_in = Component.Input()
    p_in = Component.Input()
    h_in = Component.Input()

    m_dot_out = Component.Output()
    vol_dot_out = Component.Output()
    p_out = Component.Output()
    h_out = Component.Output()
    t_out = Component.Output()
    delta_p = Component.Output()
    ff = Component.Output()

    def calculate(self):
        rho = 1000.0
        mu = 1.0e-3
        d = max(self.pipe_id.v, 1.0e-6)
        area = math.pi * d**2 / 4.0
        vel = self.m_dot_in.v / max(rho * area, 1.0e-9)
        re = abs(rho * vel * d / mu)
        ff_guess = max(self.ff.v if self.ff.v == self.ff.v else 0.1, 0.1)
        ff = FricFactor_IC(max(self.roughness.v, 1.0e-8) / d, max(re, 1.0), ff_guess)
        ff = max(ff or ff_guess, 0.05)

        k_t = (8.0 * ff * max(self.pipe_length.v, 0.0)) / (math.pi**2 * d**5 * rho)
        delta_p_fric = k_t * self.m_dot_in.v**2
        delta_p_elev = self.elevation_change.v * 9.81 * rho if abs(self.m_dot_in.v) > 0.01 else 0.0
        delta_p = delta_p_fric + delta_p_elev

        self.m_dot_out.v = self.m_dot_in.v
        self.vol_dot_out.v = self.m_dot_in.v / rho
        self.p_out.v = self.p_in.v - delta_p
        self.h_out.v = self.h_in.v
        self.t_out.v = self.t_out.v if self.t_out.v == self.t_out.v else 300.0
        self.delta_p.v = delta_p
        self.ff.v = ff
