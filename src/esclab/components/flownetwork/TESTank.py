"""Type 4100 TES tank converted from Fortran."""

import math

from esclab.simulate import Component


class TESTank(Component):
    """
    TRNSYS Type 4100: ESOL4100-TESTank.

    Parameters
    ----------
    d_in, height, ins_th, k_iso, emiss, t0, l0, id_fluid : float

    Inputs
    ------
    t_in, m_in, m_out, t_env, v_env, v_air : float

    Outputs
    -------
    output_1..output_9 : float
        Pump draw, state, wall temperature, and heat loss diagnostics.
    """

    d_in = Component.Parameter()
    height = Component.Parameter()
    ins_th = Component.Parameter()
    k_iso = Component.Parameter()
    emiss = Component.Parameter()
    t0 = Component.Parameter()
    l0 = Component.Parameter()
    id_fluid = Component.Parameter()

    t_in = Component.Input()
    m_in = Component.Input()
    m_out = Component.Input()
    t_env = Component.Input()
    v_env = Component.Input()
    v_air = Component.Input()

    for _idx in range(1, 10):
        locals()[f"output_{_idx}"] = Component.Output()

    _tank_temp = 300.0
    _tank_level = 0.0

    def calculate(self):
        area = math.pi * max(self.d_in.v, 1.0e-6) ** 2 / 4.0
        rho = 1800.0
        cp = 1500.0
        ts = max(self.model.settings.timestep, 1.0)

        if self.model.is_first_step:
            self._tank_temp = self.t0.v
            self._tank_level = self.l0.v

        mass = max(area * self._tank_level * rho, 1.0)
        q_loss = max(self.k_iso.v, 0.0) * area * max(self._tank_temp - self.t_env.v, 0.0)
        d_m = (self.m_in.v - self.m_out.v) * ts
        d_u = (self.m_in.v * cp * self.t_in.v - self.m_out.v * cp * self._tank_temp - q_loss) * ts
        mass_new = max(mass + d_m, 1.0)
        self._tank_temp = (mass * cp * self._tank_temp + d_u) / (mass_new * cp)
        self._tank_level = max(min(mass_new / (rho * area), self.height.v), 0.0)
        p_out = 101325.0 + rho * 9.81 * self._tank_level

        self.output_1.v = self.m_out.v
        self.output_2.v = self.m_out.v / rho
        self.output_3.v = self._tank_temp
        self.output_4.v = p_out
        self.output_5.v = self._tank_temp
        self.output_6.v = self._tank_level
        self.output_7.v = mass_new
        self.output_8.v = self._tank_temp - 5.0
        self.output_9.v = q_loss

