"""Type 4012 expansion system converted from Fortran."""

import math

from esclab.simulate import Component


class ExpansionSystem(Component):
    """
    TRNSYS Type 4012: ESOL4012-ExpansionSystem.

    Parameters
    ----------
    parameter_1..parameter_21 : float
        Fortran parameter-index mapping (fluid, pressures, geometry, and controls).

    Inputs
    ------
    input_1..input_8 : float
        Valve positions and process state as in Fortran input order.

    Outputs
    -------
    output_1..output_15 : float
        Tank levels, outlet state, branch flow diagnostics, and alarms.
    """

    for _idx in range(1, 22):
        locals()[f"parameter_{_idx}"] = Component.Parameter()
    for _idx in range(1, 9):
        locals()[f"input_{_idx}"] = Component.Input()
    for _idx in range(1, 16):
        locals()[f"output_{_idx}"] = Component.Output()

    _last_level_ev = 0.0
    _last_level_of = 0.0

    def _tank_level(self, mass, diameter, height, n_tanks, rho):
        area = math.pi * (max(diameter, 1.0e-6) / 2.0) ** 2
        denom = max(area * max(height, 1.0e-6) * max(n_tanks, 1.0), 1.0e-9)
        return max(mass / max(rho, 1.0) / denom, 0.0)

    def calculate(self):
        rho = 1000.0
        fluid_id = self.parameter_1.v
        p_ev = self.parameter_2.v
        p_of = self.parameter_3.v
        n_of = self.parameter_8.v
        n_ev = self.parameter_9.v
        h_of = self.parameter_10.v
        h_ev = self.parameter_11.v
        d_of = self.parameter_13.v
        d_ev = self.parameter_14.v
        total_mass = self.parameter_15.v
        only_exp = self.parameter_16.v
        t_of = self.parameter_18.v
        t_ev = self.parameter_19.v

        m_counter = self.input_8.v
        m_dot_in = self.input_5.v
        t_in = self.input_6.v

        if self.model.is_first_step:
            if only_exp == 1.0:
                mass_ev = max(total_mass - m_counter, 0.0)
                level_ev = self._tank_level(mass_ev, d_ev, h_ev, n_ev, rho)
                level_of = 0.0
            else:
                mass_of = math.pi * (max(d_of, 1.0e-6) / 2.0) ** 2 * max(h_of, 0.0) * max(n_of, 0.0) * rho
                mass_ev = max(total_mass - mass_of - m_counter, 0.0)
                level_of = self._tank_level(mass_of, d_of, h_of, n_of, rho)
                level_ev = self._tank_level(mass_ev, d_ev, h_ev, n_ev, rho)

            self._last_level_ev = level_ev
            self._last_level_of = level_of
            self.output_1.v = level_ev
            self.output_2.v = level_of
            self.output_3.v = m_dot_in
            self.output_4.v = t_in
            self.output_5.v = p_ev + rho * 9.81 * level_ev * max(h_ev, 0.0)
            self.output_6.v = 0.0
            self.output_7.v = 0.0
            self.output_8.v = 0.0
            self.output_9.v = 0.0
            self.output_10.v = t_of
            self.output_11.v = t_ev
            for idx in range(12, 16):
                getattr(self, f"output_{idx}").v = 0.0
            return

        # TODO: Map full hydraulic matrix solver from Fortran (SF_piping_functions / CV_data dependencies).
        self.output_1.v = self._last_level_ev
        self.output_2.v = self._last_level_of
        self.output_3.v = m_dot_in
        self.output_4.v = t_in
        self.output_5.v = max(self.output_5.v, p_of)
