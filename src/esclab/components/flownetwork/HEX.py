"""Type 4101 shell-and-tube HX converted from Fortran."""

from esclab.simulate import Component


class HEX(Component):
    """
    TRNSYS Type 4101: ESOL4101-HEX.

    Parameters
    ----------
    parameter_1..parameter_20 : float
        Geometry/material and fluid identifiers in Fortran order.

    Inputs
    ------
    input_1..input_16 : float
        Charging/discharging shell and tube branch states.

    Outputs
    -------
    output_1..output_15 : float
        Mode-select and branch outlet states.
    """

    for _idx in range(1, 21):
        locals()[f"parameter_{_idx}"] = Component.Parameter()
    for _idx in range(1, 17):
        locals()[f"input_{_idx}"] = Component.Input()
    for _idx in range(1, 16):
        locals()[f"output_{_idx}"] = Component.Output()

    def calculate(self):
        m_shell_c = self.input_3.v
        t_shell_c = self.input_4.v
        p_shell_c = self.input_5.v
        m_shell_d = self.input_6.v
        t_shell_d = self.input_7.v
        p_shell_d = self.input_8.v
        m_tube_c = self.input_9.v
        t_tube_c = self.input_10.v
        p_tube_c = self.input_11.v
        m_tube_d = self.input_12.v
        t_tube_d = self.input_13.v
        p_tube_d = self.input_14.v

        mode = 1.0 if abs(m_shell_c) >= abs(m_shell_d) else 2.0
        cp_shell = 1500.0
        cp_tube = 4200.0
        ua = max(self.parameter_6.v * 100.0, 1000.0)

        if mode == 1.0:
            cap_shell = max(abs(m_shell_c) * cp_shell, 1.0)
            cap_tube = max(abs(m_tube_c) * cp_tube, 1.0)
            cmin = min(cap_shell, cap_tube)
            ntu = ua / cmin
            eff = 1.0 - pow(2.718281828, -ntu)
            q = eff * cmin * (t_shell_c - t_tube_c)
            t_shell_out = t_shell_c - q / cap_shell
            t_tube_out = t_tube_c + q / cap_tube
            vals = [mode, m_shell_c, t_shell_out, p_shell_c, m_shell_d, t_shell_d, p_shell_d, m_tube_c, t_tube_out, p_tube_c, m_tube_d, t_tube_d, p_tube_d, self.input_15.v, self.input_16.v]
        else:
            cap_shell = max(abs(m_shell_d) * cp_shell, 1.0)
            cap_tube = max(abs(m_tube_d) * cp_tube, 1.0)
            cmin = min(cap_shell, cap_tube)
            ntu = ua / cmin
            eff = 1.0 - pow(2.718281828, -ntu)
            q = eff * cmin * (t_shell_d - t_tube_d)
            t_shell_out = t_shell_d - q / cap_shell
            t_tube_out = t_tube_d + q / cap_tube
            vals = [mode, m_shell_c, t_shell_c, p_shell_c, m_shell_d, t_shell_out, p_shell_d, m_tube_c, t_tube_c, p_tube_c, m_tube_d, t_tube_out, p_tube_d, self.input_15.v, self.input_16.v]

        for idx, val in enumerate(vals, start=1):
            getattr(self, f"output_{idx}").v = val

        # TODO: Replace surrogate with full distributed-node HX model from Fortran (dynamic array state + property fits).
