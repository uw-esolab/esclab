"""Type 4102 HEX display mapper converted from Fortran."""

from esclab.simulate import Component


class HEXDisplay(Component):
    """
    TRNSYS Type 4102: ESOL4102-HEX-Display.

    Parameters
    ----------
    None.

    Inputs
    ------
    input_1..input_25 : float
        TES mode plus charging/discharging shell/tube in/out streams.

    Outputs
    -------
    output_1..output_16 : float
        Top/bottom display streams and finite-difference temperature derivatives.
    """

    tes_mode = Component.Input()

    charge_shell_m_in = Component.Input()
    charge_shell_t_in = Component.Input()
    charge_shell_p_in = Component.Input()
    discharge_shell_m_out = Component.Input()
    discharge_shell_t_out = Component.Input()
    discharge_shell_p_out = Component.Input()
    charge_tube_m_in = Component.Input()
    charge_tube_t_in = Component.Input()
    charge_tube_p_in = Component.Input()
    discharge_tube_m_in = Component.Input()
    discharge_tube_t_in = Component.Input()
    discharge_tube_p_in = Component.Input()
    charge_shell_m_out = Component.Input()
    charge_shell_t_out = Component.Input()
    charge_shell_p_out = Component.Input()
    discharge_shell_m_in = Component.Input()
    discharge_shell_t_in = Component.Input()
    discharge_shell_p_in = Component.Input()
    charge_tube_m_out = Component.Input()
    charge_tube_t_out = Component.Input()
    charge_tube_p_out = Component.Input()
    discharge_tube_m_out = Component.Input()
    discharge_tube_t_out = Component.Input()
    discharge_tube_p_out = Component.Input()

    display_1 = Component.Output()
    tfd_1 = Component.Output()
    display_3 = Component.Output()
    display_4 = Component.Output()
    tfd_2 = Component.Output()
    display_6 = Component.Output()
    display_7 = Component.Output()
    tfd_3 = Component.Output()
    display_9 = Component.Output()
    display_10 = Component.Output()
    tfd_4 = Component.Output()
    display_12 = Component.Output()
    dtfd_1 = Component.Output()
    dtfd_2 = Component.Output()
    dtfd_3 = Component.Output()
    dtfd_4 = Component.Output()

    _prev = [0.0, 0.0, 0.0, 0.0]

    def calculate(self):
        mode = 1.0 if self.tes_mode.v < 2.0 else 2.0
        if mode == 1.0:
            vals = [
                self.charge_shell_m_in.v, self.charge_shell_t_in.v, self.charge_shell_p_in.v,
                self.charge_tube_m_out.v, self.charge_tube_t_out.v, self.charge_tube_p_out.v,
                self.charge_shell_m_out.v, self.charge_shell_t_out.v, self.charge_shell_p_out.v,
                self.charge_tube_m_in.v, self.charge_tube_t_in.v, self.charge_tube_p_in.v,
            ]
        else:
            vals = [
                self.discharge_shell_m_in.v, self.discharge_shell_t_in.v, self.discharge_shell_p_in.v,
                self.discharge_tube_m_in.v, self.discharge_tube_t_in.v, self.discharge_tube_p_in.v,
                self.discharge_shell_m_out.v, self.discharge_shell_t_out.v, self.discharge_shell_p_out.v,
                self.discharge_tube_m_out.v, self.discharge_tube_t_out.v, self.discharge_tube_p_out.v,
            ]

        self.display_1.v = vals[0]
        self.tfd_1.v = vals[1]
        self.display_3.v = vals[2]
        self.display_4.v = vals[3]
        self.tfd_2.v = vals[4]
        self.display_6.v = vals[5]
        self.display_7.v = vals[6]
        self.tfd_3.v = vals[7]
        self.display_9.v = vals[8]
        self.display_10.v = vals[9]
        self.tfd_4.v = vals[10]
        self.display_12.v = vals[11]

        ts = max(self.model.settings.timestep, 1.0)
        curr = [self.tfd_1.v, self.tfd_2.v, self.tfd_3.v, self.tfd_4.v]
        self.dtfd_1.v = (curr[0] - self._prev[0]) / ts
        self.dtfd_2.v = (curr[1] - self._prev[1]) / ts
        self.dtfd_3.v = (curr[2] - self._prev[2]) / ts
        self.dtfd_4.v = (curr[3] - self._prev[3]) / ts
        self._prev = curr
