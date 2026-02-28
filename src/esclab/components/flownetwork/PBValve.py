"""Type 6001 power-block valve converted from Fortran."""

from esclab.simulate import Component


class PBValve(Component):
    """
    TRNSYS Type 6001: ESOL6001-PB_Valve.

    Parameters
    ----------
    valve_diameter, valve_speed, valve_type : float

    Inputs
    ------
    m_dot_in, p_in, h_in, vp_input : float

    Outputs
    -------
    m_dot_out, vol_dot_out, p_out, h_out, t_out, delta_p, vp_output : float
    """

    valve_diameter = Component.Parameter()
    valve_speed = Component.Parameter()
    valve_type = Component.Parameter()

    m_dot_in = Component.Input()
    p_in = Component.Input()
    h_in = Component.Input()
    vp_input = Component.Input()

    m_dot_out = Component.Output()
    vol_dot_out = Component.Output()
    p_out = Component.Output()
    h_out = Component.Output()
    t_out = Component.Output()
    delta_p = Component.Output()
    vp_output = Component.Output()

    def _cv(self, vp):
        # TODO: Replace with PB_CV_data map from ESOL6015_myfunctions.
        return max(1.0e-4, (self.valve_diameter.v * 39.3701) ** 2 * max(vp, 1.0e-4) * (20.0 + self.valve_type.v))

    def calculate(self):
        vp_req = max(min(self.vp_input.v, 1.0), 0.0)
        if self.model.is_first_step:
            vp = vp_req
        elif self.model.is_first_iteration:
            vp_prev = self.vp_output.v if self.vp_output.v == self.vp_output.v else vp_req
            ts = self.model.settings.timestep
            step = self.valve_speed.v * ts / 90.0
            if vp_req > vp_prev:
                vp = min(vp_prev + step, vp_req)
            else:
                vp = max(vp_prev - step, vp_req)
        else:
            vp = self.vp_output.v if self.vp_output.v == self.vp_output.v else vp_req

        rho = 1000.0
        q_gpm = (self.m_dot_in.v / rho) * 15850.3
        cv = self._cv(vp)
        delta_p_psi = (q_gpm / max(cv, 1.0e-6)) ** 2
        delta_p = delta_p_psi * 6894.76

        self.m_dot_out.v = self.m_dot_in.v
        self.vol_dot_out.v = self.m_dot_in.v / rho
        self.p_out.v = self.p_in.v - delta_p
        self.h_out.v = self.h_in.v
        self.t_out.v = self.t_out.v if self.t_out.v == self.t_out.v else 300.0
        self.delta_p.v = delta_p
        self.vp_output.v = vp