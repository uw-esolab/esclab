"""Type 4034 solar-field sector surrogate converted from Fortran."""

from esclab.simulate import Component


class SolarFieldSector(Component):
    """
    TRNSYS Type 4034: ESOL4034-SolarFieldSector.

    Parameters
    ----------
    parameter_1..parameter_41 : float
        Optical/thermal configuration terms from Fortran input order.

    Inputs
    ------
    mass_flow, pressure, temperature, control_signal, mass_counter, ani, t_amb, t_sky, wind, theta, phi, t_tracking : float
        Sector hydraulic, thermal, and weather conditions.

    Outputs
    -------
    output_1..output_22 : float
        Main outlet state and sector diagnostics in Fortran order.
    """

    for _idx in range(1, 42):
        locals()[f"parameter_{_idx}"] = Component.Parameter()

    mass_flow = Component.Input()
    pressure = Component.Input()
    temperature = Component.Input()
    control_signal = Component.Input()
    mass_counter = Component.Input()
    ani = Component.Input()
    t_amb = Component.Input()
    t_sky = Component.Input()
    wind = Component.Input()
    theta = Component.Input()
    phi = Component.Input()
    t_tracking = Component.Input()

    for _idx in range(1, 23):
        locals()[f"output_{_idx}"] = Component.Output()

    def calculate(self):
        eta_defocus = max(min(self.parameter_1.v, 1.0), 0.0)
        eta_tracking = max(min(self.parameter_8.v, 1.0), 0.0)
        eta_soil = max(min(self.parameter_9.v, 1.0), 0.0)
        eta_reflect = max(min(self.parameter_10.v, 1.0), 0.0)
        sf_avail = max(min(self.parameter_11.v, 1.0), 0.0)
        w_ap = max(self.parameter_14.v, 0.0)
        n_loop = max(self.parameter_18.v, 1.0)
        cp = 2200.0

        eta_tot = eta_defocus * eta_tracking * eta_soil * eta_reflect * sf_avail
        q_dot = max(self.ani.v, 0.0) * eta_tot * w_ap * n_loop
        m_dot = max(self.mass_flow.v, 1.0e-6)
        t_out = self.temperature.v + q_dot / (m_dot * cp)

        self.output_1.v = self.mass_flow.v
        self.output_2.v = self.pressure.v
        self.output_3.v = t_out
        self.output_4.v = t_out
        self.output_5.v = self.mass_counter.v
        self.output_6.v = q_dot
        self.output_7.v = eta_tot
        self.output_8.v = self.theta.v
        self.output_9.v = self.phi.v
        self.output_10.v = self.wind.v
        for idx in range(11, 23):
            getattr(self, f"output_{idx}").v = getattr(self, f"output_{idx}").v if idx > 12 else 0.0

        # TODO: Map neural-network optics and detailed loop/header thermal-hydraulic logic from Fortran modules.