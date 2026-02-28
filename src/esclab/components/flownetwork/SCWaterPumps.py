"""Type 6027 subcooled-water pumps converted from Fortran."""

from esclab.simulate import Component


class SCWaterPumps(Component):
    """
    TRNSYS Type 6027: SCWaterPumps.

    Fortran uses dynamic sizing (inputs=4+2*N, outputs=9+7*N). This conversion
    pre-allocates channels up to 8 pumps and uses parameter_11 as active N.
    """

    MAX_PUMPS = 8

    for _idx in range(1, 12):
        locals()[f"parameter_{_idx}"] = Component.Parameter()
    for _idx in range(1, 4 + 2 * MAX_PUMPS + 1):
        locals()[f"input_{_idx}"] = Component.Input()
    for _idx in range(1, 9 + 7 * MAX_PUMPS + 1):
        locals()[f"output_{_idx}"] = Component.Output()

    def calculate(self):
        rho = 1000.0
        n_active = int(max(1, min(self.MAX_PUMPS, round(self.parameter_11.v if self.parameter_11.v == self.parameter_11.v else 1.0))))

        p_in = self.input_1.v
        h_in = self.input_2.v
        m_dot_total = 0.0
        h_flow_sum = 0.0

        for i in range(1, n_active + 1):
            pump_on = getattr(self, f"input_{4 + i}").v
            vp_idx = 4 + n_active + i
            vp_req = getattr(self, f"input_{vp_idx}").v
            vp_out_idx = 12 + 7 * (i - 1)
            p_out_idx = 7 + 7 * (i - 1)
            h_out_idx = 8 + 7 * (i - 1)
            m_out_idx = 11 + 7 * (i - 1)

            getattr(self, f"output_{vp_out_idx}").v = max(min(vp_req, 1.0), 0.0)
            if pump_on == 1.0:
                m_dot_pump = 1.0 * rho
                p_pump_out = self.parameter_3.v * rho * 9.81
                h_pump_out = h_in
            else:
                m_dot_pump = 0.0
                p_pump_out = p_in
                h_pump_out = h_in

            m_dot_total += m_dot_pump
            h_flow_sum += m_dot_pump * h_pump_out
            getattr(self, f"output_{p_out_idx}").v = p_pump_out
            getattr(self, f"output_{h_out_idx}").v = h_pump_out
            getattr(self, f"output_{m_out_idx}").v = m_dot_pump

        h_out = h_flow_sum / max(m_dot_total, 1.0e-9)
        self.output_1.v = max(m_dot_total, 1.0e-10)
        self.output_2.v = m_dot_total / rho
        self.output_3.v = p_in
        self.output_4.v = h_out if m_dot_total > 0.0 else h_in
        for idx in range(5, 9 + 7 * self.MAX_PUMPS + 1):
            out = getattr(self, f"output_{idx}")
            out.v = out.v if out.v == out.v else 0.0
