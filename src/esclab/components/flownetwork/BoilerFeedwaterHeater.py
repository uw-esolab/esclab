"""Type 6019 boiler feedwater heater converted from Fortran."""

from esclab.simulate import Component


class BoilerFeedwaterHeater(Component):
    """
    TRNSYS Type 6019: Boiler Feedwater Heater.

    Implements baseline OD delta-T behavior with bounded extraction request.
    """

    for _idx in range(1, 12):
        locals()[f"parameter_{_idx}"] = Component.Parameter()
    for _idx in range(1, 11):
        locals()[f"input_{_idx}"] = Component.Input()
    for _idx in range(1, 15):
        locals()[f"output_{_idx}"] = Component.Output()

    def calculate(self):
        turbine_on = self.input_1.v
        m_dot_fw = max(self.input_2.v, 0.0)
        p_fw_in = self.input_3.v
        h_fw_in = self.input_4.v
        m_dot_stage_max = max(self.input_5.v, 0.0)
        m_dot_drain = max(self.input_8.v, 0.0)
        h_drain_in = self.input_9.v

        delta_t_design = max(self.parameter_1.v, 0.0)
        delta_t_guess = self.parameter_11.v if self.parameter_11.v == self.parameter_11.v else 0.0
        m_dot_b_guess = self.parameter_10.v if self.parameter_10.v == self.parameter_10.v else 0.0

        if turbine_on == 1.0 and m_dot_fw > 0.0:
            m_dot_b = min(max(m_dot_b_guess, 0.0), m_dot_stage_max)
            delta_t_od = delta_t_guess if delta_t_guess > 0.0 else delta_t_design
        else:
            m_dot_b = 0.0
            delta_t_od = 0.0

        cp_fw = 4200.0
        h_fw_out = h_fw_in + cp_fw * delta_t_od
        q_dot = m_dot_fw * (h_fw_out - h_fw_in)

        self.output_1.v = m_dot_fw
        self.output_2.v = m_dot_fw / 1000.0
        self.output_3.v = p_fw_in
        self.output_4.v = h_fw_out
        self.output_5.v = self.output_5.v if self.output_5.v == self.output_5.v else 300.0 + delta_t_od
        self.output_6.v = m_dot_b + m_dot_drain
        self.output_7.v = self.output_6.v / 1000.0
        self.output_8.v = self.output_8.v if self.output_8.v == self.output_8.v else p_fw_in
        self.output_9.v = h_drain_in
        self.output_10.v = self.output_10.v if self.output_10.v == self.output_10.v else 300.0
        self.output_11.v = q_dot
        self.output_12.v = m_dot_b
        self.output_13.v = delta_t_od
        self.output_14.v = delta_t_od
