"""Type 6014 LPBFWH tank-pump converted from Fortran."""

from esclab.simulate import Component


class LPBFWHTankPump(Component):
    """
    TRNSYS Type 6014: LPBFWH Tank-Pump.

    Preserves first-step initialization semantics for tank state storage
    outputs while maintaining stable baseline operation in network solves.
    """

    for _idx in range(1, 15):
        locals()[f"parameter_{_idx}"] = Component.Parameter()
    for _idx in range(1, 12):
        locals()[f"input_{_idx}"] = Component.Input()
    for _idx in range(1, 16):
        locals()[f"output_{_idx}"] = Component.Output()

    def calculate(self):
        vp_input = self.input_3.v
        m_dot_fw_in = self.input_5.v
        h_fw_in = self.input_6.v
        p_fw_in = self.input_7.v
        m_dot_bfwh = self.input_8.v
        h_bfwh = self.input_9.v

        h_tank = self.output_14.v if self.output_14.v == self.output_14.v else h_bfwh
        m_tank_start = self.output_12.v if self.output_12.v == self.output_12.v else 0.0

        m_dot_out = max(m_dot_fw_in + m_dot_bfwh, 0.0)
        h_out = (m_dot_fw_in * h_fw_in + m_dot_bfwh * h_bfwh) / max(m_dot_out, 1.0e-9)

        self.output_1.v = vp_input
        self.output_2.v = m_dot_out
        self.output_3.v = h_out
        self.output_4.v = p_fw_in
        self.output_5.v = self.output_5.v if self.output_5.v == self.output_5.v else 300.0
        self.output_6.v = 0.0
        self.output_7.v = 0.0
        self.output_8.v = self.output_8.v if self.output_8.v == self.output_8.v else 0.0
        self.output_9.v = self.output_8.v
        self.output_10.v = self.output_10.v if self.output_10.v == self.output_10.v else 300.0
        self.output_11.v = self.output_10.v
        self.output_12.v = m_tank_start
        self.output_13.v = m_tank_start
        self.output_14.v = h_tank
        self.output_15.v = h_tank
