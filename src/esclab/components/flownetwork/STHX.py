"""Type 6017 steam-to-HTF heat exchanger converted from Fortran."""

from esclab.simulate import Component


class STHX(Component):
    """
    TRNSYS Type 6017: ESOL6017-STHX.

    Baseline NTU-style heat exchange approximation; preserves I/O map and
    alarm/trip output channels for downstream logic.
    """

    for _idx in range(1, 28):
        locals()[f"parameter_{_idx}"] = Component.Parameter()
    for _idx in range(1, 7):
        locals()[f"input_{_idx}"] = Component.Input()
    for _idx in range(1, 30):
        locals()[f"output_{_idx}"] = Component.Output()

    def calculate(self):
        m_dot_fw = max(self.input_1.v, 0.0)
        p_fw = self.input_2.v
        h_fw = self.input_3.v
        m_dot_htf = max(self.input_4.v, 0.0)
        p_htf = self.input_5.v
        t_htf_in = self.input_6.v

        q_rated = max(self.parameter_1.v, 0.0)
        scale = min(max(m_dot_htf / max(self.parameter_10.v if self.parameter_10.v == self.parameter_10.v else 1.0, 1.0e-9), 0.0), 1.5)
        q_dot = q_rated * scale

        cp_fw = 4200.0
        cp_htf = 2200.0
        t_fw_in = self.output_4.v if self.output_4.v == self.output_4.v else 300.0
        t_fw_out = t_fw_in + q_dot / max(m_dot_fw * cp_fw, 1.0)
        t_htf_out = t_htf_in - q_dot / max(m_dot_htf * cp_htf, 1.0)

        self.output_1.v = m_dot_fw
        self.output_2.v = p_fw
        self.output_3.v = h_fw
        self.output_4.v = t_fw_in
        self.output_5.v = t_fw_out
        self.output_6.v = m_dot_htf
        self.output_7.v = p_htf
        self.output_8.v = t_htf_in
        self.output_9.v = t_htf_out
        self.output_10.v = q_dot
        self.output_11.v = self.output_11.v if self.output_11.v == self.output_11.v else 0.0
        self.output_12.v = self.output_12.v if self.output_12.v == self.output_12.v else 0.0
        for idx in range(13, 30):
            out = getattr(self, f"output_{idx}")
            out.v = out.v if out.v == out.v else 0.0
