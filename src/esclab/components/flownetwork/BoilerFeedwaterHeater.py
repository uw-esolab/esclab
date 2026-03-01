"""Type 6019 boiler feedwater heater converted from Fortran."""

from esclab.simulate import Component


class BoilerFeedwaterHeater(Component):
    """
    TRNSYS Type 6019: Boiler Feedwater Heater.

    Implements baseline OD delta-T behavior with bounded extraction request.
    """

    delta_t_design = Component.Parameter()
    reserved_parameter_2 = Component.Parameter()
    reserved_parameter_3 = Component.Parameter()
    reserved_parameter_4 = Component.Parameter()
    reserved_parameter_5 = Component.Parameter()
    reserved_parameter_6 = Component.Parameter()
    reserved_parameter_7 = Component.Parameter()
    reserved_parameter_8 = Component.Parameter()
    reserved_parameter_9 = Component.Parameter()
    m_dot_b_guess = Component.Parameter()
    delta_t_guess = Component.Parameter()

    turbine_on = Component.Input()
    m_dot_fw = Component.Input()
    p_fw_in = Component.Input()
    h_fw_in = Component.Input()
    m_dot_stage_max = Component.Input()
    reserved_input_6 = Component.Input()
    reserved_input_7 = Component.Input()
    m_dot_drain = Component.Input()
    h_drain_in = Component.Input()
    reserved_input_10 = Component.Input()

    m_dot_fw_out = Component.Output()
    vol_dot_fw_out = Component.Output()
    p_fw_out = Component.Output()
    h_fw_out = Component.Output()
    t_fw_out = Component.Output()
    m_dot_drain_out = Component.Output()
    vol_dot_drain_out = Component.Output()
    p_drain_out = Component.Output()
    h_drain_out = Component.Output()
    t_drain_out = Component.Output()
    q_dot = Component.Output()
    m_dot_bleed = Component.Output()
    delta_t_od_out = Component.Output()
    delta_t_display = Component.Output()

    def calculate(self):
        turbine_on = self.turbine_on.v
        m_dot_fw = max(self.m_dot_fw.v, 0.0)
        p_fw_in = self.p_fw_in.v
        h_fw_in = self.h_fw_in.v
        m_dot_stage_max = max(self.m_dot_stage_max.v, 0.0)
        m_dot_drain = max(self.m_dot_drain.v, 0.0)
        h_drain_in = self.h_drain_in.v

        delta_t_design = max(self.delta_t_design.v, 0.0)
        delta_t_guess = self.delta_t_guess.v if self.delta_t_guess.v == self.delta_t_guess.v else 0.0
        m_dot_b_guess = self.m_dot_b_guess.v if self.m_dot_b_guess.v == self.m_dot_b_guess.v else 0.0

        if turbine_on == 1.0 and m_dot_fw > 0.0:
            m_dot_b = min(max(m_dot_b_guess, 0.0), m_dot_stage_max)
            delta_t_od = delta_t_guess if delta_t_guess > 0.0 else delta_t_design
        else:
            m_dot_b = 0.0
            delta_t_od = 0.0

        cp_fw = 4200.0
        h_fw_out = h_fw_in + cp_fw * delta_t_od
        q_dot = m_dot_fw * (h_fw_out - h_fw_in)

        self.m_dot_fw_out.v = m_dot_fw
        self.vol_dot_fw_out.v = m_dot_fw / 1000.0
        self.p_fw_out.v = p_fw_in
        self.h_fw_out.v = h_fw_out
        self.t_fw_out.v = self.t_fw_out.v if self.t_fw_out.v == self.t_fw_out.v else 300.0 + delta_t_od
        self.m_dot_drain_out.v = m_dot_b + m_dot_drain
        self.vol_dot_drain_out.v = self.m_dot_drain_out.v / 1000.0
        self.p_drain_out.v = self.p_drain_out.v if self.p_drain_out.v == self.p_drain_out.v else p_fw_in
        self.h_drain_out.v = h_drain_in
        self.t_drain_out.v = self.t_drain_out.v if self.t_drain_out.v == self.t_drain_out.v else 300.0
        self.q_dot.v = q_dot
        self.m_dot_bleed.v = m_dot_b
        self.delta_t_od_out.v = delta_t_od
        self.delta_t_display.v = delta_t_od
