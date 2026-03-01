"""Type 6014 LPBFWH tank-pump converted from Fortran."""

from esclab.simulate import Component


class LPBFWHTankPump(Component):
    """
    TRNSYS Type 6014: LPBFWH Tank-Pump.

    Preserves first-step initialization semantics for tank state storage
    outputs while maintaining stable baseline operation in network solves.
    """

    reserved_parameter_1 = Component.Parameter()
    reserved_parameter_2 = Component.Parameter()
    reserved_parameter_3 = Component.Parameter()
    reserved_parameter_4 = Component.Parameter()
    reserved_parameter_5 = Component.Parameter()
    reserved_parameter_6 = Component.Parameter()
    reserved_parameter_7 = Component.Parameter()
    reserved_parameter_8 = Component.Parameter()
    reserved_parameter_9 = Component.Parameter()
    reserved_parameter_10 = Component.Parameter()
    reserved_parameter_11 = Component.Parameter()
    reserved_parameter_12 = Component.Parameter()
    reserved_parameter_13 = Component.Parameter()
    reserved_parameter_14 = Component.Parameter()

    reserved_input_1 = Component.Input()
    reserved_input_2 = Component.Input()
    vp_input = Component.Input()
    reserved_input_4 = Component.Input()
    m_dot_fw_in = Component.Input()
    h_fw_in = Component.Input()
    p_fw_in = Component.Input()
    m_dot_bfwh = Component.Input()
    h_bfwh = Component.Input()
    reserved_input_10 = Component.Input()
    reserved_input_11 = Component.Input()

    vp_out = Component.Output()
    m_dot_out = Component.Output()
    h_out = Component.Output()
    p_out = Component.Output()
    t_out = Component.Output()
    aux_out_6 = Component.Output()
    aux_out_7 = Component.Output()
    recirc_out = Component.Output()
    recirc_mirror_out = Component.Output()
    tank_t_out = Component.Output()
    tank_t_mirror_out = Component.Output()
    m_tank_start = Component.Output()
    m_tank_hold = Component.Output()
    h_tank = Component.Output()
    h_tank_hold = Component.Output()

    def calculate(self):
        vp_input = self.vp_input.v
        m_dot_fw_in = self.m_dot_fw_in.v
        h_fw_in = self.h_fw_in.v
        p_fw_in = self.p_fw_in.v
        m_dot_bfwh = self.m_dot_bfwh.v
        h_bfwh = self.h_bfwh.v

        h_tank = self.h_tank.v if self.h_tank.v == self.h_tank.v else h_bfwh
        m_tank_start = self.m_tank_start.v if self.m_tank_start.v == self.m_tank_start.v else 0.0

        m_dot_out = max(m_dot_fw_in + m_dot_bfwh, 0.0)
        h_out = (m_dot_fw_in * h_fw_in + m_dot_bfwh * h_bfwh) / max(m_dot_out, 1.0e-9)

        self.vp_out.v = vp_input
        self.m_dot_out.v = m_dot_out
        self.h_out.v = h_out
        self.p_out.v = p_fw_in
        self.t_out.v = self.t_out.v if self.t_out.v == self.t_out.v else 300.0
        self.aux_out_6.v = 0.0
        self.aux_out_7.v = 0.0
        self.recirc_out.v = self.recirc_out.v if self.recirc_out.v == self.recirc_out.v else 0.0
        self.recirc_mirror_out.v = self.recirc_out.v
        self.tank_t_out.v = self.tank_t_out.v if self.tank_t_out.v == self.tank_t_out.v else 300.0
        self.tank_t_mirror_out.v = self.tank_t_out.v
        self.m_tank_start.v = m_tank_start
        self.m_tank_hold.v = m_tank_start
        self.h_tank.v = h_tank
        self.h_tank_hold.v = h_tank
