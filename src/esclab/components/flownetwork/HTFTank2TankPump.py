"""Type 6032 HTF tank-to-tank pump converted from Fortran."""

from esclab.simulate import Component


class HTFTank2TankPump(Component):
    """
    TRNSYS Type 6032: HTFTank2TankPump.
    """

    l_down = Component.Parameter()
    pcoef_a = Component.Parameter()
    pcoef_b = Component.Parameter()
    pcoef_c = Component.Parameter()
    reserved_parameter_5 = Component.Parameter()
    reserved_parameter_6 = Component.Parameter()
    reserved_parameter_7 = Component.Parameter()
    reserved_parameter_8 = Component.Parameter()
    reserved_parameter_9 = Component.Parameter()
    reserved_parameter_10 = Component.Parameter()

    pump_on = Component.Input()
    pump_speed_i = Component.Input()
    h_in = Component.Input()
    p_tank1 = Component.Input()
    l_tank1 = Component.Input()
    reserved_input_6 = Component.Input()
    reserved_input_7 = Component.Input()

    m_dot_out = Component.Output()
    vol_dot_out = Component.Output()
    p_out = Component.Output()
    h_out = Component.Output()
    aux_out_1 = Component.Output()
    aux_out_2 = Component.Output()
    m_dot_recycle = Component.Output()
    head_prev = Component.Output()
    vol_dot_pump = Component.Output()
    head_out = Component.Output()
    speed_out = Component.Output()

    def calculate(self):
        pump_on = self.pump_on.v
        pump_speed_i = self.pump_speed_i.v
        p_tank1 = self.p_tank1.v
        l_tank1 = self.l_tank1.v

        rho = 1000.0
        g = 9.81
        pcoef_a = self.pcoef_a.v
        pcoef_b = self.pcoef_b.v
        pcoef_c = self.pcoef_c.v
        l_down = self.l_down.v

        if pump_on == 1.0:
            q = max(pump_speed_i, 1.0e-6)
            m_dot = q * rho
            p_in = p_tank1 + (l_tank1 + l_down) * rho * g
            head = pcoef_a * q * q + pcoef_b * q + pcoef_c * pump_speed_i * pump_speed_i
            p_out = p_in + max(head, 0.0) * rho * g
        else:
            m_dot = 0.0
            p_out = p_tank1
            head = 0.0

        self.m_dot_out.v = max(m_dot, 1.0e-10)
        self.vol_dot_out.v = m_dot / rho
        self.p_out.v = p_out
        self.h_out.v = self.h_in.v
        self.aux_out_1.v = 0.0
        self.aux_out_2.v = 0.0
        self.m_dot_recycle.v = self.vol_dot_pump.v if self.vol_dot_pump.v == self.vol_dot_pump.v else 0.0
        self.head_prev.v = self.head_out.v if self.head_out.v == self.head_out.v else 0.0
        self.vol_dot_pump.v = m_dot / rho
        self.head_out.v = head
        self.speed_out.v = pump_speed_i
