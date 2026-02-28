"""Type 6032 HTF tank-to-tank pump converted from Fortran."""

from esclab.simulate import Component


class HTFTank2TankPump(Component):
    """
    TRNSYS Type 6032: HTFTank2TankPump.
    """

    for _idx in range(1, 11):
        locals()[f"parameter_{_idx}"] = Component.Parameter()
    for _idx in range(1, 8):
        locals()[f"input_{_idx}"] = Component.Input()
    for _idx in range(1, 12):
        locals()[f"output_{_idx}"] = Component.Output()

    def calculate(self):
        pump_on = self.input_1.v
        pump_speed_i = self.input_2.v
        p_tank1 = self.input_4.v
        l_tank1 = self.input_5.v

        rho = 1000.0
        g = 9.81
        pcoef_a = self.parameter_2.v
        pcoef_b = self.parameter_3.v
        pcoef_c = self.parameter_4.v
        l_down = self.parameter_1.v

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

        self.output_1.v = max(m_dot, 1.0e-10)
        self.output_2.v = m_dot / rho
        self.output_3.v = p_out
        self.output_4.v = self.input_3.v
        self.output_5.v = 0.0
        self.output_6.v = 0.0
        self.output_7.v = self.output_9.v if self.output_9.v == self.output_9.v else 0.0
        self.output_8.v = self.output_10.v if self.output_10.v == self.output_10.v else 0.0
        self.output_9.v = m_dot / rho
        self.output_10.v = head
        self.output_11.v = pump_speed_i
