"""Type 6031 Solana hydraulic model converted from Fortran."""

from esclab.simulate import Component


class SolanaHydraulicModel(Component):
    """
    TRNSYS Type 6031: SolanaHydraulicModel.

    Uses the Fortran first-step flow-fraction structure and a stable pressure
    ratio fallback for iterative updates.
    """

    reserved_parameter_1 = Component.Parameter()
    reserved_parameter_2 = Component.Parameter()
    reserved_parameter_3 = Component.Parameter()
    reserved_parameter_4 = Component.Parameter()
    reserved_parameter_5 = Component.Parameter()
    reserved_parameter_6 = Component.Parameter()
    reserved_parameter_7 = Component.Parameter()
    m_pump_guess = Component.Parameter()
    ff1_guess = Component.Parameter()
    ff2_guess = Component.Parameter()
    ff3_guess = Component.Parameter()

    m_pump_input = Component.Input()
    reserved_input_2 = Component.Input()
    reserved_input_3 = Component.Input()
    reserved_input_4 = Component.Input()
    reserved_input_5 = Component.Input()
    reserved_input_6 = Component.Input()
    reserved_input_7 = Component.Input()
    reserved_input_8 = Component.Input()
    reserved_input_9 = Component.Input()
    reserved_input_10 = Component.Input()
    reserved_input_11 = Component.Input()
    reserved_input_12 = Component.Input()
    reserved_input_13 = Component.Input()
    p14 = Component.Input()
    p15 = Component.Input()
    p16 = Component.Input()
    p17 = Component.Input()
    p18 = Component.Input()
    p19 = Component.Input()
    p20 = Component.Input()
    p21 = Component.Input()
    reserved_input_22 = Component.Input()
    reserved_input_23 = Component.Input()
    reserved_input_24 = Component.Input()
    reserved_input_25 = Component.Input()
    recirculation_signal = Component.Input()
    reserved_input_27 = Component.Input()
    reserved_input_28 = Component.Input()

    m_pump_out = Component.Output()
    m0_out = Component.Output()
    m1_out = Component.Output()
    m2_out = Component.Output()
    m3_out = Component.Output()
    m4_out = Component.Output()
    m5_out = Component.Output()
    m6_out = Component.Output()
    m7_out = Component.Output()
    m8_out = Component.Output()
    m9_out = Component.Output()
    signal_out = Component.Output()

    @staticmethod
    def _frac(a, b, default):
        d = abs(a) + abs(b)
        if d <= 1.0e-12:
            return default
        return max(min(abs(a) / d, 1.0), 0.0)

    def calculate(self):
        m_pump = self.m_pump_input.v if self.m_pump_input.v == self.m_pump_input.v else self.m_pump_guess.v
        ff1 = self.ff1_guess.v
        ff2 = self.ff2_guess.v
        ff3 = self.ff3_guess.v
        ff4 = 0.5

        if not self.model.is_first_step:
            ff1 = self._frac(self.p14.v, self.p15.v, ff1)
            ff2 = self._frac(self.p16.v, self.p17.v, ff2)
            ff3 = self._frac(self.p18.v, self.p19.v, ff3)
            ff4 = self._frac(self.p20.v, self.p21.v, ff4)

        m1 = ff1 * m_pump
        m0 = (1.0 - ff1) * m_pump
        m2 = ff2 * m1
        m3 = (1.0 - ff2) * m1
        m4 = ff3 * m2
        m5 = (1.0 - ff3) * m2
        m6 = ff4 * m3
        m7 = (1.0 - ff4) * m3
        m8 = m4 + m6 + m0
        m9 = m5 + m7

        self.m_pump_out.v = m_pump
        self.m0_out.v = m0
        self.m1_out.v = m1
        self.m2_out.v = m2
        self.m3_out.v = m3
        self.m4_out.v = m4
        self.m5_out.v = m5
        self.m6_out.v = m6
        self.m7_out.v = m7
        self.m8_out.v = m8
        self.m9_out.v = m9
        self.signal_out.v = self.recirculation_signal.v
