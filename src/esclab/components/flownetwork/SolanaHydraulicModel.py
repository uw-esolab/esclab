"""Type 6031 Solana hydraulic model converted from Fortran."""

from esclab.simulate import Component


class SolanaHydraulicModel(Component):
    """
    TRNSYS Type 6031: SolanaHydraulicModel.

    Uses the Fortran first-step flow-fraction structure and a stable pressure
    ratio fallback for iterative updates.
    """

    for _idx in range(1, 12):
        locals()[f"parameter_{_idx}"] = Component.Parameter()
    for _idx in range(1, 29):
        locals()[f"input_{_idx}"] = Component.Input()
    for _idx in range(1, 13):
        locals()[f"output_{_idx}"] = Component.Output()

    @staticmethod
    def _frac(a, b, default):
        d = abs(a) + abs(b)
        if d <= 1.0e-12:
            return default
        return max(min(abs(a) / d, 1.0), 0.0)

    def calculate(self):
        m_pump = self.input_1.v if self.input_1.v == self.input_1.v else self.parameter_8.v
        ff1 = self.parameter_9.v
        ff2 = self.parameter_10.v
        ff3 = self.parameter_11.v
        ff4 = 0.5

        if not self.model.is_first_step:
            ff1 = self._frac(self.input_14.v, self.input_15.v, ff1)
            ff2 = self._frac(self.input_16.v, self.input_17.v, ff2)
            ff3 = self._frac(self.input_18.v, self.input_19.v, ff3)
            ff4 = self._frac(self.input_20.v, self.input_21.v, ff4)

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

        self.output_1.v = m_pump
        self.output_2.v = m0
        self.output_3.v = m1
        self.output_4.v = m2
        self.output_5.v = m3
        self.output_6.v = m4
        self.output_7.v = m5
        self.output_8.v = m6
        self.output_9.v = m7
        self.output_10.v = m8
        self.output_11.v = m9
        self.output_12.v = self.input_26.v
