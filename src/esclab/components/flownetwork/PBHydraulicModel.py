"""Type 6030 power-block hydraulic split solver converted from Fortran."""

from esclab.simulate import Component


class PBHydraulicModel(Component):
    """
    TRNSYS Type 6030: ESOL6030-PB-HydraulicModel.

    Parameters
    ----------
    fluid_id, m_guess, ff1_guess, ff2_guess, ff3_guess, ff4_guess, lr : float

    Inputs
    ------
    input_1..input_17 : float
        Inlet flow, branch flows, branch pressures, and inlet temperature.

    Outputs
    -------
    output_1..output_12 : float
        Branch flows and split fractions.
    """

    fluid_id = Component.Parameter()
    m_guess = Component.Parameter()
    ff1_guess = Component.Parameter()
    ff2_guess = Component.Parameter()
    ff3_guess = Component.Parameter()
    ff4_guess = Component.Parameter()
    lr = Component.Parameter(0.5)

    m_in_input = Component.Input()
    reserved_input_2 = Component.Input()
    reserved_input_3 = Component.Input()
    reserved_input_4 = Component.Input()
    reserved_input_5 = Component.Input()
    reserved_input_6 = Component.Input()
    reserved_input_7 = Component.Input()
    reserved_input_8 = Component.Input()
    reserved_input_9 = Component.Input()
    reserved_input_10 = Component.Input()
    p2_in = Component.Input()
    p3_in = Component.Input()
    p4_in = Component.Input()
    p5_in = Component.Input()
    p6_in = Component.Input()
    p7_in = Component.Input()
    reserved_input_17 = Component.Input()

    m2_out = Component.Output()
    m3_out = Component.Output()
    m4_out = Component.Output()
    m5_out = Component.Output()
    m6_out = Component.Output()
    m7_out = Component.Output()
    m8_out = Component.Output()
    m9_out = Component.Output()
    ff1_out = Component.Output()
    ff2_out = Component.Output()
    ff3_out = Component.Output()
    ff4_out = Component.Output()

    def _clipf(self, x):
        return max(min(x, 1.0), 0.0)

    def calculate(self):
        m_in = self.m_in_input.v if self.m_in_input.v == self.m_in_input.v else max(self.m_guess.v, 0.0)
        p2, p3, p4, p5, p6, p7 = self.p2_in.v, self.p3_in.v, self.p4_in.v, self.p5_in.v, self.p6_in.v, self.p7_in.v

        ff1 = self._clipf(self.ff1_guess.v if self.model.is_first_step else (p3 / max(p2 + p3, 1.0e-9) if p2 + p3 != 0 else self.ff1_guess.v))
        m2 = m_in * ff1
        m3 = m_in - m2
        ff2 = self._clipf(self.ff2_guess.v if self.model.is_first_step else (p5 / max(p4 + p5, 1.0e-9) if p4 + p5 != 0 else self.ff2_guess.v))
        m4 = m2 * ff2
        m5 = m2 - m4
        ff3 = self._clipf(self.ff3_guess.v if self.model.is_first_step else (p7 / max(p6 + p7, 1.0e-9) if p6 + p7 != 0 else self.ff3_guess.v))
        m6 = m4 * ff3
        m7 = m4 - m6
        ff4 = self._clipf(self.ff4_guess.v if self.model.is_first_step else 0.5)
        m8 = m4 * ff4
        m9 = m4 - m8

        self.m2_out.v = m2
        self.m3_out.v = m3
        self.m4_out.v = m4
        self.m5_out.v = m5
        self.m6_out.v = m6
        self.m7_out.v = m7
        self.m8_out.v = m8
        self.m9_out.v = m9
        self.ff1_out.v = ff1
        self.ff2_out.v = ff2
        self.ff3_out.v = ff3
        self.ff4_out.v = ff4

