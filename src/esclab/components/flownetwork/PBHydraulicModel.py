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

    for _idx in range(1, 18):
        locals()[f"input_{_idx}"] = Component.Input()
    for _idx in range(1, 13):
        locals()[f"output_{_idx}"] = Component.Output()

    def _clipf(self, x):
        return max(min(x, 1.0), 0.0)

    def calculate(self):
        m_in = self.input_1.v if self.input_1.v == self.input_1.v else max(self.m_guess.v, 0.0)
        p_in, p2, p3, p4, p5, p6, p7 = self.input_10.v, self.input_11.v, self.input_12.v, self.input_13.v, self.input_14.v, self.input_15.v, self.input_16.v

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

        vals = [m2, m3, m4, m5, m6, m7, m8, m9, ff1, ff2, ff3, ff4]
        for idx, val in enumerate(vals, start=1):
            getattr(self, f"output_{idx}").v = val

