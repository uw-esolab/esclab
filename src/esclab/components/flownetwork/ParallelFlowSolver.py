"""Type 4050 parallel-flow solver converted from Fortran."""

from esclab.simulate import Component


class ParallelFlowSolver(Component):
    """
    TRNSYS Type 4050: ESOL4050-Parallel-Flow-Solver.

    Parameters
    ----------
    solver, fluid_id, coef_a, coef_b, coef_c, n_pumps : float

    Inputs
    ------
    input_1..input_24 : float
        Branch flows, node pressures, pump state, and temperature.

    Outputs
    -------
    f_23, f_45, f_67 : float
        Flow fractions for tee branches.
    """

    solver = Component.Parameter()
    fluid_id = Component.Parameter()
    coef_a = Component.Parameter()
    coef_b = Component.Parameter()
    coef_c = Component.Parameter()
    n_pumps = Component.Parameter()

    for _idx in range(1, 25):
        locals()[f"input_{_idx}"] = Component.Input()

    f_23 = Component.Output()
    f_45 = Component.Output()
    f_67 = Component.Output()

    def _frac(self, a, b):
        denom = abs(a) + abs(b)
        if denom < 1.0e-12:
            return 0.5
        return max(min(abs(a) / denom, 1.0), 0.0)

    def calculate(self):
        m2, m3 = self.input_2.v, self.input_3.v
        m4, m5 = self.input_4.v, self.input_5.v
        m6, m7 = self.input_6.v, self.input_7.v

        if self.solver.v in (1.0, 2.0):
            p2, p3 = self.input_11.v, self.input_12.v
            p4, p5 = self.input_13.v, self.input_14.v
            p6, p7 = self.input_15.v, self.input_16.v
            self.f_23.v = self._frac(max(p2, 0.0), max(p3, 0.0))
            self.f_45.v = self._frac(max(p4, 0.0), max(p5, 0.0))
            self.f_67.v = self._frac(max(p6, 0.0), max(p7, 0.0))
        else:
            self.f_23.v = self._frac(m2, m3)
            self.f_45.v = self._frac(m4, m5)
            self.f_67.v = self._frac(m6, m7)

