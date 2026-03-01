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

    solver_mode_signal = Component.Input()
    m2_in = Component.Input()
    m3_in = Component.Input()
    m4_in = Component.Input()
    m5_in = Component.Input()
    m6_in = Component.Input()
    m7_in = Component.Input()
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
    reserved_input_18 = Component.Input()
    reserved_input_19 = Component.Input()
    reserved_input_20 = Component.Input()
    reserved_input_21 = Component.Input()
    reserved_input_22 = Component.Input()
    reserved_input_23 = Component.Input()
    reserved_input_24 = Component.Input()

    f_23 = Component.Output()
    f_45 = Component.Output()
    f_67 = Component.Output()

    def _frac(self, a, b):
        denom = abs(a) + abs(b)
        if denom < 1.0e-12:
            return 0.5
        return max(min(abs(a) / denom, 1.0), 0.0)

    def calculate(self):
        m2, m3 = self.m2_in.v, self.m3_in.v
        m4, m5 = self.m4_in.v, self.m5_in.v
        m6, m7 = self.m6_in.v, self.m7_in.v

        if self.solver.v in (1.0, 2.0):
            p2, p3 = self.p2_in.v, self.p3_in.v
            p4, p5 = self.p4_in.v, self.p5_in.v
            p6, p7 = self.p6_in.v, self.p7_in.v
            self.f_23.v = self._frac(max(p2, 0.0), max(p3, 0.0))
            self.f_45.v = self._frac(max(p4, 0.0), max(p5, 0.0))
            self.f_67.v = self._frac(max(p6, 0.0), max(p7, 0.0))
        else:
            self.f_23.v = self._frac(m2, m3)
            self.f_45.v = self._frac(m4, m5)
            self.f_67.v = self._frac(m6, m7)

