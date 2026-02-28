"""Type 6022 FIT-TPtoH surrogate converted from Fortran."""

from esclab.simulate import Component


class FITTPtoH(Component):
    """
    TRNSYS Type 6022: ESOL6022-FIT-TPtoH.

    Parameters
    ----------
    None.

    Inputs
    ------
    temperature, pressure, if_sat, use_temp, x_value : float

    Outputs
    -------
    enthalpy : float
        Estimated specific enthalpy [J/kg].
    """

    temperature = Component.Input()
    pressure = Component.Input()
    if_sat = Component.Input()
    use_temp = Component.Input()
    x_value = Component.Input()

    enthalpy = Component.Output()

    def calculate(self):
        t = self.temperature.v
        p = self.pressure.v
        x = max(min(self.x_value.v, 1.0), 0.0)
        cp_liq = 4200.0
        h_fg = 2.257e6
        if self.if_sat.v == 1.0:
            # TODO: Replace with FIT_TQ/FIT_PQ property calls when available in Python stack.
            h = cp_liq * max(t - 273.15, 0.0) + x * h_fg
        else:
            h = cp_liq * max(t - 273.15, 0.0)
        if p < 1000.0:
            h = self.enthalpy.v if self.enthalpy.v == self.enthalpy.v else h
        self.enthalpy.v = h
