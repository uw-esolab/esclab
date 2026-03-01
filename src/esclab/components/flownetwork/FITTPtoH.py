"""Type 6022 FIT-TPtoH converted from Fortran."""

from eeslib import fluid_properties as fp

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
        Specific enthalpy [J/kg].
    """

    temperature = Component.Input()
    pressure = Component.Input()
    if_sat = Component.Input()
    use_temp = Component.Input()
    x_value = Component.Input()

    enthalpy = Component.Output()

    @staticmethod
    def _safe(value, default):
        return value if value == value else default

    def calculate(self):
        # Read the Inputs
        temperature = max(self._safe(self.temperature.v, 273.15), 273.15)
        pressure = max(self._safe(self.pressure.v, 101325.0), 1.0)
        if_sat = self._safe(self.if_sat.v, 0.0)
        use_temp = self._safe(self.use_temp.v, 1.0)
        x_value = min(max(self._safe(self.x_value.v, 0.0), 0.0), 1.0)

        try:
            if if_sat == 1.0:
                # Wants an enthalpy value in the saturation region
                if use_temp == 1.0:
                    # Use Temp and Quality to find enthalpy (FIT_TQ equivalent)
                    h = float(fp.enthalpy("water", T=temperature, Q=x_value))
                else:
                    # Use pressure and quality to find enthalpy (FIT_PQ equivalent)
                    h = float(fp.enthalpy("water", P=pressure, Q=x_value))
            else:
                # Enthalpy value is not in the saturation region (FIT_TP equivalent)
                h = float(fp.enthalpy("water", T=temperature, P=pressure))
        except Exception:
            # Conservative fallback if property backend cannot resolve a state.
            cp_liq = 4200.0
            h_fg = 2.257e6
            h = cp_liq * max(temperature - 273.15, 0.0)
            if if_sat == 1.0:
                h += x_value * h_fg

        # Set the Outputs from this Model
        self.enthalpy.v = h
