"""Water enthalpy lookup component model (Type 6022)."""

from eeslib import fluid_properties as fp

from esclab.simulate import Component


class WaterEnthalpyLookup(Component):
    """
    Object: FIT-TPtoH
    Simulation Studio Model: ESOL6022-FIT-TPtoH

    Computes water/steam enthalpy given temperature, pressure, and an optional
    saturation quality. Supports superheated, compressed-liquid, and saturated
    (two-phase) states.

    Inputs
    ------
    temperature : float
        Temperature used to find enthalpy [K].
    pressure : float
        Pressure used to find enthalpy [Pa].
    if_sat : float
        1 = use saturation correlations; 0 = use T-P (non-saturated) state.
    use_temp : float
        If 1 and if_sat==1, use temperature + quality to find enthalpy;
        otherwise use pressure + quality.
    x_value : float
        Steam quality desired when in saturation region [-].

    Outputs
    -------
    enthalpy : float
        Fluid enthalpy [J/kg].
    """

    #    INPUTS
    temperature = Component.Input()  # Temperature used to find enthalpy
    pressure = Component.Input()     # Pressure used to find enthalpy
    if_sat = Component.Input()       # if you want saturated condition
    use_temp = Component.Input()     # If value inputted is 1, temperature will be used to find enthalpy for saturated water, else pressure will be used
    x_value = Component.Input()      # Quality value desired if saturated

    #    OUTPUTS
    enthalpy = Component.Output()  # enthalpy [J/kg]

    def _compute_enthalpy(self):
        if self.if_sat.v == 1.0:  # Wants a enthalpy value in the saturation region
            if self.use_temp.v == 1.0:  # use Temp and Quality to find enthalpy
                # TODO-NEEDS LIBRARY: eeslib equivalent of FIT_TQ("water", T, Q, enth) for saturation-region enthalpy
                enthalpy = fp.enthalpy("water", T=self.temperature.v, Q=self.x_value.v)
                # TODO-NEEDS UNITS CHECK: FIT_TQ returns kJ/kg; multiplied by 1000 to convert to J/kg
                enthalpy = enthalpy * 1000.0
            else:  # use pressure and quality to find enthalpy
                # TODO-NEEDS LIBRARY: eeslib equivalent of FIT_PQ("water", P, Q, enth) for saturation-region enthalpy
                # TODO-NEEDS UNITS CHECK: FIT returns kJ/kg; incoming pressure was divided by 1000 in Fortran (Pa->kPa)
                enthalpy = fp.enthalpy("water", P=self.pressure.v / 1000.0, Q=self.x_value.v)
                # TODO-NEEDS UNITS CHECK: FIT_PQ returns kJ/kg; multiplied by 1000 to convert to J/kg
                enthalpy = enthalpy * 1000.0
        else:  # enthalpy value is not in the saturation region
            # TODO-NEEDS UNITS CHECK: FIT_TP("water", T, P/1000, enth); incoming pressure converted Pa->kPa
            enthalpy = fp.enthalpy("water", T=self.temperature.v, P=self.pressure.v / 1000.0)
            # TODO-NEEDS UNITS CHECK: FIT_TP returns kJ/kg; multiplied by 1000 to convert to J/kg
            enthalpy = enthalpy * 1000.0
        return enthalpy

    def presim_setup(self, **kwargs):
        # Set the Initial Values of the Outputs (#,Value)
        self.enthalpy.v = self._compute_enthalpy()  # enthalpy

    def calculate(self):
        if self.model.is_first_step:
            # Set the Initial Values of the Outputs (#,Value)
            self.enthalpy.v = self._compute_enthalpy()  # enthalpy
            return

        # Set the Outputs from this Model (#,Value)
        self.enthalpy.v = self._compute_enthalpy()  # enthalpy
