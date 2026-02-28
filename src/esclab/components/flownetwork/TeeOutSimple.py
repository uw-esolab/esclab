"""Type 4015 tee passthrough converted from Fortran."""

from esclab.simulate import Component


class TeeOutSimple(Component):
    """
    TRNSYS Type 4015: ESOL4015-TeeOut-Simple.

    Parameters
    ----------
    None.

    Inputs
    ------
    temperature, pressure, mass_counter : float
        Through variables.

    Outputs
    -------
    temperature_out, pressure_out, mass_counter_out : float
        Direct passthrough outputs.
    """

    temperature = Component.Input()
    pressure = Component.Input()
    mass_counter = Component.Input()

    temperature_out = Component.Output()
    pressure_out = Component.Output()
    mass_counter_out = Component.Output()

    def calculate(self):
        self.temperature_out.v = self.temperature.v
        self.pressure_out.v = self.pressure.v
        self.mass_counter_out.v = self.mass_counter.v
