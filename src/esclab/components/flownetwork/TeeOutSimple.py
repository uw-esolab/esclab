"""Simple tee-out junction pass-through component model (Type 4015)."""

from esclab.simulate import Component


class TeeOutSimple(Component):
    """
    Object: ESOL4015-TeeOut-Simple
    Simulation Studio Model: ESOL4015-TeeOut-Simple

    Author: Matt Tuman
    Date: January 29, 2024
    last modified: January 29, 2024

    A minimal pass-through junction that forwards Temperature, Pressure, and
    Mass Counter unchanged. No parameters; no hydraulic or thermal calculation.
    """

    #    INPUTS
    Temperature = Component.Input()   # Temperature [-]
    Pressure = Component.Input()      # Pressure [-]
    Mass_Counter = Component.Input()  # Mass Counter [-]

    #    OUTPUTS
    Temperature_out = Component.Output()   # Temperature [-]
    Pressure_out = Component.Output()      # Pressure [-]
    Mass_Counter_out = Component.Output()  # Mass Counter [-]

    def presim_setup(self, **kwargs):
        # Set the Initial Values of the Outputs (#,Value)
        self.Temperature_out.v = self.Temperature.v   # Temperature
        self.Pressure_out.v = self.Pressure.v         # Pressure
        self.Mass_Counter_out.v = self.Mass_Counter.v # Mass Counter

    def calculate(self):
        if self.model.is_first_step:
            # Set the Initial Values of the Outputs (#,Value)
            self.Temperature_out.v = self.Temperature.v   # Temperature
            self.Pressure_out.v = self.Pressure.v         # Pressure
            self.Mass_Counter_out.v = self.Mass_Counter.v # Mass Counter
            return

        # Set the Outputs from this Model (#,Value)
        self.Temperature_out.v = self.Temperature.v   # Temperature
        self.Pressure_out.v = self.Pressure.v         # Pressure
        self.Mass_Counter_out.v = self.Mass_Counter.v # Mass Counter
