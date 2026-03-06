"""Power block piping component model (Type 6016)."""

import numpy as np
from eeslib import fluid_properties as fp

from esclab.simulate import Component
from esclab.components.flownetwork.SimplePipe import FricFactor_IC


class PowerBlockPiping(Component):
    """
    Object: Power Block Piping
    Simulation Studio Model: ESOL6016-PB_Piping
    """
    trnsys_type = "6016"

    # PARAMETERS
    Pipe_ID = Component.Parameter()           # pipe inner diameter
    Pipe_length = Component.Parameter()       # pipe length
    Roughness = Component.Parameter()         # pipe roughness
    elevation_change = Component.Parameter()  # elevation change

    # INPUTS
    m_dot_in = Component.Input()   # mass flow rate in
    P_in = Component.Input()       # inlet pressure
    h_in = Component.Input()       # inlet enthalpy

    # OUTPUTS & VARIABLES
    m_dot_out = Component.Output()    # 1: Mass Flow Rate leaving the pipe
    vol_dot_out = Component.Output()  # 2: Volumetric Flow Rate leaving the pipe
    P_out = Component.Output()        # 3: P_out
    h_out = Component.Output()        # 4: h_out
    T_out = Component.Output()        # 5: T_out
    DELTA_P = Component.Output()      # 6: Pressure_loss
    ff_guess = Component.Output()     # 7: Friction Factor guess for next iteration

    def presim_setup(self, **kwargs):
        super().presim_setup(**kwargs)

        # Do All of the First Timestep Manipulations Here - There Are No Iterations at the Initial Time
        T_in = fp.temperature("water", P=self.P_in.v, h=self.h_in.v)  # CONVERTED-NEEDS UNITS CHECK: removed /1000 Pa->kPa and /1000 J/kg->kJ/kg; eeslib uses Pa and J/kg
        rho = fp.density("water", P=self.P_in.v, h=self.h_in.v)       # CONVERTED-NEEDS UNITS CHECK: removed /1000 Pa->kPa and /1000 J/kg->kJ/kg; eeslib uses Pa and J/kg
        guess = 10.  # setting friction factor guess value for next iteration

        if rho > 0.:
            self.vol_dot_out.v = self.m_dot_in.v / rho
        else:
            self.vol_dot_out.v = self.m_dot_in.v / 1000.

        # Set the Initial Values of the Outputs (#,Value)
        self.m_dot_out.v = self.m_dot_in.v    # Mass Flow Rate leaving the pipe
        # vol_dot_out already set above        # Volumetric Flow Rate leaving the pipe
        self.P_out.v = self.P_in.v            # P_out
        self.h_out.v = self.h_in.v            # h_out
        self.T_out.v = T_in                   # T_out
        self.DELTA_P.v = 0.                   # Pressure_loss
        self.ff_guess.v = guess               # Friction Factor guess for previous iterations

    def calculate(self):
        super().calculate()

        # Perform Any "After Convergence" Manipulations That May Be Required at the End of Each Timestep
        if self.model.is_converged:
            ff = 28.1
            self.ff_guess.v = ff  # Friction Factor guess for next iteration
            return

        # Read the Inputs
        guess = max(self.ff_guess.v, 0.1)

        # mass balance
        self.m_dot_out.v = self.m_dot_in.v

        # Calculating Friction Pressure Drop
        if self.P_in.v > 10000.:
            rho = fp.density("water", P=self.P_in.v, h=self.h_in.v)  # CONVERTED-NEEDS UNITS CHECK: removed /1000 Pa->kPa and /1000 J/kg->kJ/kg; eeslib uses Pa and J/kg
            if rho == 0.:
                rho = 1000.  # Pressure entering was negative, assume density is 1000 kg/s for this iteration
        else:
            rho = 1000.

        vel = self.m_dot_in.v / rho / (3.14 / 4. * self.Pipe_ID.v ** 2.)
        Re = rho * vel * self.Pipe_ID.v / 0.001  # Assuming viscosity is 0.001 [Pa-s] WILL REPLACE LATER

        ff = FricFactor_IC(self.Roughness.v / self.Pipe_ID.v, Re, guess)
        K_T = (8. * ff * self.Pipe_length.v) / ((3.14 ** 2.) * (self.Pipe_ID.v ** 5.) * rho)
        DELTA_P_fric = K_T * self.m_dot_in.v ** 2.

        # Calculating pressure drop based on elevation change
        if self.m_dot_in.v > 0.01:  # account for a pressure drop in the system
            DELTA_P_elevation = self.elevation_change.v * 9.81 * rho
        else:
            DELTA_P_elevation = 0.

        # Overall Pressure Drop
        self.DELTA_P.v = DELTA_P_fric + DELTA_P_elevation
        self.P_out.v = self.P_in.v - self.DELTA_P.v

        # energy balance - assuming no heat transfer with environment
        self.h_out.v = self.h_in.v
        if self.P_in.v > 10000.:  # make sure it can solve with FIT
            self.T_out.v = fp.temperature("water", P=self.P_out.v, h=self.h_out.v)  # CONVERTED-NEEDS UNITS CHECK: removed /1000 Pa->kPa and /1000 J/kg->kJ/kg; eeslib uses Pa and J/kg
        else:
            self.T_out.v = self.T_out.v  # use previous output value

        # volumetric flow rate
        self.vol_dot_out.v = self.m_dot_out.v / rho

        # Set the Outputs from this Model (#,Value)
        # self.m_dot_out.v  # 1: Mass Flow Rate Leaving the pipe
        # self.vol_dot_out.v  # 2: Volumetric Flow Rate leaving the pipe
        # self.P_out.v        # 3: P_out
        # self.h_out.v        # 4: h_out
        # self.T_out.v        # 5: T_out
        # self.DELTA_P.v      # 6: Pressure_loss
        self.ff_guess.v = ff  # 7: Friction Factor guess for next iteration
