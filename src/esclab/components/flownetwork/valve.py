"""Valve component model (Type 4007) and CV interpolation helper."""

import numpy as np
from scipy.interpolate import RectBivariateSpline

from esclab.simulate import Component
from esclab.components.esol_properties import Incompressible as Inc


class Valve(Component):
    """
    # Object: ESOL4007-Valve
    # Simulation Studio Model: ESOL4007-Valve
    #

    # Author: Matt Tuman
    # Editor:
    # Date:     October 27, 2022
    # last modified: October 27, 2022
    """
    trnsys_type = "4007"

    #    PARAMETERS
    Diameter = Component.Parameter()
    Fluid_ID = Component.Parameter()
    Valve_Type = Component.Parameter()
    Valve_speed = Component.Parameter()

    #    INPUTS
    m_dot = Component.Input()
    Pressure = Component.Input()
    Temperature = Component.Input()
    fraction_open = Component.Input()
    mass_counter = Component.Input()
    Cv = Component.Input()

    #    OUTPUTS
    m_dot_out = Component.Output()
    Pressure_out = Component.Output()
    Temperature_out = Component.Output()
    mass_counter_out = Component.Output()
    Cv_out = Component.Output()
    VP_output = Component.Output()

    # stored data
    fraction_open = float('nan')

    def calculate(self):

        if self.model.is_first_step:
            ## COMPUTE PRESSURE DROP
            ##########################
            # Compute Volumetric Flow Rate
            rho_fluid = Inc.density(self.Fluid_ID.v, self.Temperature.v, 0.0)
            Q = self.m_dot.v/rho_fluid
            # Convert flowrate to gpm
            Q = Q*15850.323140625002
            # Compute specific gravity of fluid
            SG = rho_fluid/1000.0
            # Compute Cv of valve
            Cv = CV_data(self.Valve_Type.v, self.Diameter.v, self.fraction_open)
            # Compute pressure drop [Pa]
            dP = SG*Q**2/(Cv**2) * 6894.76

            #Set the Initial Values of the Outputs (#,Value)
            self.m_dot_out.v = self.m_dot.v     #1
            self.Pressure_out.v = self.Pressure.v - dP      #2
            self.Temperature_out.v = self.Temperature.v     #3
            self.mass_counter_out.v = self.mass_counter.v       #4
            self.Cv_out.v = Cv      #5
            self.VP_output.v = self.fraction_open        #6

            # Set Initial Values of Dynamic Storage
            # Call SetDynamicArrayValueThisIteration(1, fraction_open) # Store valve position

            return

        # Iteration calculations
        # -------------------------------------------------------------------------------------------------------
        # VP_output = getOutputValue(6)

        # -----------------------------------------------------------------------------------------------------------------------
        if self.model.iteration == 0: # do not update flow rate
            if (not self.model.is_first_step): # not the first timestep
                fraction_open_d = self.fraction_open * 90.0         # Valve Position Requested, converting from percent open to degrees
                VP_output_d = self.VP_output.v * 90.0                 # Last Timesteps Valve Position, converting from percent open to degrees
                Timestep_s = self.model.timestep * 3600.0                 # convert timestep to seconds instead of hours
                if (self.VP_output.v == self.fraction_open):            # valve position is at the requested input
                    self.VP_output.v = self.fraction_open
                elif (VP_output_d > fraction_open_d): # current valve position is greater than requested input, close valve based on valve speed
                    VP_output_d = max(VP_output_d - self.Valve_speed.v * Timestep_s, fraction_open_d)
                    VP_output = VP_output_d/90.0
                else: # current valve position is less than requested input, open valve based on valve speed
                    VP_output_d = min(VP_output_d + self.Valve_speed.v * Timestep_s, fraction_open_d)
                    VP_output = VP_output_d/90.0
            else: # First Timestep of the simulation, set equal to the input value rather than intial value (needed if using forcing function types)
                VP_output = self.fraction_open

        #  Compute Volumetric Flow Rate
        rho_fluid = Inc.density(self.Fluid_ID.v, self.Temperature.v, 0.0)
        Q = self.m_dot.v/rho_fluid
        #  Convert flowrate to gpm
        Q = Q*15850.323140625002
        #  Compute specific gravity of fluid
        SG = rho_fluid/1000.0
        #  Compute Cv of valve (Check if valve position has changed from last timestep)
        if self.model.timestep_iteration == 0:
            if self.Pos_last.v != VP_output:
                Cv = CV_data(self.Valve_Type.v, self.Diameter.v, VP_output)
        else:
            Cv = self.Cv_out.v
        #  Compute pressure drop [Pa]
        dP = SG*Q**2/(Cv**2) * 6894.76

        # -----------------------------------------------------------------------------------------------------------------------
        # Set the Outputs from this Model (#,Value)
        self.m_dot_out.v = self.m_dot.v     #1
        self.Pressure_out.v = self.Pressure.v - dP      #2
        self.Temperature_out.v = self.Temperature.v     #3
        self.mass_counter_out.v = self.mass_counter.v       #4
        self.Cv_out.v = Cv      #5
        self.VP_output.v = VP_output        #6
        # --------------------------------------------------------------------------------------------------------

    def converged(self):
        # Store valve position for next iteration
        self.fraction_open = self.VP_output.v
        return


def CV_data(Valve_Type, D_in, Valve_position):
    """
    Compute C_v for various valve types using 2D interpolation.

    Parameters
    ----------
    Valve_Type : int
        1 = Concentric Butterfly Valve
        2 = Triple Offset Butterfly Valve
    D_in : float
        Valve diameter [m]
    Valve_position : float
        Valve opening position [0-1], where 1 is fully open

    Returns
    -------
    float
        Flow coefficient C_v

    References
    ----------
    https://www.valvesonline.com.au/references/flow-rates/butterfly-valves/
    """
    C_v_min = 0.0001  # minimum CV allowed

    # Return minimum flow to prevent hydraulic model crashes
    if Valve_position == 0:
        return C_v_min

    D = D_in * 39.3701  # Convert m -> inches

    # Define valve data: diameters, positions, and CV values
    valve_data = {
        1: {  # Concentric Butterfly
            'D': np.array([4.0, 8.0, 12.0, 16.0, 20.0, 24.0, 32.0]),
            'Pos': np.array([0.0, 1/9, 2/9, 3/9, 4/9, 5/9, 6/9, 7/9, 8/9, 1.0]),
            'CV': np.array([
                [0.0, 0.5, 17.0, 36.0, 78.0, 139.0, 230.0, 364.0, 546.0, 600.0],
                [0.0, 3.0, 89.0, 188.0, 408.0, 727.0, 1202.0, 1903.0, 2854.0, 3136.0],
                [0.0, 5.0, 234.0, 495.0, 1072.0, 1911.0, 3162.0, 5005.0, 7507.0, 8250.0],
                [0.0, 8.0, 464.0, 983.0, 2130.0, 3797.0, 6282.0, 9942.0, 14913.0, 16388.0],
                [0.0, 14.0, 791.0, 1674.0, 3628.0, 6465.0, 10698.0, 16931.0, 25396.0, 27908.0],
                [0.0, 22.0, 1222.0, 2587.0, 5605.0, 9989.0, 16528.0, 26157.0, 39236.0, 43116.0],
                [0.0, 45.0, 2387.0, 4791.0, 8736.0, 13788.0, 20613.0, 31395.0, 48117.0, 68250.0],
            ])
        },
        2: {  # Triple Offset Butterfly
            'D': np.array([4.0, 8.0, 12.0, 16.0, 20.0, 24.0, 32.0]),
            'Pos': np.array([0.0, 1/9, 2/9, 3/9, 4/9, 5/9, 6/9, 7/9, 8/9, 1.0]),
            'CV': np.array([
                [0.0, 8.4, 29.3, 58.6, 92.0, 140.0, 200.0, 330.0, 370.0, 420.0],
                [0.0, 38.2, 140.0, 270.0, 420.0, 640.0, 900.0, 1500.0, 1690.0, 1920.0],
                [0.0, 88.4, 310.0, 620.0, 980.0, 1460.0, 2080.0, 3450.0, 3890.0, 4420.0],
                [0.0, 150.0, 530.0, 1060.0, 1660.0, 2490.0, 3540.0, 5870.0, 6620.0, 7520.0],
                [0.0, 270.0, 930.0, 1850.0, 2900.0, 4350.0, 6190.0, 10300.0, 11600.0, 13200.0],
                [0.0, 420.0, 1450.0, 2890.0, 4530.0, 6800.0, 9680.0, 16100.0, 18100.0, 20600.0],
                [0.0, 810.0, 2820.0, 5630.0, 8840.0, 13300.0, 18900.0, 31400.0, 35400.0, 40200.0],
            ])
        }
    }

    if Valve_Type not in valve_data:
        raise ValueError(f"Valve_Type {Valve_Type} not supported (1=Concentric, 2=Triple Offset)")

    data = valve_data[Valve_Type]
    D_vals = data['D']
    Pos_vals = data['Pos']
    CV_vals = data['CV']

    # Create 2D bivariate spline interpolator (cubic by default)
    # RectBivariateSpline expects data with x as columns, y as rows
    spl = RectBivariateSpline(D_vals, Pos_vals, CV_vals, kx=1, ky=1)

    # Evaluate at the requested diameter and position
    C_v = float(spl(D, Valve_position)[0, 0])

    # Ensure minimum CV
    return max(C_v, C_v_min)
