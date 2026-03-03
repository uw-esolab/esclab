"""Power Block Valve component model (Type 6001)."""

from eeslib import fluid_properties as fp
from esclab.components.flownetwork.esol6015_helpers import PB_CV_data

from esclab.simulate import Component


class PowerBlockValve(Component):
    """
    # Object: Valve
    # Simulation Studio Model: ESOL6001
    """

    #    PARAMETERS
    Valve_diameter = Component.Parameter()  # diameter of valve [m]
    Valve_speed = Component.Parameter()     # valve speed [deg/s]
    Valve_type = Component.Parameter()

    #    INPUTS
    m_dot_in = Component.Input()   # Expected in kg/s
    P_in = Component.Input()       # Expected in Pa
    h_in = Component.Input()       # Expected in J/kg
    VP_input = Component.Input()   # Percent Open; 0 = Fully Closed, 1 = Fully Open

    #    OUTPUTS
    m_dot_out = Component.Output()   # Mass flow leaving valve [kg/s]
    Vol_dot_out = Component.Output() # Volumetric flow rate leaving the valve [m^3/s]
    P_out = Component.Output()       # Pressure leaving valve [Pa]
    h_out = Component.Output()       # Enthalpy Leaving valve [J/kg]
    T_out = Component.Output()       # Temperature leaving valve [K]
    DELTA_P = Component.Output()     # Pressure drop across valve [Pa]
    VP_output = Component.Output()   # output valve position, only changes between timestep transitions

    def calculate(self):

        # Do All of the "First Timestep Manipulations" Here - There Are No Iterations at the Initial Time
        if self.model.is_first_step:

            VP_input = self.VP_input.v

            if VP_input <= 0.00001:
                VP_input = 0.00001

            # Find pressure drop across valve
            # Flow Coefficient based on valve type, valve diameter and valve position
            CV = PB_CV_data(Valve_type=int(self.Valve_type.v), D_in=self.Valve_diameter.v, Valve_position=VP_input)

            # USED FOR VALIDATION OF FEEDWATER PUMPS
            # if (Valve_type /= 5) Then !Call Valve Function
            # else !6 inch Globe Valve
            #      if (VP_input >= .8) Then
            #       CV = 569.81*VP_input - 341.9
            #       !CV = 9293.4d0*VP_output**6.d0-23966.d0*VP_output**5.d0+24158.d0*VP_output**4.d0-11717.d0*VP_output**3.d0+2632.d0*VP_output**2.d0-70.393d0*VP_output+8.3633
            #
            #     else if (VP_input >= 0.025d0) Then
            #         CV = 9293.4d0*VP_input**6.d0-23966.d0*VP_input**5.d0+24158.d0*VP_input**4.d0-11717.d0*VP_input**3.d0+2632.d0*VP_input**2.d0-70.393d0*VP_input+8.3633
            #     else
            #         CV = VP_input * 7.8d0/0.025d0
            #     endif
            # endif

            # CONVERTED-NEEDS UNITS CHECK: P=101325.0 is already in Pa; eeslib uses Pa
            v_ref = fp.specific_volume("water", T=288.5, P=101325.0)  # density of water at reference state (atm pressure and 15.5 deg C)
            # CONVERTED-NEEDS UNITS CHECK: removed /1000 Pa->kPa and /1000 J/kg->kJ/kg; eeslib uses SI
            v_in = fp.specific_volume("water", P=self.P_in.v, h=self.h_in.v)  # density of water entering valve

            if v_in > 0.0:  # ERROR with FIT, set spec_grav to 1.d0 for this iteration
                spec_grav = v_in / v_ref  # specific gravity (actual density over reference density)
            else:  # specific gravity value not
                v_in = 1.0 / 1000.0
                spec_grav = 1.0
            Vol_in = self.m_dot_in.v * v_in                        # Volumetric flow entering the valve [m^3/s]
            Vol_in_gpm = Vol_in * 15850.3                           # Volumetric flow entering valve [GPM]
            delta_P_psi = spec_grav / (CV / Vol_in_gpm) ** 2.0     # pressure drop in psi across valve
            DELTA_P = delta_P_psi * 6894.76                         # pressure drop in Pa across valve
            P_out = self.P_in.v - DELTA_P

            if VP_input == 0.00001:
                VP_input = 0.0

            # TODO-NEEDS CONVERSION REVIEW: Original Fortran assigns h_in = h_out here, but h_out has not
            # yet been assigned at this point in the code. This appears to be a bug in the original Fortran
            # (h_out is uninitialized). Ported as-is.
            h_in_local = self.h_out.v
            # CONVERTED-NEEDS UNITS CHECK: removed /1000 Pa->kPa and /1000 J/kg->kJ/kg; eeslib uses SI
            rho = fp.density("water", P=P_out, h=self.h_out.v)
            if rho > 0.0:
                Vol_dot_out = self.m_dot_in.v / rho
            else:
                Vol_dot_out = self.m_dot_in.v / 1000.0

            # Set the Initial Values of the Outputs (#,Value)
            self.m_dot_out.v = self.m_dot_in.v  # 1  mass leaving the valve [kg/s]
            # TODO-NEEDS CONVERSION REVIEW: Output index 2 appears twice in the original Fortran source.
            # The first SetOutputValue(2, Vol_dot_out) is immediately overwritten by SetOutputValue(2, P_out).
            # This is a bug in the original Fortran (Type 6001). Both assignments are ported as-is below.
            self.Vol_dot_out.v = Vol_dot_out    # 2  Volumetric Flow rate leaving the valve [m^3/s]  (duplicate index 2 in Fortran)
            self.P_out.v = P_out                # 2  Pressure leaving the valve [Pa]  (duplicate index 2 in Fortran - overwrites Vol_dot_out above)
            self.h_out.v = 0.0                  # 3  Enthalpy leaving the valve [J/kg]  NOTE: Fortran sets output 3 = h_out (uninitialized)
            self.T_out.v = 0.0                  # 4  Temperature leaving the valve [K]
            self.DELTA_P.v = DELTA_P            # 5  Pressure Drop Across Valve [Pa]
            self.VP_output.v = VP_input         # 6  set initial outlet VP as desired input

            return

        # *** PERFORM ALL THE CALCULATION HERE FOR THIS MODEL. ***
        # Convert Valve Diameter from inches to m
        if self.model.timestep_iteration >= 1:  # Any other iteration besides the first in a timestep
            VP_output = self.VP_output.v
            if VP_output <= 0.00001:
                VP_output = 0.00001

            # Calculate pressure drop across valve
            # Flow Coefficient based on valve type, valve diameter and valve position
            CV = PB_CV_data(Valve_type=int(self.Valve_type.v), D_in=self.Valve_diameter.v, Valve_position=VP_output)

            # USED FOR VALIDATION OF FEEDWATER PUMPS
            # if (Valve_type /= 3) Then !Call Valve Function
            # else !6 inch Globe Valve
            #     if (VP_output >= .8) Then
            #         CV = 569.81*VP_output - 341.9
            #         !CV = 9293.4d0*VP_output**6.d0-23966.d0*VP_output**5.d0+24158.d0*VP_output**4.d0-11717.d0*VP_output**3.d0+2632.d0*VP_output**2.d0-70.393d0*VP_output+8.3633
            #
            #       else if (VP_output >= 0.025d0) Then
            #           CV = 9293.4d0*VP_output**6.d0-23966.d0*VP_output**5.d0+24158.d0*VP_output**4.d0-11717.d0*VP_output**3.d0+2632.d0*VP_output**2.d0-70.393d0*VP_output+8.3633
            #       else
            #           CV = VP_output * 7.8d0/0.025d0
            #       endif
            # endif

            if self.P_in.v > 1000.0:  # Pressure is possible, calculate specific gravity using FIT
                v_ref = 0.0009600113091621002  # density of water at reference state (atm pressure and 15.5 deg C)
                # CONVERTED-NEEDS UNITS CHECK: removed /1000 Pa->kPa and /1000 J/kg->kJ/kg; eeslib uses SI
                v_in = fp.specific_volume("water", P=self.P_in.v, h=self.h_in.v)  # density of water entering valve
                spec_grav = v_in / v_ref  # specific gravity (actual density over reference density)
            else:
                spec_grav = 1.0  # guess until there is a high enough pressure to compute with FIT Model
            Vol_in = self.m_dot_in.v * v_in      # Volumetric flow entering the valve [m^3/s]
            Vol_in_gpm = Vol_in * 15850.3         # Volumetric flow entering valve [GPM]
            delta_P_psi = spec_grav / (CV / Vol_in_gpm) ** 2.0  # pressure drop in psi across valve
            DELTA_P = delta_P_psi * 6894.76       # pressure drop in Pa across valve

            # Pressure, Enthalpy, & Temp Exiting Valve
            P_out = self.P_in.v - DELTA_P
            h_out = self.h_in.v  # No work or heat is traveling through system boundary, enthalpy does not change
            if P_out > 1000.0:  # will work with FIT
                # CONVERTED-NEEDS UNITS CHECK: removed /1000 Pa->kPa and /1000 J/kg->kJ/kg; eeslib uses SI
                T_out = fp.temperature("water", P=P_out, h=h_out)
                rho = fp.density("water", P=P_out, h=h_out)
            else:
                T_out = 0.0  # set temperature to zero, won't be entered as input to any type so okay to do

            if VP_output == 0.00001:
                VP_output = 0.0

        elif self.model.timestep_iteration == 0:  # New Timestep, move valve position if requested by operator
            if self.model.time != self.model.settings.timestep:
                self.Valve_speed.v = self.Valve_speed.v               # received in units of deg/s
                VP_output = self.VP_output.v                          # previous timestep's valve position [%]
                VP_input_d = self.VP_input.v * 90.0                   # Valve Position Requested, converting from percent open to degrees
                VP_output_d = VP_output * 90.0                        # Last Timestep's Valve Position, converting from percent open to degrees
                ts = self.model.timestep * 3600.0                     # Convert timestep from hr to s
                if VP_input_d == VP_output_d:                         # If input requested matched last timestep's input then do nothing
                    VP_output = self.VP_input.v
                elif VP_input_d > VP_output_d:                        # If input requested is greater than last timestep's value, open valve position more
                    VP_output = min(VP_output_d + (self.Valve_speed.v * ts), VP_input_d) / 90.0
                else:                                                  # If input requested is lower than last timestep's value, close valve more
                    VP_output = max(VP_output_d - (self.Valve_speed.v * ts), VP_input_d) / 90.0
            else:  # First iteration of each timestep
                VP_output = self.VP_input.v  # set it as the first input entered

            if VP_output <= 0.00001:  # Make sure valve input is not at zero because it will cause divide by zero errors
                VP_output = 0.00001

            # Calculate pressure drop across valve
            # Flow Coefficient based on valve type, valve diameter and valve position
            CV = PB_CV_data(Valve_type=int(self.Valve_type.v), D_in=self.Valve_diameter.v, Valve_position=VP_output)

            if self.P_in.v > 1000.0:  # Pressure is possible, compute specific gravity
                v_ref = 0.0009600113091621002  # density of water at reference state (atm pressure and 15.5 deg C)
                # CONVERTED-NEEDS UNITS CHECK: removed /1000 Pa->kPa and /1000 J/kg->kJ/kg; eeslib uses SI
                v_in = fp.specific_volume("water", P=self.P_in.v, h=self.h_in.v)  # density of water entering valve
                spec_grav = v_in / v_ref   # specific gravity (actual density over reference density)
            else:
                spec_grav = 1.0  # assume until pressure is valid to use fit
            Vol_in = self.m_dot_in.v * v_in      # Volumetric flow entering the valve [m^3/s]
            Vol_in_gpm = Vol_in * 15850.3         # Volumetric flow entering valve [GPM]
            delta_P_psi = spec_grav / (CV / Vol_in_gpm) ** 2.0  # pressure drop in psi across valve
            DELTA_P = delta_P_psi * 6894.76       # pressure drop in Pa across valve

            # Pressure, Enthalpy, & Temp Exiting Valve
            P_out = self.P_in.v - DELTA_P
            h_out = self.h_in.v  # No work or heat is traveling through system boundary, enthalpy does not change
            if P_out > 1000.0:  # will work with FIT
                # CONVERTED-NEEDS UNITS CHECK: removed /1000 Pa->kPa and /1000 J/kg->kJ/kg; eeslib uses SI
                T_out = fp.temperature("water", P=P_out, h=h_out)
                rho = fp.density("water", P=P_out, h=h_out)
            else:
                T_out = 0.0  # set temperature to zero, won't be entered as input to any type so okay to do

            if VP_output == 0.00001:  # Return valve output back to 0
                VP_output = 0.0

        # calculating volumetric flow rates
        if rho > 0.0:
            Vol_dot_out = self.m_dot_in.v / rho
        else:
            Vol_dot_out = self.m_dot_in.v / 1000.0

        # Set the Outputs from this Model (#,Value)
        self.m_dot_out.v = self.m_dot_in.v   # 1  Mass flow leaving valve [kg/s]
        self.Vol_dot_out.v = Vol_dot_out      # 2  Volumetric flow rate leaving the valve [m^3/s]
        self.P_out.v = P_out                  # 3  Pressure leaving valve [Pa]
        self.h_out.v = h_out                  # 4  Enthalpy Leaving valve [J/kg]
        self.T_out.v = T_out                  # 5  Temperature leaving valve [K]
        self.DELTA_P.v = DELTA_P              # 6  Pressure drop across valve [Pa]
        self.VP_output.v = VP_output          # 7  output valve position, only changes between timestep transitions
