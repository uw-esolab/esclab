"""Expansion system component model (Type 4012)."""

import math

import numpy as np

from esclab.simulate import Component
from esclab.components.esol_properties import Incompressible as Inc

from esclab.components.flownetwork.sf_piping_helpers import PressureDrop
from esclab.components.flownetwork.valve import CV_data
from esclab.components.flownetwork.sergio_scripts import matrixinv  # SergioScripts module

from eeslib.functions import convert

class ExpansionSystem(Component):
    """
    Object: ESOL4012-ExpansionSystem
    Simulation Studio Model: ESOL4012-ExpansionSystem

    Author: Matt Tuman
    Date:    January 25, 2023
    last modified: January 25, 2023

    Parameters
    ----------
    Fluid_ID : str
        Fluid identifier (e.g., 21 for Thermanol).
    P_ev : float
        Pressure in vapor space of expansion tank [Pa].
    P_of : float
        Pressure in vapor space of overflow tank [Pa].
    PC_a : float
        Pump curve coefficient a  (H_l = a + b*q + c*q^2).
    PC_b : float
        Pump curve coefficient b.
    PC_c : float
        Pump curve coefficient c.
    roughness : float
        Roughness of all pipes [-].
    N_of : float
        Number of overflow tanks [-].
    N_ev : float
        Number of expansion tanks [-].
    H_of : float
        Height of overflow tank [m].
    H_ev : float
        Height of expansion tank [m].
    H_of_init : float
        Initial fill fraction of HTF in overflow tank [-].
    D_of : float
        Diameter of overflow tank [m].
    D_ev : float
        Diameter of expansion tank [m].
    total_mass : float
        Total mass of HTF in the entire plant [kg].
    only_Exp : float
        1 = expansion tanks only; 0 = expansion and overflow tanks.
    valve_speed : float
        Speed at which valves open/close [deg/s].
    T_of : float
        Initial temperature of HTF in overflow tank [K or °C – see Fortran convention].
    T_ev : float
        Initial temperature of HTF in expansion tanks [K or °C – see Fortran convention].
    CV_orifice : float
        CV of orifice in bypass line.
    T_chill_out : float
        Temperature of HTF leaving the chiller.

    Inputs
    ------
    HV_204_09A : float
        Valve position [0-1] where 1 is fully open.
    LV_205_14 : float
        Valve position [0-1] where 1 is fully open.
    LV_204_05 : float
        Valve position [0-1] where 1 is fully open.
    LV_204_09B : float
        Valve position [0-1] where 1 is fully open.
    m_dot_in : float
        Mass flow into expansion vessels [kg/s].
    T_in : float
        Temperature of mass flow into expansion vessel [°C or K – Fortran convention].
    P_in : float
        Pressure of the HTF entering the expansion vessel [Pa].
    m_counter : float
        Mass counter of system [kg].

    Outputs
    -------
    level_ev : float
        Level in expansion vessel [fraction 0-1].
    level_of : float
        Level in overflow tank [fraction 0-1].
    m_dot_in_out : float
        Mass flow in (passed through) [kg/s].
    T_in_out : float
        Temperature of inflow (passed through).
    P_bot_ev : float
        Pressure at bottom of expansion tank [Pa].
    net_of_flow : float
        Net mass flow to overflow tank (m_in_of - m_out_of) [kg/s].
    m_dot_exp : float
        Mass flow to/from expansion tank [kg/s].
    m_dot_recirc : float
        Recirculation mass flow rate [kg/s].
    m_dot_chiller : float
        Chiller mass flow rate [kg/s].
    T_of_out : float
        Temperature of HTF in overflow tank.
    T_ev_out : float
        Temperature of HTF in expansion tanks.
    alarm_high_high_pressure : float
        High-high pressure alarm flag (0 or 1).
    """
    trnsys_type = "4012"

    pi = 3.1415927

    # *** Model Parameters ***
    Fluid_ID = Component.Parameter()      # Fluid identifier (21 for Thermanol)
    P_ev = Component.Parameter()          # Pressure in vapor space of expansion tank [Pa]
    P_of = Component.Parameter()          # Pressure in vapor space of overflow tank [Pa]
    PC_a = Component.Parameter()          # Pump curve coefficient a
    PC_b = Component.Parameter()          # Pump curve coefficient b
    PC_c = Component.Parameter()          # Pump curve coefficient c
    roughness = Component.Parameter()     # Roughness of all pipes [-]
    N_of = Component.Parameter()          # Number of overflow tanks [-]
    N_ev = Component.Parameter()          # Number of expansion tanks [-]
    H_of = Component.Parameter()          # Height of overflow tank [m]
    H_ev = Component.Parameter()          # Height of expansion tank [m]
    H_of_init = Component.Parameter()     # Initial fill fraction of HTF in overflow
    D_of = Component.Parameter()          # Diameter of overflow tank [m]
    D_ev = Component.Parameter()          # Diameter of expansion tank [m]
    total_mass = Component.Parameter()    # Total mass of HTF in entire plant [kg]
    only_Exp = Component.Parameter()      # 1 = expansion only; 0 = expansion and overflow
    valve_speed = Component.Parameter()   # Speed at which valves open/close [deg/s]
    T_of = Component.Parameter()          # Initial temperature of HTF in overflow
    T_ev = Component.Parameter()          # Initial temperature of HTF in expansion tanks
    CV_orifice = Component.Parameter()    # CV of orifice in bypass line
    T_chill_out = Component.Parameter()   # Temperature of HTF leaving the chiller

    # *** Model Inputs ***
    HV_204_09A = Component.Input()        # Valve position [0-1] fully open
    LV_205_14 = Component.Input()         # Valve position [0-1] fully open
    LV_204_05 = Component.Input()         # Valve position [0-1] fully open
    LV_204_09B = Component.Input()        # Valve position [0-1] fully open
    m_dot_in = Component.Input()          # Mass flow into expansion vessels [kg/s]
    T_in = Component.Input()              # Temperature of mass flow into expansion vessel
    P_in = Component.Input()              # Pressure of HTF entering expansion vessel [Pa]
    m_counter = Component.Input()         # Mass counter of system [kg]

    # *** Model Outputs ***
    level_ev = Component.Output()              # Level in expansion vessel [0-1]
    level_of = Component.Output()              # Level in overflow tank [0-1]
    m_dot_in_out = Component.Output()          # Mass flow (passed through) [kg/s]
    T_in_out = Component.Output()              # Temperature of inflow (passed through)
    P_bot_ev = Component.Output()              # Pressure at bottom of expansion tank [Pa]
    net_of_flow = Component.Output()           # Net mass flow to overflow (m_in_of - m_out_of) [kg/s]
    m_dot_exp = Component.Output()             # Mass flow to expansion tank [kg/s]
    m_dot_recirc = Component.Output()          # Recirculation mass flow rate [kg/s]
    m_dot_chiller = Component.Output()         # Chiller mass flow rate [kg/s]
    T_of_out = Component.Output()              # Temperature of HTF in overflow tank
    T_ev_out = Component.Output()              # Temperature of HTF in expansion tanks
    alarm_high_high_pressure = Component.Output()  # High-high pressure alarm (0 or 1)

    def presim_setup(self, **kwargs):
        super().presim_setup(**kwargs)

        self._m1_happen = 0
        self._m2_happen = 0
        self._m3_happen = 0

        self._n_pump = 2.0

        # Dynamic storage variables (replace getDynamicArrayValueLastTimestep / setDynamicArrayValueThisIteration)
        self._T_of_state = self.T_of.v          # Temperature of the overflow tank
        self._T_ev_state = self.T_ev.v          # Temperature of the expansion tank
        self._LV_204_05_prev = 0.0              # Control valve HX
        self._LV_204_09B_prev = 0.0             # Control valve OF to Exp
        self._HV_204_09A_prev = 0.0             # Control valve Exp to OF
        self._LV_205_14_prev = 0.0              # Control valve HX bypass
        self._P_bot_of = self.P_of.v            # Pressure at bottom of the overflow tank
        self._m_count_prev = 0.0

        # Per-mode flow rate state
        self._m1_4012 = np.zeros(3)
        self._m2_4012 = np.zeros(3)
        self._m3_4012 = np.zeros(4)

        # Tank mass state
        self._mass_of = 0.0
        self._mass_ev = 0.0

        # Piping geometry arrays (read from geometry file at start time)
        self._D = np.zeros(9)   # Pipe diameters [m]
        self._L = np.zeros(9)   # Pipe lengths [m]

        # Valve CV values (computed from CV_data lookup)
        self._CV_HV_204_09 = 0.0
        self._CV_LV_204_09 = 0.0
        self._CV_LV_204_05 = 0.0
        self._CV_LV_205_14 = 0.0

        # Gravity
        self._g = 9.81

    def calculate(self):
        super().calculate()

        # Do All of the First Timestep Manipulations Here - There Are No Iterations at the Initial Time
        if self.model.is_first_step:

            # Check if only using expansion tanks
            if self.only_Exp.v == 1.0:
                # Compute mass in expansion tanks
                mass_ev = self.total_mass.v - self.m_counter.v
                # Compute level in expansion system
                self.level_ev.v = mass_ev / Inc.density(self.Fluid_ID.v, T=self.T_ev.v, P=self.P_ev.v) / ((self.pi * (self.D_ev.v / 2) ** 2 * self.H_ev.v) * self.N_ev.v)

                # Set overflow variables to zero
                self._mass_of = 0.0
                self.level_of.v = 0.0

            else:
                # Compute mass in overflow system
                self._mass_of = ((self.pi * (self.D_of.v / 2) ** 2) * self.H_of_init.v * self.H_of.v * self.N_of.v) * Inc.density(self.Fluid_ID.v, T=self.T_of.v, P=self.P_of.v)
                self.level_of.v = self.H_of_init.v

                # Compute mass in expansion system
                mass_ev = self.total_mass.v - self._mass_of - self.m_counter.v

                # Compute level in expansion system
                self.level_ev.v = mass_ev / Inc.density(self.Fluid_ID.v, T=self.T_ev.v, P=self.P_ev.v) / ((self.pi * (self.D_ev.v / 2) ** 2 * self.H_ev.v) * self.N_ev.v)

                # Read in Piping Geometry
                # Get specified geometries
                # TODO-NEEDS CONVERSION REVIEW: getLabel(CurrentUnit, 1) maps to a geometry file path label.
                # No direct esclab equivalent; GeomFile must be supplied as an instance attribute or parameter.
                GeomFile = getattr(self, '_geom_file', None)
                if GeomFile is not None:
                    with open(GeomFile) as f_geom:
                        f_geom.readline()
                        f_geom.readline()
                        # Read in diameters
                        for cc in range(9):
                            f_geom.readline()
                            self._D[cc] = float(f_geom.readline().strip())
                        # Read in lengths
                        f_geom.readline()
                        for cc in range(9):
                            f_geom.readline()
                            self._L[cc] = float(f_geom.readline().strip())

                # Compute Cv values for valves
                self._CV_HV_204_09 = CV_data(2, self._D[4], self.HV_204_09A.v)
                self._CV_LV_204_09 = CV_data(2, self._D[2], self.LV_204_09B.v)
                self._CV_LV_204_05 = CV_data(2, self._D[6], self.LV_204_05.v) 
                self._CV_LV_205_14 = CV_data(2, self._D[7], self.LV_205_14.v) 

            self.m_dot_in_out.v = self.m_dot_in.v
            self.T_in_out.v = self.T_in.v
            self.P_bot_ev.v = self.P_ev.v + Inc.density(self.Fluid_ID.v, T=self.T_ev.v, P=self.P_of.v) * self._g * self.level_ev.v * self.H_ev.v
            self.net_of_flow.v = 0.0  # m_in_of - m_out_of (uninitialized at start)
            self.m_dot_exp.v = 0.0    # m_dot_exp uninitialized at start
            self.m_dot_recirc.v = 0.0
            self.m_dot_chiller.v = 0.0
            self.T_of_out.v = self.T_of.v
            self.T_ev_out.v = self.T_ev.v

            # Store initial dynamic values
            self._T_of_state = self.T_of.v              # Temperature of the overflow tank
            self._T_ev_state = self.T_ev.v              # Temperature of the expansion tank
            self._LV_204_05_prev = self.LV_204_05.v     # Control valve HX
            self._LV_204_09B_prev = self.LV_204_09B.v   # Control valve OF to Exp
            self._HV_204_09A_prev = self.HV_204_09A.v   # Control valve Exp to OF
            self._LV_205_14_prev = self.LV_205_14.v     # Control valve HX bypass
            # Pressure at bottom of overflow
            self._P_bot_of = self.P_of.v + Inc.density(self.Fluid_ID.v, T=self.T_of.v, P=self.P_of.v) * self._g * self.level_of.v * self.H_of.v
            self._m_count_prev = self.m_counter.v
            return

        # Load state at the start of each new timestep
        if self.model.timestep_iteration == 0:
            T_of = self._T_of_state                  # Temperature of HTF in overflow tanks
            T_ev = self._T_ev_state                  # Temperature of HTF in expansion tanks
            LV_204_05_prev = self._LV_204_05_prev    # Control valve HX
            LV_204_09B_prev = self._LV_204_09B_prev  # Control valve OF to Exp
            HV_204_09A_prev = self._HV_204_09A_prev  # Control valve Exp to OF
            LV_205_14_prev = self._LV_205_14_prev    # Control valve HX bypass
            P_bot_of = self._P_bot_of                # Pressure at bottom of the overflow tank
            P_bot_ev = self.P_bot_ev.v               # Pressure at bottom of the expansion tank
            m_count_prev = self._m_count_prev

        # Check if Expansion Tanks Only
        if self.only_Exp.v == 1.0:

            # Update Dynamic Storage and temperatures at the end of each timestep
            if self.model.is_converged:
                # Reload in case they somehow got overwritten by another type
                T_ev = self._T_ev_state
                m_count_prev = self._m_count_prev

                # Compute mass flow into or out of expansion tank
                m_dot_SF_to_exp = 0.0
                if self.model.time > self.model.settings.timestep:
                    m_dot_SF_to_exp = (m_count_prev - self.m_counter.v) / self.model.settings.timestep

                # Update tank temperature
                if m_dot_SF_to_exp > 0:
                    h_plant = Inc.enthalpy(self.Fluid_ID.v, T=self.T_in.v)
                else:
                    h_plant = Inc.enthalpy(self.Fluid_ID.v, T=T_ev)
                h = Inc.enthalpy(self.Fluid_ID.v, T=T_ev)
                c = Inc.specheat(self.Fluid_ID.v, T=T_ev, P=0.0)
                # TODO-NEEDS CONVERSION REVIEW: mass_ev is a local variable not yet computed at this point;
                # it is computed below, but the temperature update uses the previous mass_ev.
                mass_ev = self.total_mass.v - self.m_counter.v
                dT_dt = 1 / (mass_ev * c) * (-h * m_dot_SF_to_exp + m_dot_SF_to_exp * h_plant)
                T_ev = T_ev + dT_dt * self.model.settings.timestep

                self._T_ev_state = T_ev
                self._m_count_prev = self.m_counter.v

            # Compute mass in expansion tanks
            mass_ev = self.total_mass.v - self.m_counter.v

            # Compute level in expansion tanks
            self.level_ev.v = mass_ev / Inc.density(self.Fluid_ID.v, T=T_ev, P=self.P_ev.v) / ((self.pi * (self.D_ev.v / 2) ** 2 * self.H_ev.v) * self.N_ev.v)

            # Send outputs
            self.m_dot_in_out.v = self.m_dot_in.v
            self.T_in_out.v = self.T_in.v
            self.P_bot_ev.v = self.P_ev.v + Inc.density(self.Fluid_ID.v, T=T_ev, P=self.P_of.v) * self._g * self.level_ev.v * self.H_ev.v
            self.net_of_flow.v = 0.0
            self.m_dot_exp.v = 0.0
            self.m_dot_recirc.v = 0.0
            self.m_dot_chiller.v = 0.0
            self.T_of_out.v = self.T_of.v
            self.T_ev_out.v = T_ev

        else:
            # Update tank temperatures, levels, and flowrates
            if self.model.is_converged:
                # Load in relevant variables
                P_bot_of = self._P_bot_of
                P_bot_ev = self.P_bot_ev.v

                level_ev = self.level_ev.v
                level_of = self.level_of.v

                T_of = self.T_of_out.v
                T_ev = self.T_ev_out.v

                # Reload m_count_prev from instance state (local var only set at timestep_iteration==0)
                m_count_prev = self._m_count_prev

                ############################################################################
                ######  1) COMPUTE FLOW RATES TO AND FROM OVERFLOW/EXPANSION  ##############
                ########  note: m_dot_exp = (+) for mass flow to expansion tank
                ########        m_dot_of = mass flow entering of tank - mass flow leaving of tank
                ########        m_dot_recirc = mass flow recirculating

                m_dot_exp = 0.0
                T_to_exp = 0.0
                m_in_of = 0.0
                m_out_of = 0.0
                T_to_of = 0.0
                m_dot_recirc = 0.0
                m_dot_chiller = 0.0

                alarm_high_high_pressure = 0.0

                ### MODE 1: Moving HTF from the overflow tanks to the expansion tanks
                if self.LV_204_09B.v > 0:
                    # Initialize mass flow guesses if this mode has never happened before
                    if self._m1_happen == 0:
                        self._m1_4012[0] = 40.0
                        self._m1_4012[1] = 20.0
                        self._m1_4012[2] = 20.0
                        self._m1_happen = 1

                    # Iterate until tolerance is achieved
                    tol_it = 0.1
                    error_tot = 100.0

                    while error_tot > tol_it:
                        # Compute P_1
                        dp1 = PressureDrop(self.Fluid_ID.v, self._m1_4012[0], T_of, 1.0, self._D[0],
                                           self.roughness.v, self._L[0], 0.0, 0.0, 4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)

                        dp_pump = self._g * (Inc.density(self.Fluid_ID.v, T=T_of, P=self.P_of.v) * self.PC_a.v + self.PC_b.v * self._m1_4012[0] / self._n_pump + self.PC_c.v * self._m1_4012[0] ** 2 / Inc.density(self.Fluid_ID.v, T=T_of, P=self.P_of.v) / self._n_pump ** 2)

                        P_1 = P_bot_of - math.copysign(dp1, self._m1_4012[0]) + dp_pump

                        # Compute P_34
                        # Compute Volumetric Flow Rate
                        rho_fluid = Inc.density(self.Fluid_ID.v, T=T_of, P=0.0)
                        Q = self._m1_4012[2] / rho_fluid * convert('m^3/s', 'gpm')  # Convert flowrate to gpm
                        # Compute specific gravity of fluid
                        SG = rho_fluid / 1000.0
                        dp3 = PressureDrop(self.Fluid_ID.v, self._m1_4012[2], T_of, 1.0, self._D[2],
                                           self.roughness.v, self._L[2], 0.0, 0.0, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0) \
                                           + SG * Q ** 2 / (self._CV_LV_204_09 ** 2) * convert('psi','Pa')

                        dp34 = PressureDrop(self.Fluid_ID.v, self._m1_4012[2], T_of, 1.0, self._D[3],
                                            self.roughness.v, self._L[3], 0.0, 0.0, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)

                        P_34 = P_1 - math.copysign(dp3, self._m1_4012[2]) - math.copysign(dp34, self._m1_4012[2])

                        # Compute P_8
                        # Compute Volumetric Flow Rate
                        Q = self._m1_4012[1] / rho_fluid * convert('m^3/s', 'gpm')  # Convert flowrate to gpm
                        
                        dp_orifice = SG * Q ** 2 / (self.CV_orifice.v ** 2) * convert('psi','Pa')

                        dp2 = PressureDrop(self.Fluid_ID.v, self._m1_4012[1], T_of, 1.0, self._D[1],
                                           self.roughness.v, self._L[1], 0.0, 0.0, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)

                        dp8 = PressureDrop(self.Fluid_ID.v, self._m1_4012[1], T_of, 1.0, self._D[4],
                                           self.roughness.v, self._L[5] + self._L[7] + self._L[8], 0.0, 0.0, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0) \
                                            + SG * Q ** 2 / (self._CV_LV_205_14 ** 2) * convert('psi','Pa')

                        dp_elevation = self.H_of.v * rho_fluid * self._g

                        P_8 = P_1 - math.copysign(dp2, self._m1_4012[1]) - math.copysign(dp_orifice, self._m1_4012[1]) - math.copysign(dp8, self._m1_4012[1]) - dp_elevation

                        # Compute Error Values
                        error_tot = abs(P_8 - self.P_of.v) + abs(P_34 - P_bot_ev)

                        # Compute K values
                        K1 = math.copysign(dp1 / self._m1_4012[0] ** 2, self._m1_4012[0])

                        K3 = math.copysign((dp34 + dp3) / self._m1_4012[2] ** 2, self._m1_4012[2])

                        K8 = math.copysign((dp_orifice + dp2 + dp8) / self._m1_4012[1] ** 2, self._m1_4012[1])

                        # Update matrix
                        A1 = np.zeros((3, 3))

                        # Continuity
                        A1[0, 0] = -1.0
                        A1[0, 1] = 1.0
                        A1[0, 2] = 1.0
                        b1 = np.zeros(3)
                        b1[0] = 0.0

                        # Pressure Constraint: P_34 = P_bot_ev
                        A1[1, 0] = self._g * (self.PC_b.v / self._n_pump + self._m1_4012[0] * self.PC_c.v / (Inc.density(self.Fluid_ID.v, T=T_of, P=self.P_of.v) * self._n_pump ** 2)) - self._m1_4012[0] * K1
                        A1[1, 2] = -self._m1_4012[2] * K3
                        b1[1] = P_bot_ev - P_bot_of - self._g * rho_fluid * self.PC_a.v

                        # Pressure Constraint: P_8 = P_of
                        A1[2, 0] = self._g * (self.PC_b.v / self._n_pump + self._m1_4012[0] * self.PC_c.v / (rho_fluid * self._n_pump ** 2)) - self._m1_4012[0] * K1
                        A1[2, 1] = -self._m1_4012[1] * K8
                        b1[2] = self.P_of.v + dp_elevation - P_bot_of - self._g * rho_fluid * self.PC_a.v

                        # Invert matrix to obtain new flow rates
                        # matrixinv available from sergio_scripts; np.linalg.inv used here as equivalent
                        A1inv = np.linalg.inv(A1)
                        mdots_new = A1inv @ b1

                        ## Update flow rates
                        # If all valves are closed
                        if mdots_new[0] < 0.001 and mdots_new[0] > 0:
                            m_hold = math.log10(mdots_new[0]) * 0.5 + math.log10(self._m1_4012[0]) * (1 - 0.5)
                            self._m1_4012[0] = 10 ** m_hold
                        else:
                            self._m1_4012[0] = self._m1_4012[0] * 0.5 + mdots_new[0] * 0.5

                        # If all valves are closed
                        if mdots_new[2] < 0.001 and mdots_new[2] > 0:
                            m_hold = math.log10(mdots_new[2]) * 0.5 + math.log10(self._m1_4012[2]) * (1 - 0.5)
                            self._m1_4012[2] = 10 ** m_hold
                        else:
                            self._m1_4012[2] = self._m1_4012[2] * 0.5 + mdots_new[2] * 0.5

                        # If LV-204-09 is barely cracked open
                        if mdots_new[1] < 0.001 and mdots_new[1] > 0:
                            m_hold = math.log10(mdots_new[1]) * 0.5 + math.log10(self._m1_4012[1]) * (1 - 0.5)
                            self._m1_4012[1] = 10 ** m_hold
                        else:
                            self._m1_4012[1] = self._m1_4012[1] * 0.5 + mdots_new[1] * 0.5

                    # Set mass flow rates
                    m_out_of = self._m1_4012[0]
                    m_in_of = self._m1_4012[1]
                    m_dot_recirc = self._m1_4012[1]
                    m_dot_exp = self._m1_4012[2]

                    # Set temperatures
                    T_to_exp = T_of
                    T_to_of = T_of

                ### MODE 2: Moving HTF from expansion tanks to overflow tanks ###
                elif self.HV_204_09A.v > 0:
                    # Initialize mass flow guesses if this mode has never happened before
                    if self._m2_happen == 0:
                        self._m2_4012[0] = 32.0
                        self._m2_4012[1] = 100.0
                        self._m2_4012[2] = self._m2_4012[0] + self._m2_4012[1]
                        self._m2_happen = 1

                    # Iterate until tolerance is achieved
                    tol_it = 0.1
                    error_tot = 100.0

                    while error_tot > tol_it:

                        # Compute Pressure drops for fluid from expansion to overflow
                        dp34 = PressureDrop(self.Fluid_ID.v, self._m2_4012[1], T_ev, 1.0, self._D[3],
                                            self.roughness.v, self._L[3], 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)

                        # Compute Volumetric Flow Rate
                        rho_fluid = Inc.density(self.Fluid_ID.v, T=T_ev, P=0.0)
                        Q = self._m2_4012[1] / rho_fluid * convert('m^3/s', 'gpm')  # Convert flowrate to gpm
                        
                        # Compute specific gravity of fluid
                        SG = rho_fluid / 1000.0
                        dp4 = PressureDrop(self.Fluid_ID.v, self._m2_4012[1], T_ev, 1.0, self._D[4],
                                           self.roughness.v, self._L[4], 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0) \
                                            + SG * Q ** 2 / (self._CV_HV_204_09 ** 2) * convert('psi','Pa')

                        P_4 = P_bot_ev - math.copysign(dp34, self._m2_4012[1]) - math.copysign(dp4, self._m2_4012[1])

                        # Compute Volumetric Flow Rate
                        rho_fluid = Inc.density(self.Fluid_ID.v, T=(T_ev * self._m2_4012[1] + T_of * self._m2_4012[0]) / (self._m2_4012[1] + self._m2_4012[0]), P=0.0)
                        Q = self._m2_4012[2] / rho_fluid * convert('m^3/s', 'gpm')  # Convert flowrate to gpm
                        
                        # Compute specific gravity of fluid
                        SG = rho_fluid / 1000.0
                        T_mix = (T_ev * self._m2_4012[1] + T_of * self._m2_4012[0]) / (self._m2_4012[1] + self._m2_4012[0])
                        dp8 = PressureDrop(self.Fluid_ID.v, self._m2_4012[2], T_mix, 1.0, self._D[4],
                                           self.roughness.v, self._L[5] + self._L[7] + self._L[8], 0.0, 0.0, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0) \
                                            + SG * Q ** 2 / (self._CV_LV_205_14 ** 2) * convert('psi','Pa')

                        dp_elevation = self.H_of.v * Inc.density(self.Fluid_ID.v, T=T_mix, P=self.P_of.v) * self._g

                        P_8 = P_bot_ev - math.copysign(dp34, self._m2_4012[1]) - math.copysign(dp4, self._m2_4012[1]) - math.copysign(dp8, self._m2_4012[2]) - dp_elevation

                        # Compute pressure drops from overflow to overflow
                        dp1 = PressureDrop(self.Fluid_ID.v, self._m2_4012[0], T_of, 1.0, self._D[0],
                                           self.roughness.v, self._L[0] + self._L[1], 0.0, 0.0, 4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)

                        rho_of = Inc.density(self.Fluid_ID.v, T=T_of, P=self.P_of.v)
                        dp_pump = self._g * (rho_of * self.PC_a.v + self.PC_b.v * self._m2_4012[0] / self._n_pump + self.PC_c.v * self._m2_4012[0] ** 2 / rho_of / self._n_pump ** 2)

                        # Compute Volumetric Flow Rate
                        rho_fluid = Inc.density(self.Fluid_ID.v, T=T_of, P=0.0)
                        Q = self._m2_4012[0] / rho_fluid * convert('m^3/s', 'gpm')  # Convert flowrate to gpm
                        
                        # Compute specific gravity of fluid
                        SG = rho_fluid / 1000.0

                        dp_orifice = SG * Q ** 2 / (self.CV_orifice.v ** 2) * convert('psi','Pa')

                        dp1 = dp1 + dp_orifice

                        P_2 = P_bot_of - math.copysign(dp1, self._m2_4012[0]) + dp_pump

                        # Compute Error Values
                        error_tot = abs(P_8 - self.P_of.v) + abs(P_2 - P_4)

                        # Compute K values
                        K1 = math.copysign(dp1 / self._m2_4012[0] ** 2, self._m2_4012[0])

                        K4 = math.copysign((dp34 + dp4) / self._m2_4012[1] ** 2, self._m2_4012[1])

                        K8 = dp8 / self._m2_4012[2] ** 2

                        # Update matrix
                        A2 = np.zeros((3, 3))

                        # Continuity
                        A2[0, 0] = 1.0
                        A2[0, 1] = 1.0
                        A2[0, 2] = -1.0
                        b2 = np.zeros(3)
                        b2[0] = 0.0

                        # Pressure Constraint: P_8 = P_of
                        A2[1, 1] = -self._m2_4012[1] * K4
                        A2[1, 2] = -self._m2_4012[2] * K8
                        b2[1] = self.P_of.v + dp_elevation - P_bot_ev

                        # Pressure Constraint: P_2 = P_4
                        A2[2, 0] = self._g * (self.PC_b.v / self._n_pump + self._m2_4012[0] * self.PC_c.v / (Inc.density(self.Fluid_ID.v, T=T_of, P=self.P_of.v) * self._n_pump ** 2)) - self._m2_4012[0] * K1
                        A2[2, 1] = self._m2_4012[1] * K4
                        b2[2] = P_bot_ev - P_bot_of - self._g * Inc.density(self.Fluid_ID.v, T=T_of, P=self.P_of.v) * self.PC_a.v

                        # Invert matrix to obtain new flow rates
                        # matrixinv available from sergio_scripts; np.linalg.inv used here as equivalent
                        A2inv = np.linalg.inv(A2)
                        mdots_new = A2inv @ b2

                        ## Update flow rates
                        # If all valves are closed
                        if mdots_new[0] < 0.001 and mdots_new[0] > 0:
                            m_hold = math.log10(mdots_new[0]) * 0.5 + math.log10(self._m2_4012[0]) * (1 - 0.5)
                            self._m2_4012[0] = 10 ** m_hold
                        else:
                            self._m2_4012[0] = self._m2_4012[0] * 0.5 + mdots_new[0] * 0.5

                        # If all valves are closed
                        if mdots_new[2] < 0.001 and mdots_new[2] > 0:
                            m_hold = math.log10(mdots_new[2]) * 0.5 + math.log10(self._m2_4012[2]) * (1 - 0.5)
                            self._m2_4012[2] = 10 ** m_hold
                        else:
                            self._m2_4012[2] = self._m2_4012[2] * 0.5 + mdots_new[2] * 0.5

                        # If HV-204-09 is barely cracked open
                        if mdots_new[1] < 0.001 and mdots_new[1] > 0:
                            m_hold = math.log10(mdots_new[1]) * 0.5 + math.log10(self._m2_4012[1]) * (1 - 0.5)
                            self._m2_4012[1] = 10 ** m_hold
                        else:
                            self._m2_4012[1] = self._m2_4012[1] * 0.5 + mdots_new[1] * 0.5

                    # Set mass flow rates
                    m_out_of = self._m2_4012[0]
                    m_in_of = self._m2_4012[2]
                    m_dot_recirc = self._m2_4012[0]
                    m_dot_exp = -self._m2_4012[1]

                    # Set temperatures
                    T_to_exp = 0.0
                    T_to_of = T_mix

                ### MODE 3: HTF Re-circulation/Cooling
                else:
                    # Initialize mass flow guesses if this mode has never happened before
                    if self._m3_happen == 0:
                        self._m3_4012[0] = 100.0
                        self._m3_4012[1] = 50.0
                        self._m3_4012[2] = 50.0
                        self._m3_4012[3] = 100.0
                        self._m3_happen = 1

                    # Iterate until tolerance is achieved
                    tol_it = 0.1
                    error_tot = 100.0

                    while error_tot > tol_it:
                        ## Compute P_5
                        # Compute Volumetric Flow Rate
                        rho_fluid = Inc.density(self.Fluid_ID.v, T=T_of, P=0.0)
                        Q = self._m3_4012[0] / rho_fluid * convert('m^3/s', 'gpm')  # Convert flowrate to gpm
                        
                        # Compute specific gravity of fluid
                        SG = rho_fluid / 1000.0

                        dp1 = PressureDrop(self.Fluid_ID.v, self._m3_4012[0], T_of, 1.0, self._D[0],
                                           self.roughness.v, self._L[0] + self._L[1], 0.0, 0.0, 4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)

                        dp_pump = self._g * (rho_fluid * self.PC_a.v + self.PC_b.v * self._m3_4012[0] / self._n_pump + self.PC_c.v * self._m3_4012[0] ** 2 / rho_fluid / self._n_pump ** 2)

                        dp_orifice = SG * Q ** 2 / (self.CV_orifice.v ** 2) * convert('psi','Pa')

                        dp5 = PressureDrop(self.Fluid_ID.v, self._m3_4012[0], T_of, 1.0, self._D[5],
                                           self.roughness.v, self._L[5], 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)

                        P_2 = P_bot_of - dp1 + dp_pump - dp_orifice - dp5

                        # Compute P_7
                        Q = self._m3_4012[2] / rho_fluid * convert('m^3/s', 'gpm')  # Convert flowrate to gpm
                        dp7 = PressureDrop(self.Fluid_ID.v, self._m3_4012[2], T_of, 1.0, self._D[7],
                                           self.roughness.v, self._L[7], 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0) + SG * Q ** 2 / (self._CV_LV_205_14 ** 2) * convert('psi','Pa')

                        P_7 = P_2 - dp7

                        ## Compute P_6
                        # Compute Volumetric Flow Rate
                        rho_fluid = Inc.density(self.Fluid_ID.v, T=self.T_chill_out.v, P=0.0)
                        Q = self._m3_4012[1] / rho_fluid * convert('m^3/s', 'gpm')  # Convert flowrate to gpm
                        
                        # Compute specific gravity of fluid
                        SG = rho_fluid / 1000.0

                        dp6_1 = PressureDrop(self.Fluid_ID.v, self._m3_4012[1], T_of, 1.0, self._D[6],
                                             self.roughness.v, self._L[6] / 2, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)

                        dp6_2 = PressureDrop(self.Fluid_ID.v, self._m3_4012[1], self.T_chill_out.v, 1.0, self._D[6],
                                             self.roughness.v, self._L[6] / 2, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0) + SG * Q ** 2 / (self._CV_LV_204_05 ** 2) * convert('psi','Pa')

                        K4 = PressureDrop(self.Fluid_ID.v, self._m3_4012[2], self.T_chill_out.v, 1.0, self._D[6],
                                          self.roughness.v, self._L[6] / 2, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)

                        P_6 = P_2 - dp6_1 - dp6_2

                        ## Compute P_8
                        T_8 = (self.T_chill_out.v * self._m3_4012[1] + T_of * self._m3_4012[2]) / (self._m3_4012[0])
                        # Compute Volumetric Flow Rate
                        rho_fluid = Inc.density(self.Fluid_ID.v, T=T_8, P=0.0)

                        dp8 = PressureDrop(self.Fluid_ID.v, self._m3_4012[3], T_8, 1.0, self._D[4],
                                           self.roughness.v, self._L[8], 0.0, 0.0, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)

                        dp_elevation = self.H_of.v * rho_fluid * self._g

                        P_8 = P_6 - dp8 - dp_elevation

                        # Compute Error Values
                        error_tot = abs(P_8 - self.P_of.v) + abs(P_6 - P_7)

                        # Compute K values
                        K1 = math.copysign((dp1 + dp_orifice + dp5) / self._m3_4012[0] ** 2, self._m3_4012[0])

                        K6 = math.copysign((dp6_1 + dp6_2) / self._m3_4012[1] ** 2, self._m3_4012[1])

                        K7 = math.copysign(dp7 / self._m3_4012[2] ** 2, self._m3_4012[2])

                        K8 = math.copysign(dp8 / self._m3_4012[3] ** 2, self._m3_4012[3])

                        # Update matrix
                        A3 = np.zeros((4, 4))

                        # Continuity
                        A3[0, 0] = 1.0
                        A3[0, 1] = -1.0
                        A3[0, 2] = -1.0

                        A3[1, 1] = 1.0
                        A3[1, 2] = 1.0
                        A3[1, 3] = -1.0

                        b3 = np.zeros(4)
                        b3[0] = 0.0
                        b3[1] = 0.0

                        # Pressure Constraint: P_7 = P_8
                        A3[2, 1] = self._m3_4012[1] * K6
                        A3[2, 2] = -self._m3_4012[2] * K7
                        b3[2] = 0.0

                        # Pressure Constraint: P_8 = P_of
                        rho_fluid = Inc.density(self.Fluid_ID.v, T=T_of, P=0.0)
                        A3[3, 0] = self._g * (self.PC_b.v / self._n_pump + self._m3_4012[0] * self.PC_c.v / (rho_fluid * self._n_pump ** 2)) - self._m3_4012[0] * K1
                        A3[3, 2] = -self._m3_4012[2] * K7
                        A3[3, 3] = -self._m3_4012[3] * K8
                        b3[3] = self.P_of.v + dp_elevation - P_bot_of - self._g * rho_fluid * self.PC_a.v

                        # Invert matrix to obtain new flow rates
                        # matrixinv available from sergio_scripts; np.linalg.inv used here as equivalent
                        A3inv = np.linalg.inv(A3)
                        mdots_new = A3inv @ b3

                        ## Update flow rates
                        # If all valves are closed
                        if mdots_new[0] < 0.001 and mdots_new[0] > 0:
                            m_hold = math.log10(mdots_new[0]) * 0.5 + math.log10(self._m3_4012[0]) * (1 - 0.5)
                            self._m3_4012[0] = 10 ** m_hold
                        else:
                            self._m3_4012[0] = self._m3_4012[0] * 0.5 + mdots_new[0] * 0.5

                        # If all valves are closed
                        if mdots_new[3] < 0.001 and mdots_new[3] > 0:
                            m_hold = math.log10(mdots_new[3]) * 0.5 + math.log10(self._m3_4012[3]) * (1 - 0.5)
                            self._m3_4012[3] = 10 ** m_hold
                        else:
                            self._m3_4012[3] = self._m3_4012[3] * 0.5 + mdots_new[3] * 0.5

                        # If Chiller line is closed
                        if mdots_new[1] < 0.001 and mdots_new[1] > 0:
                            m_hold = math.log10(mdots_new[1]) * 0.5 + math.log10(self._m3_4012[1]) * (1 - 0.5)
                            self._m3_4012[1] = 10 ** m_hold
                        else:
                            self._m3_4012[1] = self._m3_4012[1] * 0.5 + mdots_new[1] * 0.5

                        # If bypass line is closed
                        if mdots_new[2] < 0.001 and mdots_new[2] > 0:
                            m_hold = math.log10(mdots_new[2]) * 0.5 + math.log10(self._m3_4012[2]) * (1 - 0.5)
                            self._m3_4012[2] = 10 ** m_hold
                        else:
                            self._m3_4012[2] = self._m3_4012[2] * 0.5 + mdots_new[2] * 0.5

                    # Set mass flow rates
                    m_out_of = self._m3_4012[0]
                    m_in_of = self._m3_4012[3]
                    m_dot_chiller = self._m3_4012[1]
                    m_dot_recirc = self._m3_4012[0]
                    m_dot_exp = 0.0

                    # Set temperatures
                    T_to_exp = 0.0
                    T_to_of = T_8

                ############################################################################
                ######  2) Update tank levels  #############################################
                m_dot_SF_to_exp = 0.0
                if self.model.timeTime > self.model.settings.timestep:
                    m_dot_SF_to_exp = (m_count_prev - self.m_counter.v) / self.model.settings.timestep

                # Compute mass and level in overflow system
                self._mass_of = self._mass_of + (m_in_of - m_out_of) * self.model.settings.timestep
                self.level_of.v = self._mass_of / Inc.density(self.Fluid_ID.v, T=T_of, P=self.P_of.v) / ((self.pi * (self.D_of.v / 2) ** 2 * self.H_of.v) * self.N_of.v)

                # Compute mass and level in expansion system
                mass_ev = self.total_mass.v - self._mass_of - self.m_counter.v
                self.level_ev.v = mass_ev / Inc.density(self.Fluid_ID.v, T=T_ev, P=self.P_ev.v) / ((self.pi * (self.D_ev.v / 2) ** 2 * self.H_ev.v) * self.N_ev.v)

                ############################################################################
                ######  2) Update tank temperatures  #######################################
                # Overflow tank (Eulers Method)
                c = Inc.specheat(self.Fluid_ID.v, T=T_of, P=0.0)
                h = Inc.enthalpy(self.Fluid_ID.v, T=T_of)
                h_in = Inc.enthalpy(self.Fluid_ID.v, T=T_to_of)
                h_out = Inc.enthalpy(self.Fluid_ID.v, T=T_of)
                dT_dt = 1 / (self._mass_of * c) * (m_in_of * (h_in - h) - m_out_of * (h - h_out))
                T_of = T_of + dT_dt * self.model.settings.timestep

                # Expansion Tank (Eulers Method)
                if m_dot_exp > 0:
                    m_in_ev = m_dot_exp
                    h_in = Inc.enthalpy(self.Fluid_ID.v, T=T_to_exp)
                    m_out_ev = 0.0
                    h_out = 0.0
                else:
                    m_out_ev = -m_dot_exp
                    h_out = Inc.enthalpy(self.Fluid_ID.v, T=T_ev)
                    m_in_ev = 0.0
                    h_in = 0.0
                c = Inc.specheat(self.Fluid_ID.v, T=T_ev, P=0.0)
                h_plant = Inc.enthalpy(self.Fluid_ID.v, T=self.T_in.v)
                h = Inc.enthalpy(self.Fluid_ID.v, T=T_ev)
                dT_dt = 1 / (mass_ev * c) * (-h * (m_dot_SF_to_exp + m_in_ev - m_out_ev) + m_dot_SF_to_exp * h_plant + m_in_ev * h_in - m_out_ev * h_out)
                T_ev = T_ev + dT_dt * self.model.settings.timestep

                ############################################################################
                ######  3) Update pressure values at bottom of tanks
                P_bot_of = self.level_of.v * self.H_of.v * Inc.density(self.Fluid_ID.v, T=T_of, P=self.P_of.v) * self._g + self.P_of.v
                P_bot_ev = self.level_ev.v * self.H_ev.v * Inc.density(self.Fluid_ID.v, T=T_ev, P=self.P_ev.v) * self._g + self.P_ev.v

                if P_bot_ev > 1.172e6:
                    alarm_high_high_pressure = 1.0
                else:
                    alarm_high_high_pressure = 0.0

                # Set outputs
                self.m_dot_in_out.v = self.m_dot_in.v
                self.T_in_out.v = self.T_in.v
                self.P_bot_ev.v = P_bot_ev
                self.net_of_flow.v = m_in_of - m_out_of
                self.m_dot_exp.v = m_dot_exp
                self.m_dot_recirc.v = m_dot_recirc
                self.m_dot_chiller.v = m_dot_chiller
                self.T_of_out.v = T_of
                self.T_ev_out.v = T_ev
                self.alarm_high_high_pressure.v = alarm_high_high_pressure

                # Update dynamic storage
                self._T_of_state = T_of              # Temperature of the overflow tank
                self._T_ev_state = T_ev              # Temperature of the expansion tank
                self._P_bot_of = P_bot_of            # Pressure at bottom of overflow
                self._m_count_prev = self.m_counter.v

                # call find_66_outputs(Time, CurrentUnit, getTimestepIteration(), end - start, [kappa(A1, A1inv, 3), kappa(A2, A2inv, 3), kappa(A3, A3inv, 4)])
                # call find_66_iter(Time, CurrentUnit, getTimestepIteration())

                return

            # Check the Valve Positions for Problems
            # TODO-NEEDS CONVERSION REVIEW: FoundBadInput / fatal error checking; raising ValueError here
            if self.HV_204_09A.v > 0 and self.LV_204_09B.v > 0:
                raise ValueError("Cant have HV_204_09A and LV_204_09B open at the same time")
            if self.HV_204_09A.v > 0 and self.LV_204_05.v > 0:
                raise ValueError("Cant have HV_204_09A and LV_204_05 open at the same time")
            if self.LV_204_09B.v > 0 and self.LV_204_05.v > 0:
                raise ValueError("Cant have LV_204_09B and LV_204_05 open at the same time")
            if self.HV_204_09A.v > 0 and self.LV_205_14.v == 0:
                raise ValueError("Cant have LV_205_14 closed with HV_204_09A open at the same time")

            if self.HV_204_09A.v == 0.0:
                self.m_counter.v = 1.0

            # Update Valve Positions and CV Values
            if self.model.timestep_iteration == 0:

                # CV Expansion to OF
                if self.HV_204_09A.v != HV_204_09A_prev:
                    if self.HV_204_09A.v > HV_204_09A_prev:
                        HV_204_09A_cur = min(self.HV_204_09A.v, HV_204_09A_prev + self.valve_speed.v / 90.0 * self.model.settings.timestep)
                    else:
                        HV_204_09A_cur = max(self.HV_204_09A.v, HV_204_09A_prev - self.valve_speed.v / 90.0 * self.model.settings.timestep)
                    self._CV_HV_204_09 = CV_data(1, self._D[4], HV_204_09A_cur)
                    self._HV_204_09A_prev = HV_204_09A_cur

                # CV OF to Expansion
                if self.LV_204_09B.v != LV_204_09B_prev:
                    if self.LV_204_09B.v > LV_204_09B_prev:
                        LV_204_09B_cur = min(self.LV_204_09B.v, LV_204_09B_prev + self.valve_speed.v / 90.0 * self.model.settings.timestep)
                    else:
                        LV_204_09B_cur = max(self.LV_204_09B.v, LV_204_09B_prev - self.valve_speed.v / 90.0 * self.model.settings.timestep)
                    self._CV_LV_204_09 = CV_data(1, self._D[2], LV_204_09B_cur)
                    self._LV_204_09B_prev = LV_204_09B_cur

                # CV to HX
                if self.LV_204_05.v != LV_204_05_prev:
                    if self.LV_204_05.v > LV_204_05_prev:
                        LV_204_05_cur = min(self.LV_204_05.v, LV_204_05_prev + self.valve_speed.v / 90.0 * self.model.settings.timestep)
                    else:
                        LV_204_05_cur = max(self.LV_204_05.v, LV_204_05_prev - self.valve_speed.v / 90.0 * self.model.settings.timestep)
                    self._CV_LV_204_05 = CV_data(1, self._D[6], LV_204_05_cur)
                    self._LV_204_05_prev = LV_204_05_cur

                # CV to HX bypass
                if self.LV_205_14.v != LV_205_14_prev:
                    if self.LV_205_14.v > LV_205_14_prev:
                        LV_205_14_cur = min(self.LV_205_14.v, LV_205_14_prev + self.valve_speed.v / 90.0 * self.model.settings.timestep)
                    else:
                        LV_205_14_cur = max(self.LV_205_14.v, LV_205_14_prev - self.valve_speed.v / 90.0 * self.model.settings.timestep)
                    self._CV_LV_205_14 = CV_data(1, self._D[7], LV_205_14_cur)
                    self._LV_205_14_prev = LV_205_14_cur

            # Set outputs
            # Note: most of the outputs are just recycled and aren't changed until the last timestep
            self.m_dot_in_out.v = self.m_dot_in.v
            self.T_in_out.v = self.T_in.v
            # alarm_high_high_pressure is preserved from the previous converged timestep value
            # (In Fortran the integer variable defaulted to 0; here we keep whatever was last set.)
            self.alarm_high_high_pressure.v = self.alarm_high_high_pressure.v
