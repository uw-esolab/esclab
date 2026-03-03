"""Deaerator-Pump component model (Type 6011)."""

import math

import eeslib.fluid_properties as fp

from esclab.simulate import Component
from esclab.components.flownetwork.esol6015_helpers import drhodhcp, drhodpch, dudhcp, dudpch, tank_level


class DeaeratorPump(Component):
    """
    Object: Deaerator-Pump
    Simulation Studio Model: ESOL6011-DA-PUMP

    Three-pump deaerator system (one variable-speed pump, two constant-speed pumps)
    with a horizontal cylindrical tank model. Tracks tank pressure, enthalpy, level,
    and mass using a 4th-order Runge-Kutta integration scheme, and determines LP
    turbine steam extraction for deaeration.

    Parameters (46 total)
    ---------------------
    PUMP CURVE:      Coef_A*Flow^2 + Coef_B*Flow + Coef_C = Pump Head
    EFFICIENCY CURVE: Eta_Coef_A*Flow^4 + Eta_Coef_B*Flow^3 + Eta_Coef_C*Flow^2 + Eta_Coef_D*Flow
    NPSH CURVE:      NPSH_A*Flow^3 + NPSH_B*Flow^2 + NPSH_C*Flow + NPSH_D

    Inputs (18 total)
    -----------------
    See Input declarations below.

    Outputs (43 total)
    ------------------
    See Output declarations below.
    """

    # -------------------------------------------------------------------------
    # Parameters
    # -------------------------------------------------------------------------

    # Variable Speed Pump 1 - Pump Curve Coefficients
    P1_Coef_A = Component.Parameter()   # Variable Speed Pump 1 - Pump Curve Coefficient A
    P1_Coef_B = Component.Parameter()   # Variable Speed Pump 1 - Pump Curve Coefficient B
    P1_Coef_C = Component.Parameter()   # Variable Speed Pump 1 - Pump Curve Coefficient C
    # Variable Speed Pump 1 - Efficiency Curve Coefficients
    P1_Eta_A = Component.Parameter()    # Variable Speed Pump 1 - Efficiency Curve Coefficient A
    P1_Eta_B = Component.Parameter()    # Variable Speed Pump 1 - Efficiency Curve Coefficient B
    P1_Eta_C = Component.Parameter()    # Variable Speed Pump 1 - Efficiency Curve Coefficient C
    P1_Eta_D = Component.Parameter()    # Variable Speed Pump 1 - Efficiency Curve Coefficient D
    # Variable Speed Pump 1 - NPSH Coefficients
    P1_NPSH_A = Component.Parameter()   # Variable Speed Pump 1 - Net Positive Suction Head Coefficient A
    P1_NPSH_B = Component.Parameter()   # Variable Speed Pump 1 - Net Positive Suction Head Coefficient B
    P1_NPSH_C = Component.Parameter()   # Variable Speed Pump 1 - Net Positive Suction Head Coefficient C
    P1_NPSH_D = Component.Parameter()   # Variable Speed Pump 1 - Net Positive Suction Head Coefficient D

    # Constant Speed Pump 2 - Pump Curve Coefficients
    P2_Coef_A = Component.Parameter()   # Constant Speed Pump 2 - Pump Curve Coefficient A
    P2_Coef_B = Component.Parameter()   # Constant Speed Pump 2 - Pump Curve Coefficient B
    P2_Coef_C = Component.Parameter()   # Constant Speed Pump 2 - Pump Curve Coefficient C
    # Constant Speed Pump 2 - Efficiency Curve Coefficients
    P2_Eta_A = Component.Parameter()    # Constant Speed Pump 2 - Efficiency Curve Coefficient A
    P2_Eta_B = Component.Parameter()    # Constant Speed Pump 2 - Efficiency Curve Coefficient B
    P2_Eta_C = Component.Parameter()    # Constant Speed Pump 2 - Efficiency Curve Coefficient C
    P2_Eta_D = Component.Parameter()    # Constant Speed Pump 2 - Efficiency Curve Coefficient D
    # Constant Speed Pump 2 - NPSH Coefficients
    P2_NPSH_A = Component.Parameter()   # Constant Speed Pump 2 - Net Positive Suction Head Coefficient A
    P2_NPSH_B = Component.Parameter()   # Constant Speed Pump 2 - Net Positive Suction Head Coefficient B
    P2_NPSH_C = Component.Parameter()   # Constant Speed Pump 2 - Net Positive Suction Head Coefficient C
    P2_NPSH_D = Component.Parameter()   # Constant Speed Pump 2 - Net Positive Suction Head Coefficient D

    # Constant Speed Pump 3 - Pump Curve Coefficients
    P3_Coef_A = Component.Parameter()   # Constant Speed Pump 3 - Pump Curve Coefficient A
    P3_Coef_B = Component.Parameter()   # Constant Speed Pump 3 - Pump Curve Coefficient B
    P3_Coef_C = Component.Parameter()   # Constant Speed Pump 3 - Pump Curve Coefficient C
    # Constant Speed Pump 3 - Efficiency Curve Coefficients
    P3_Eta_A = Component.Parameter()    # Constant Speed Pump 3 - Efficiency Curve Coefficient A
    P3_Eta_B = Component.Parameter()    # Variable Speed Pump 3 - Efficiency Curve Coefficient B
    P3_Eta_C = Component.Parameter()    # Variable Speed Pump 3 - Efficiency Curve Coefficient C
    P3_Eta_D = Component.Parameter()    # Variable Speed Pump 3 - Efficiency Curve Coefficient D
    # Constant Speed Pump 3 - NPSH Coefficients
    P3_NPSH_A = Component.Parameter()   # Constant Speed Pump 3 - Net Positive Suction Head Coefficient A
    P3_NPSH_B = Component.Parameter()   # Constant Speed Pump 3 - Net Positive Suction Head Coefficient B
    P3_NPSH_C = Component.Parameter()   # Constant Speed Pump 3 - Net Positive Suction Head Coefficient C
    P3_NPSH_D = Component.Parameter()   # Constant Speed Pump 3 - Net Positive Suction Head Coefficient D

    # Tank geometry and initial conditions
    D_tank = Component.Parameter()          # Diameter of the Deaerator Tank
    Length_tank = Component.Parameter()     # Length of the Deaerator Tank
    Length_tank2pump = Component.Parameter()  # Length of the downcomer from the tank to the pumps
    P_tank_ini = Component.Parameter()      # Initial Pressure of the tank at t = 0 [s]
    L_tank_ini = Component.Parameter()      # Initial Level of the tank at t = 0[s]

    # Flow and extraction limits
    m_dot_vent_frac = Component.Parameter()     # percent of incoming water that is vented out of deaerator tank
    m_dot_LPB1_max = Component.Parameter()  # Maximum flow rate coming from Low Pressure Turbine Stage [kg/s]
    DA_ss_LPB1 = Component.Parameter()      # max increase/decrease in extraction flow rate into the deaerator
    extraction_tol = Component.Parameter()  # Tolerance for Boiler Feedwater Heater extraction flow convergence

    # Level alarm/trip setpoints
    LL_Alarm_sp = Component.Parameter()    # Level in tank when the low level alarm switches on
    LL_Trip_sp = Component.Parameter()     # Level in the tank when the low level trip begins
    HL_Alarm_sp = Component.Parameter()    # Level in the tank when the high level alarm switches on
    HL_Trip_sp = Component.Parameter()     # Level in the tank when the high level trip begins

    # -------------------------------------------------------------------------
    # Inputs
    # -------------------------------------------------------------------------

    # Turbine Power Inputs (1 == ON, 0 == OFF)
    Turbine_ON = Component.Input()
    # Power Input for Variable Speed Pump (PUMP 1) (1==ON, 0 == OFF)
    Power_VarPump = Component.Input()
    # Power Input for Constant Speed Pump 2 (1==ON,0==OFF)
    Power_Pump2 = Component.Input()
    # Power Input for Constant Speed Pump 3 (1==ON, 0==OFF)
    Power_Pump3 = Component.Input()
    # Pump Speed for Variable Speed Pump
    pump_speed = Component.Input()
    # Mass flow rate of feedwater entering the deaerator
    m_dot_fw_in = Component.Input()
    # Pressure of feedwater entering the deaerator
    P_fw_in = Component.Input()
    # Enthalpy of feedwater entering the deaerator
    h_fw_in = Component.Input()
    # Pressure of the SD to determine operating point on pump curve
    P_SD = Component.Input()
    # Pressure of the piping system
    P_piping_sys = Component.Input()
    # Mass flow rate coming from the turbine bypass type
    m_dot_TB = Component.Input()
    # Pressure of the flow from the turbine bypass type
    P_TB = Component.Input()
    # enthalpy of the flow from the turbine bypass type
    h_TB = Component.Input()
    # Flow draining from the high pressure boiler feedwater heaters
    m_dot_HPFWH = Component.Input()
    # Pressure of flow draining from high pressure boiler feedwater heaters
    P_HPFWH = Component.Input()
    # Enthalpy of flow draining from high pressure boiler feedwater heaters
    h_HPFWH = Component.Input()
    # Pressure of steam extraction from low pressure turbine
    P_LPB1 = Component.Input()
    # enthalpy of steam extraction from low pressure turbine
    h_LPB1 = Component.Input()

    # -------------------------------------------------------------------------
    # Outputs
    # -------------------------------------------------------------------------

    # Combined pump discharge (outputs 1-5)
    m_dot_pump_out = Component.Output()     # output 1  - mass flow rate leaving all 3 pumps [kg/s]
    Vol_dot_pump = Component.Output()       # output 2  - Volumetric Flow Rate leaving the pumps [m^3/s]
    P_pump_out = Component.Output()         # output 3  - Pressure leaving the pumps [Pa]
    h_pump_out = Component.Output()         # output 4  - enthalpy leaving the pumps [J/kg]
    T_pump_out = Component.Output()         # output 5  - Temperature leaving the pumps [K]

    # Deaerator steam/bleed flows (outputs 6-8)
    m_dot_LPB1 = Component.Output()         # output 6  - extraction flow from low pressure turbine stage 1 [kg/s]
    m_dot_vent_out = Component.Output()     # output 7  - mass of steam released from deaerator [kg/s]
    h_vent = Component.Output()             # output 8  - enthalpy of steam released from deaerator [J/kg]

    # Tank state at START of timestep (updated at convergence, outputs 9-12)
    m_tank = Component.Output()             # output 9  - mass in the tank at the beginning of the timestep [kg]
    P_tank = Component.Output()             # output 10 - Pressure of the Deaerator during this timestep [Pa]
    L_tank = Component.Output()             # output 11 - Level of the Tank during this timestep [m]
    h_tank = Component.Output()             # output 12 - Enthalpy value of the tank [J/kg]

    # Tank state at END of timestep (iteratively computed, outputs 13-17)
    m_tank_new = Component.Output()         # output 13 - total mass in the tank for the next timestep [kg]
    T_tank_new = Component.Output()         # output 14 - Temperature of the tank for the next timestep [K]
    L_tank_new = Component.Output()         # output 15 - Level of the tank for the next timestep [m]
    P_tank_new = Component.Output()         # output 16 - Pressure of the tank for the next timestep [Pa]
    h_tank_new = Component.Output()         # output 17 - Enthalpy of the tank for the next timestep [J/kg]

    # Individual pump mass flows (outputs 18-20)
    m_dot_P1 = Component.Output()           # output 18 - Mass flow rate out of pump 1 [kg/s]
    m_dot_P2 = Component.Output()           # output 19 - Mass flow rate out of pump 2 [kg/s]
    m_dot_P3 = Component.Output()           # output 20 - Mass flow rate out of pump 3 [kg/s]

    # Pump power and efficiency (outputs 21-27)
    W_dot_total = Component.Output()        # output 21 - Total Pump Power needed [W]
    W_dot_P1 = Component.Output()           # output 22 - Power input to pump 1 [W]
    Eta_P1 = Component.Output()             # output 23 - Efficiency of Pump 1 [-]
    W_dot_P2 = Component.Output()           # output 24 - Power input to pump 2 [W]
    Eta_P2 = Component.Output()             # output 25 - Efficiency of Pump 2 [-]
    W_dot_P3 = Component.Output()           # output 26 - Power input to pump 3 [W]
    Eta_P3 = Component.Output()             # output 27 - Efficiency of Pump 3 [-]

    # Pump 1 secant-method solver state (outputs 28-30)
    P1_point_1x = Component.Output()        # output 28 - Previous Flow Rate [m^3/s]
    P1_point_1y = Component.Output()        # output 29 - Error associated with previous flow rate [m]
    P1_point_2x = Component.Output()        # output 30 - Flow Rate sent out during this iteration [m^3/s]

    # Pump 2 secant-method solver state (outputs 31-33)
    P2_point_1x = Component.Output()        # output 31 - Previous Flow Rate [m^3/s]
    P2_point_1y = Component.Output()        # output 32 - Error associated with previous flow rate [m]
    P2_point_2x = Component.Output()        # output 33 - Flow Rate sent out during this iteration [m^3/s]

    # Pump 3 secant-method solver state (outputs 34-36)
    P3_point_1x = Component.Output()        # output 34 - Previous Flow Rate [m^3/s]
    P3_point_1y = Component.Output()        # output 35 - Error associated with previous flow rate [m]
    P3_point_2x = Component.Output()        # output 36 - Flow Rate sent out during this iteration [m^3/s]

    # Pump cavitation trips (outputs 37-39)
    P1_trip = Component.Output()            # output 37 - Pump1 Trip (1 = cavitation, 0 = no problem) [-]
    P2_trip = Component.Output()            # output 38 - Pump2 Trip (1 = cavitation, 0 = no problem) [-]
    P3_trip = Component.Output()            # output 39 - Pump3 Trip (1 = cavitation, 0 = no problem) [-]

    # Level alarm/trip signals (outputs 40-43)
    LL_Alarm_out = Component.Output()       # output 40 - Low Level Alarm Signal (1 = Alarm, 0 = no problem) [-]
    LL_Trip_out = Component.Output()        # output 41 - Low Level Trip Signal (1 = Trip, 0 = no problem) [-]
    HL_Alarm_out = Component.Output()       # output 42 - High Level Alarm Signal (1 = Alarm, 0 = no problem) [-]
    HL_Trip_out = Component.Output()        # output 43 - High Level Trip Signal (1= Trip,0 = no problem) [-]

    # -------------------------------------------------------------------------

    def calculate(self):

        # -----------------------------------------------------------------------
        # Perform Any "After Convergence" Manipulations That May Be Required
        # at the End of Each Timestep
        if self.model.is_converged:
            P_tank_val = self.P_tank.v
            L_tank_val = self.L_tank.v

            P_ref = 87726.1   # vapor pressure of steam at 100 degrees C
            T_ref = 373.0
            D_pump_inlet = 0.1524

            # Check for pump cavitation through pump 1
            if self.Power_VarPump.v == 1.0:
                # CONVERTED-NEEDS UNITS CHECK: removed /1000 Pa->kPa; eeslib uses Pa
                T_v = fp.temperature("water", P=P_tank_val, Q=1.0)
                P_pump_in = P_tank_val + (L_tank_val + self.Length_tank2pump.v) * 9.81 * 1000.0
                Q_dot = self.m_dot_P1.v / 1000.0
                vel = Q_dot / (3.14 / 4.0 * D_pump_inlet ** 2.0)
                lnP1P2 = 8.314 * (1 / T_ref - 1 / T_v)   # Clausius Clapeyron Equation to find vapor pressure of water
                P_v = P_ref * math.exp(lnP1P2)
                # CONVERTED-NEEDS UNITS CHECK: removed /1000 Pa->kPa; eeslib uses Pa
                rho_v = fp.density("water", P=P_v, Q=1.0)
                # Speed-scaled NPSH curve coefficients
                A = self.P1_NPSH_A.v * self.pump_speed.v ** 3.0
                B = self.P1_NPSH_B.v * self.pump_speed.v ** 2.0
                C = self.P1_NPSH_C.v * self.pump_speed.v
                D = self.P1_NPSH_D.v
                NPSHr = A * Q_dot ** 3.0 + B * Q_dot ** 2.0 + C * Q_dot + D
                NPSHa = P_pump_in / 1000.0 / 9.81 + vel ** 2.0 / (2.0 * 9.81)  # - P_v/(rho_v * 9.81)
                if NPSHa > NPSHr:
                    cav = 0.0   # No Cavitation
                else:
                    cav = 1.0   # Cavitation
                self.P1_trip.v = cav
            else:
                self.P1_trip.v = 0.0

            # Check for Cavitation through pump 2
            if self.Power_Pump2.v == 1.0:
                # CONVERTED-NEEDS UNITS CHECK: removed /1000 Pa->kPa; eeslib uses Pa
                T_v = fp.temperature("water", P=P_tank_val, Q=1.0)
                P_pump_in = P_tank_val + (L_tank_val + self.Length_tank2pump.v) * 9.81 * 1000.0
                Q_dot = self.m_dot_P2.v / 1000.0
                vel = Q_dot / (3.14 / 4.0 * D_pump_inlet ** 2.0)
                lnP1P2 = 8.314 * (1 / T_ref - 1 / T_v)   # Clausius Clapeyron Equation to find vapor pressure of water
                P_v = P_ref * math.exp(lnP1P2)
                # CONVERTED-NEEDS UNITS CHECK: removed /1000 Pa->kPa; eeslib uses Pa
                rho_v = fp.density("water", P=P_v, Q=1.0)
                A = self.P2_NPSH_A.v * self.pump_speed.v ** 3.0
                B = self.P2_NPSH_B.v * self.pump_speed.v ** 2.0
                C = self.P2_NPSH_C.v * self.pump_speed.v
                D = self.P2_NPSH_D.v
                NPSHr = A * Q_dot ** 3.0 + B * Q_dot ** 2.0 + C * Q_dot + D
                NPSHa = P_pump_in / 1000.0 / 9.81 + vel ** 2.0 / (2.0 * 9.81)  # - P_v/(rho_v * 9.81)
                if NPSHa > NPSHr:
                    cav = 0.0   # No Cavitation
                else:
                    cav = 1.0   # Cavitation
                self.P2_trip.v = cav
            else:
                self.P2_trip.v = 0.0

            # check for cavitation through pump 3
            if self.Power_Pump3.v == 1.0:
                # CONVERTED-NEEDS UNITS CHECK: removed /1000 Pa->kPa; eeslib uses Pa
                T_v = fp.temperature("water", P=P_tank_val, Q=1.0)
                P_pump_in = P_tank_val + (L_tank_val + self.Length_tank2pump.v) * 9.81 * 1000.0
                Q_dot = self.m_dot_P3.v / 1000.0
                vel = Q_dot / (3.14 / 4.0 * D_pump_inlet ** 2.0)
                lnP1P2 = 8.314 * (1 / T_ref - 1 / T_v)   # Clausius Clapeyron Equation to find vapor pressure of water
                P_v = P_ref * math.exp(lnP1P2)
                # CONVERTED-NEEDS UNITS CHECK: removed /1000 Pa->kPa; eeslib uses Pa
                rho_v = fp.density("water", P=P_v, Q=1.0)
                # TODO-NEEDS CONVERSION REVIEW: Fortran uses P2_NPSH coefficients here instead of P3_NPSH; ported faithfully
                A = self.P2_NPSH_A.v * self.pump_speed.v ** 3.0
                B = self.P2_NPSH_B.v * self.pump_speed.v ** 2.0
                C = self.P2_NPSH_C.v * self.pump_speed.v
                D = self.P2_NPSH_D.v
                NPSHr = A * Q_dot ** 3.0 + B * Q_dot ** 2.0 + C * Q_dot + D
                NPSHa = P_pump_in / 1000.0 / 9.81 + vel ** 2.0 / (2.0 * 9.81)  # - P_v/(rho_v * 9.81)
                if NPSHa > NPSHr:
                    cav = 0.0   # No Cavitation
                else:
                    cav = 1.0   # Cavitation
                self.P3_trip.v = cav
            else:
                # TODO-NEEDS CONVERSION REVIEW: Fortran sets output 39 to 'cav' (last computed value) even in OFF case; ported faithfully
                self.P3_trip.v = cav

            # Check Level in tank for alarms or trips
            # check alarm and trip states
            # TODO-NEEDS CONVERSION REVIEW: Fortran compares alarm setpoint to L_tank_ini (initial level param) rather than current L_tank; ported faithfully
            if self.LL_Alarm_sp.v >= self.L_tank_ini.v:
                LL_Alarm = 1.0
                if self.LL_Trip_sp.v >= self.L_tank_ini.v:
                    LL_Trip = 1.0
                else:
                    LL_Trip = 0.0
            else:
                LL_Alarm = 0.0
                LL_Trip = 0.0
            if self.HL_Alarm_sp.v <= self.L_tank_ini.v:
                HL_Alarm = 1.0
                if self.HL_Trip_sp.v <= self.L_tank_ini.v:
                    HL_Trip = 1.0
                else:
                    HL_Trip = 0.0
            else:
                HL_Alarm = 0.0
                HL_Trip = 0.0
            self.LL_Alarm_out.v = LL_Alarm   # Low Level Alarm Signal (1 = Alarm, 0 = no problem)
            self.LL_Trip_out.v = LL_Trip     # Low Level Trip Signal (1 = Trip, 0 = no problem)
            self.HL_Alarm_out.v = HL_Alarm   # High Level Alarm Signal (1 = Alarm, 0 = no problem)
            # TODO-NEEDS CONVERSION REVIEW: Fortran sets output 43 (HL_Trip) to HL_Alarm value; ported faithfully
            self.HL_Trip_out.v = HL_Alarm    # High Level Trip Signal (1= Trip,0 = no problem)

            # Update Pressure, Enthalpy, Level and mass in tank for next timestep
            # mass in tank at the end of last timestep is the mass in tank the beginning of this timestep
            self.m_tank.v = self.m_tank_new.v
            # Pressure at the end of the last timestep is the pressure at the beginning of this timestep
            self.P_tank.v = self.P_tank_new.v
            # Level ...
            self.L_tank.v = self.L_tank_new.v
            # Enthalpy ...
            self.h_tank.v = self.h_tank_new.v
            return

        # -----------------------------------------------------------------------
        # Do All of the First Timestep Manipulations Here
        # There Are No Iterations at the Initial Time
        if self.model.is_first_step:

            # PUMP CURVE --> Coef_A*Flow^2 +Coef_B*Flow + Coef_C = Pump Head
            # EFFICIENCY CURVE --> Eta_Coef_A*Flow^4 +Eta_Coef_B*Flow^3+Eta_Coef_C*Flow^2+Eta_Coef_D*Flow
            # NPSH CURVE --> NPSH_A * Flow^3 + NPSH_B*Flow^2 + NPSH_C*Flow + NPSH_D

            if self.Turbine_ON.v == 1.0:
                m_dot_LPB1 = self.m_dot_LPB1_max.v
            else:
                m_dot_LPB1 = 0.0

            # Finding initial tank enthalpy
            # CONVERTED-NEEDS UNITS CHECK: removed /1000 Pa->kPa; eeslib uses Pa
            T_tank = fp.temperature("water", P=self.P_tank_ini.v, Q=1.0)
            rho_tank_g = fp.density("water", P=self.P_tank_ini.v, Q=1.0)
            rho_tank_f = fp.density("water", P=self.P_tank_ini.v, Q=0.0)
            R_tank = self.D_tank.v / 2.0
            Vol_tank = 3.14 * R_tank ** 2.0 * self.Length_tank.v
            Area_liquid = (
                math.acos((R_tank - self.L_tank_ini.v) / R_tank) * R_tank ** 2.0
                - (R_tank - self.L_tank_ini.v)
                * math.sqrt(2.0 * R_tank * self.L_tank_ini.v - self.L_tank_ini.v ** 2.0)
            )
            m_tank_f = Area_liquid * self.Length_tank.v * rho_tank_f
            m_tank_g = (Vol_tank - Area_liquid * self.Length_tank.v) * rho_tank_g
            m_tank_tot = m_tank_f + m_tank_g
            x_tank = m_tank_g / m_tank_tot
            # CONVERTED-NEEDS UNITS CHECK: removed /1000 Pa->kPa; eeslib uses Pa; removed *1000 kJ/kg->J/kg; eeslib returns J/kg
            h_tank_ini = fp.enthalpy("water", P=self.P_tank_ini.v, Q=x_tank)   # J/kg
            T_tank_new = fp.temperature("water", P=self.P_tank_ini.v, Q=x_tank)

            P_bottom = self.P_tank_ini.v + self.L_tank_ini.v * rho_tank_f * 9.81
            m_dot_pump = 0.0001
            Q_dot_pump = m_dot_pump / rho_tank_f
            # CONVERTED-NEEDS UNITS CHECK: removed /1000 Pa->kPa; eeslib uses Pa; removed *1000 kJ/kg->J/kg; eeslib returns J/kg
            h_pump_in = fp.enthalpy("water", P=P_bottom, T=T_tank)   # J/kg
            s_pump_in = fp.entropy("water", P=P_bottom, T=T_tank)    # J/(kg·K)
            # Finding enthalpy leaving the pump assuming minimal flow leaving
            P_pump_out = self.P1_Coef_C.v * rho_tank_f * 9.81 + self.P_tank_ini.v + self.L_tank_ini.v * rho_tank_f * 9.81
            s_pump_out_s = s_pump_in   # Converting from KJ/kg-K to J/kg-K
            # CONVERTED-NEEDS UNITS CHECK: removed /1000 Pa->kPa; eeslib uses Pa; removed *1000 kJ/kg->J/kg; eeslib returns J/kg
            h_pump_out_s = fp.enthalpy("water", P=P_pump_out, S=s_pump_out_s)   # J/kg
            Eta_P1 = 0.01    # pump is running at minimal efficiency
            W_dot_pump_s = m_dot_pump * (h_pump_out_s - h_pump_in)
            W_dot_P1 = W_dot_pump_s / Eta_P1
            h_pump_out = (m_dot_pump * h_pump_in + W_dot_P1) / m_dot_pump
            # CONVERTED-NEEDS UNITS CHECK: removed /1000 Pa->kPa; eeslib uses Pa; removed /1000 J/kg->kJ/kg; eeslib uses J/kg
            T_pump_out = fp.temperature("water", P=P_pump_out, H=h_pump_out)

            # check alarm and trip states
            if self.LL_Alarm_sp.v >= self.L_tank_ini.v:
                LL_Alarm = 1.0
                if self.LL_Trip_sp.v >= self.L_tank_ini.v:
                    LL_Trip = 1.0
                else:
                    LL_Trip = 0.0
            else:
                LL_Alarm = 0.0
                LL_Trip = 0.0
            if self.HL_Alarm_sp.v <= self.L_tank_ini.v:
                HL_Alarm = 1.0
                if self.HL_Trip_sp.v <= self.L_tank_ini.v:
                    HL_Trip = 1.0
                else:
                    HL_Trip = 0.0
            else:
                HL_Alarm = 0.0
                HL_Trip = 0.0

            # Set the Initial Values of the Outputs
            self.m_dot_pump_out.v = m_dot_pump      # mass leaving the pump
            self.Vol_dot_pump.v = Q_dot_pump         # Volumetric flow rate leaving the pump
            self.P_pump_out.v = P_pump_out           # Pressure leaving the pump
            self.h_pump_out.v = h_pump_out           # enthalpy leaving the pump
            self.T_pump_out.v = T_pump_out           # temperature leaving the pump
            self.m_dot_LPB1.v = m_dot_LPB1          # extraction flow from low pressure turbine stage 1
            self.m_dot_vent_out.v = 0.0              # rate of steam released from deaerator
            self.h_vent.v = 0.0                      # enthalpy of steam released from deaerator
            self.m_tank.v = m_tank_tot               # mass in the tank at the start of the simulation
            self.P_tank.v = self.P_tank_ini.v        # Pressure of tank at beginning of the simulation
            self.L_tank.v = self.L_tank_ini.v        # Level of tank at the beginning of the simulation
            self.h_tank.v = h_tank_ini               # Enthalpy value of tank at beginning of simulation
            self.m_tank_new.v = m_tank_tot           # mass in the tank that will change during iterating
            self.T_tank_new.v = T_tank_new           # Temperature of tank that will change during iterating
            self.L_tank_new.v = float('nan')         # New tank Level computed during iteration process (not yet computed)
            self.P_tank_new.v = float('nan')         # New tank pressure computed during iteration process (not yet computed)
            self.h_tank_new.v = float('nan')         # New tank enthalpy computed during iteration process (not yet computed)
            self.P1_trip.v = 0.0                     # Pump1 Trip initially set to 0
            self.P2_trip.v = 0.0                     # Pump2 Trip initially set to 0
            self.P3_trip.v = 0.0                     # Pump3 Trip initially set to 0
            self.LL_Alarm_out.v = LL_Alarm           # Low Level Alarm Signal (1 = Alarm, 0 = no problem)
            self.LL_Trip_out.v = LL_Trip             # Low Level Trip Signal (1 = Trip, 0 = no problem)
            self.HL_Alarm_out.v = HL_Alarm           # High Level Alarm Signal (1 = Alarm, 0 = no problem)
            # TODO-NEEDS CONVERSION REVIEW: Fortran sets output 43 (HL_Trip) to HL_Alarm value; ported faithfully
            self.HL_Trip_out.v = HL_Alarm            # High Level Trip Signal (1= Trip,0 = no problem)
            return

        # -----------------------------------------------------------------------
        # Read the Inputs and Parameters
        # (Parameters are class members; inputs are accessed via self.input.v)

        # Converting timestep from hr to s
        ts = self.model.timestep * 3600
        Vol_tank = 3.14 / 4.0 * self.D_tank.v ** 2.0 * (self.Length_tank.v - self.D_tank.v) + 4.0 / 3.0 * 3.14 * (self.D_tank.v / 2.0) ** 3.0
        tol = 1.0    # Curves are within 1 kPa of each other
        LR = 0.2     # Learning Rate

        # -----------------------------------------------------------------------
        # Current Tank Values (Level, Pressure, Temperature, density of liquid water)
        # Level of the Tank during this timestep
        L_tank = self.L_tank.v
        # Pressure of the Deaerator during this timestep
        P_tank = self.P_tank.v
        # CONVERTED-NEEDS UNITS CHECK: removed /1000 Pa->kPa; eeslib uses Pa
        rho_tank_f = fp.density("water", P=P_tank, Q=0.0)
        T_tank = fp.temperature("water", P=P_tank, Q=0.0)
        # Previous Pressure leaving the pumps
        P_pump_prev = self.P_pump_out.v

        # -----------------------------------------------------------------------
        # Step 1: Solve for mass flow out of Variable Speed Pump 1

        if self.Power_VarPump.v == 1.0:   # Variable Speed Pump is ON

            if self.pump_speed.v > 0.0:   # Pump Speed is Possible
                pump_speed = min(self.pump_speed.v, 1.0)
                pump_speed = max(pump_speed, 0.01)   # Lowest Pump Speed allowed for Type
                # Adjust Pump Curve Coefficients based on Pump Speed
                A = self.P1_Coef_A.v
                B = self.P1_Coef_B.v * pump_speed
                C = self.P1_Coef_C.v * pump_speed ** 2.0
                # Scale efficiency coefficients by pump speed
                P1_Eta_A = self.P1_Eta_A.v / pump_speed ** 4.0
                P1_Eta_B = self.P1_Eta_B.v / pump_speed ** 3.0
                P1_Eta_C = self.P1_Eta_C.v / pump_speed ** 2.0
                P1_Eta_D = self.P1_Eta_D.v / pump_speed

                if self.model.iteration == 0:   # First Iteration in the timestep
                    # Previous flow out of pump 1
                    m_dot_P1 = self.m_dot_P1.v

                    # Find error for previous flow rate
                    P1_point_1x = m_dot_P1 / rho_tank_f   # Previous flow rate
                    # Previous pump head used in last iteration
                    P1_pump_head = P_pump_prev / rho_tank_f / 9.81 - P_tank / rho_tank_f / 9.81 - L_tank - self.Length_tank2pump.v
                    System_headloss = (P_pump_prev - self.P_piping_sys.v) / rho_tank_f / 9.81
                    # Previous error
                    P1_point_1y = P1_pump_head - self.P_SD.v / rho_tank_f / 9.81 - System_headloss + P_tank / rho_tank_f / 9.81 + L_tank + self.Length_tank2pump.v

                    # Maximum flow allowed from pump
                    Q_max_P1 = max(
                        (-B + math.sqrt(B ** 2.0 - 4.0 * A * C)) / (2.0 * A),
                        (-B - math.sqrt(B ** 2.0 - 4.0 * A * C)) / (2.0 * A)
                    )

                    # Solve for new flow leaving pump 1
                    if abs(P1_point_1y) <= tol:
                        Q_new_P1 = P1_point_1x
                    elif P1_point_1y > 0.0:   # increase flow rate
                        Q_new_P1 = min(P1_point_1x + 0.001, Q_max_P1)   # Increase flow but do not go over maximum flow
                    else:   # decrease flow rate
                        Q_new_P1 = max(P1_point_1x - 0.001, 0.000001)   # Move to minimum flow pump can provide

                    P1_point_2x = Q_new_P1
                    P_pump_in = P_tank + rho_tank_f * 9.81 * (L_tank + self.Length_tank2pump.v)
                    # CONVERTED-NEEDS UNITS CHECK: removed /1000 Pa->kPa; eeslib uses Pa; removed *1000 kJ/kg->J/kg; eeslib returns J/kg
                    h_pump_in = fp.enthalpy("water", P=P_pump_in, T=T_tank)   # J/kg
                    s_pump_in = fp.entropy("water", P=P_pump_in, T=T_tank)    # J/(kg·K)
                    # Mass flow rate leaving pump 1
                    m_dot_P1 = max(Q_new_P1 * rho_tank_f, 0.000001)
                    P_P1_out = (A * Q_new_P1 ** 2.0 + B * Q_new_P1 + C) * rho_tank_f * 9.81 + P_pump_in

                    # Solve for pump efficiency, power and enthalpy leaving
                    # if pump speed is 0 division by zero will occur
                    pump_speed = max(pump_speed, 0.0001)
                    Eta_P1 = max(
                        P1_Eta_A * Q_new_P1 ** 4.0
                        + P1_Eta_B * Q_new_P1 ** 3.0
                        + P1_Eta_C * Q_new_P1 ** 2.0
                        + P1_Eta_D * Q_new_P1,
                        0.2
                    )
                    s_pump_out_s = s_pump_in
                    # CONVERTED-NEEDS UNITS CHECK: removed /1000 Pa->kPa; eeslib uses Pa; removed *1000 kJ/kg->J/kg; eeslib returns J/kg
                    h_pump_out_s = fp.enthalpy("water", P=P_P1_out, S=s_pump_out_s)   # J/kg
                    # Isentropic energy balance to find ideal work required
                    W_dot_pump_s = m_dot_P1 * (h_pump_out_s - h_pump_in)
                    W_dot_P1 = W_dot_pump_s / Eta_P1
                    # Actual Energy balance to find actual pump enthalpy out
                    h_P1_out = (m_dot_P1 * h_pump_in + W_dot_P1) / m_dot_P1

                    # Call Outputs
                    self.m_dot_P1.v = m_dot_P1       # Mass flow rate out of pump 1
                    self.W_dot_P1.v = W_dot_P1        # Power input to pump 1
                    self.Eta_P1.v = Eta_P1             # Efficiency of Pump 1
                    self.P1_point_1x.v = P1_point_1x  # Previous Flow Rate
                    self.P1_point_1y.v = P1_point_1y  # Error associated with previous flow rate
                    self.P1_point_2x.v = P1_point_2x  # Flow Rate sent out during this iteration

                else:
                    P1_point_1x = self.P1_point_1x.v
                    P1_point_1y = self.P1_point_1y.v
                    P1_point_2x = self.P1_point_2x.v

                    # Compute head loss in system
                    # Previous pump head used in last iteration
                    P1_pump_head = P_pump_prev / rho_tank_f / 9.81 - P_tank / rho_tank_f / 9.81 - L_tank - self.Length_tank2pump.v
                    System_headloss = (P_pump_prev - self.P_piping_sys.v) / rho_tank_f / 9.81
                    # Previous error
                    P1_point_2y = P1_pump_head - self.P_SD.v / rho_tank_f / 9.81 - System_headloss + P_tank / rho_tank_f / 9.81 + L_tank + self.Length_tank2pump.v

                    # Find maximum flow rate allowed for input pump speed
                    Q_max_P1 = max(
                        (-B + math.sqrt(B ** 2.0 - 4.0 * A * C)) / (2.0 * A),
                        (-B - math.sqrt(B ** 2.0 - 4.0 * A * C)) / (2.0 * A)
                    )

                    # Solve for new flow rate leaving pump 1
                    if P1_point_2x == P1_point_1x:
                        m = 0.0
                    else:
                        m = (P1_point_2y - P1_point_1y) / (P1_point_2x - P1_point_1x)
                    y_int = P1_point_1y - m * P1_point_1x

                    if abs(P1_point_2y) <= tol:
                        Q_new_P1 = P1_point_2x   # Keep flow rate the same
                    elif m != 0.0:
                        if P1_point_2y >= 0.0:   # If Error is Positive, pump head guess is too high so flow must increase
                            Q_new_P1 = max(-y_int / m, P1_point_2x)   # Ensure new Q_dot does not decrease
                            Q_new_P1 = min(Q_new_P1, Q_max_P1)         # Ensure Q_dot_new is not more than pump can provide based on curve
                            delta_Q = abs(Q_new_P1 - P1_point_2x)
                            Q_new_P1 = P1_point_2x + delta_Q * LR      # Increase Q_dot_new based on the learning rate
                        else:
                            Q_new_P1 = min(-y_int / m, P1_point_2x)   # Ensure new Q_dot does not increase
                            Q_new_P1 = max(Q_new_P1, 0.000001)          # Ensure new Q_dot is not negative
                            delta_Q = abs(Q_new_P1 - P1_point_2x)
                            Q_new_P1 = P1_point_2x - delta_Q * LR      # Increase Q_dot_new based on the learning rate
                    else:   # If Error is Negative, pump head guess is too low so flow must decrease
                        if P1_point_2y >= 0.0:   # Error is Positive and slope is equal to zero, pump head is too high so flow must increase
                            Q_new_P1 = min(P1_point_2x + 0.0001, Q_max_P1)
                        else:   # Error is Negative and slope is equal to zero, pump head is too low so flow must decrease
                            Q_new_P1 = max(P1_point_2x - 0.0001, 0.000001)

                    # Pressure into the pump is the pressure at the bottom of the condenser tank
                    P_pump_in = P_tank + rho_tank_f * 9.81 * (L_tank + self.Length_tank2pump.v)
                    # CONVERTED-NEEDS UNITS CHECK: removed /1000 Pa->kPa; eeslib uses Pa; removed *1000 kJ/kg->J/kg; eeslib returns J/kg
                    h_pump_in = fp.enthalpy("water", P=P_pump_in, T=T_tank)   # J/kg
                    s_pump_in = fp.entropy("water", P=P_pump_in, T=T_tank)    # J/(kg·K)
                    m_dot_P1 = max(Q_new_P1 * rho_tank_f, 0.0001)
                    P_P1_out = (A * Q_new_P1 ** 2.0 + B * Q_new_P1 + C) * rho_tank_f * 9.81 + P_pump_in

                    # Solve for pump efficiency, power and enthalpy leaving
                    # if pump speed is 0 division by zero will occur
                    pump_speed = max(pump_speed, 0.0001)
                    Eta_P1 = max(
                        (P1_Eta_A * Q_new_P1 ** 4.0
                         + P1_Eta_B * Q_new_P1 ** 3.0
                         + P1_Eta_C * Q_new_P1 ** 2.0
                         + P1_Eta_D * Q_new_P1) / 100.0,
                        0.2
                    )
                    s_pump_out_s = s_pump_in
                    # CONVERTED-NEEDS UNITS CHECK: removed /1000 Pa->kPa; eeslib uses Pa; removed *1000 kJ/kg->J/kg; eeslib returns J/kg
                    h_pump_out_s = fp.enthalpy("water", P=P_P1_out, S=s_pump_out_s)   # J/kg
                    # Isentropic energy balance to find ideal work required
                    W_dot_pump_s = m_dot_P1 * (h_pump_out_s - h_pump_in)
                    W_dot_P1 = W_dot_pump_s / Eta_P1
                    # Actual Energy balance to find actual pump enthalpy out
                    h_P1_out = (m_dot_P1 * h_pump_in + W_dot_P1) / m_dot_P1

                    P1_point_1x = P1_point_2x   # Saving 2nd point as 1st point for next iteration
                    P1_point_1y = P1_point_2y   # Saving 2nd point as 1st point for next iteration
                    P1_point_2x = Q_new_P1       # Saving new flow as 2nd point

                    # Call Outputs
                    self.m_dot_P1.v = m_dot_P1       # Mass flow rate out of pump 1
                    self.W_dot_P1.v = W_dot_P1        # Power input to pump 1
                    self.Eta_P1.v = Eta_P1             # Efficiency of Pump 1
                    self.P1_point_1x.v = P1_point_1x  # Previous Flow Rate
                    self.P1_point_1y.v = P1_point_1y  # Error associated with previous flow rate
                    self.P1_point_2x.v = P1_point_2x  # Flow Rate sent out during this iteration

            else:   # Pump Speed is 0 or not possible
                m_dot_P1 = 0.0
                W_dot_P1 = 0.0
                Eta_P1 = 0.0
                h_P1_out = 0.0
                self.m_dot_P1.v = m_dot_P1   # Mass flow rate out of pump 1
                self.W_dot_P1.v = W_dot_P1   # Power input to pump 1
                self.Eta_P1.v = Eta_P1        # Efficiency of Pump 1

        else:   # Variable Speed Pump is OFF
            m_dot_P1 = 0.0
            W_dot_P1 = 0.0
            Eta_P1 = 0.0
            h_P1_out = 0.0
            self.m_dot_P1.v = m_dot_P1   # Mass flow rate out of pump 1
            self.W_dot_P1.v = W_dot_P1   # Power input to pump 1
            self.Eta_P1.v = Eta_P1        # Efficiency of Pump 1

        # -----------------------------------------------------------------------
        # Step 2: Solve for mass flow out of Constant Speed Pump 2

        if self.Power_Pump2.v == 1:
            # Adjust Pump Curve Coefficients based on Pump Speed
            A = self.P2_Coef_A.v
            B = self.P2_Coef_B.v
            C = self.P2_Coef_C.v

            if self.model.iteration == 0:   # First Iteration in the timestep
                # Previous flow out of pump 2
                m_dot_P2 = self.m_dot_P2.v

                # Find error for previous flow rate
                P2_point_1x = m_dot_P2 / rho_tank_f   # Previous flow rate
                # Previous pump head used in last iteration
                P2_pump_head = P_pump_prev / rho_tank_f / 9.81 - P_tank / rho_tank_f / 9.81 - L_tank - self.Length_tank2pump.v
                System_headloss = (P_pump_prev - self.P_piping_sys.v) / rho_tank_f / 9.81
                # Previous error
                P2_point_1y = P2_pump_head - self.P_SD.v / rho_tank_f / 9.81 - System_headloss + P_tank / rho_tank_f / 9.81 + L_tank + self.Length_tank2pump.v

                # Maximum flow allowed from pump
                Q_max_P2 = max(
                    (-B + math.sqrt(B ** 2.0 - 4.0 * A * C)) / (2.0 * A),
                    (-B - math.sqrt(B ** 2.0 - 4.0 * A * C)) / (2.0 * A)
                )

                # Solve for new flow leaving pump 2
                if abs(P2_point_1y) <= tol:
                    Q_new_P2 = P2_point_1x
                elif P2_point_1y > 0.0:   # increase flow rate
                    Q_new_P2 = min(P2_point_1x + 0.001, Q_max_P2)   # Increase flow but do not go over maximum flow
                else:   # decrease flow rate
                    Q_new_P2 = max(P2_point_1x - 0.001, 0.000001)   # Move to minimum flow pump can provide

                P2_point_2x = Q_new_P2
                P_pump_in = P_tank + rho_tank_f * 9.81 * (L_tank + self.Length_tank2pump.v)
                # CONVERTED-NEEDS UNITS CHECK: removed /1000 Pa->kPa; eeslib uses Pa; removed *1000 kJ/kg->J/kg; eeslib returns J/kg
                h_pump_in = fp.enthalpy("water", P=P_pump_in, T=T_tank)   # J/kg
                s_pump_in = fp.entropy("water", P=P_pump_in, T=T_tank)    # J/(kg·K)
                # Mass flow rate leaving pump 2
                m_dot_P2 = max(Q_new_P2 * rho_tank_f, 0.000001)
                P_P2_out = (A * Q_new_P2 ** 2.0 + B * Q_new_P2 + C) * rho_tank_f * 9.81 + P_pump_in

                # Solve for pump efficiency, power and enthalpy leaving
                Eta_P2 = max(
                    self.P2_Eta_A.v * Q_new_P2 ** 4.0
                    + self.P2_Eta_B.v * Q_new_P2 ** 3.0
                    + self.P2_Eta_C.v * Q_new_P2 ** 2.0
                    + self.P2_Eta_D.v * Q_new_P2,
                    0.01
                )
                s_pump_out_s = s_pump_in
                # CONVERTED-NEEDS UNITS CHECK: removed /1000 Pa->kPa; eeslib uses Pa; removed *1000 kJ/kg->J/kg; eeslib returns J/kg
                h_pump_out_s = fp.enthalpy("water", P=P_P2_out, S=s_pump_out_s)   # J/kg
                # Isentropic energy balance to find ideal work required
                W_dot_pump_s = m_dot_P2 * (h_pump_out_s - h_pump_in)
                W_dot_P2 = W_dot_pump_s / Eta_P2
                # Actual Energy balance to find actual pump enthalpy out
                h_P2_out = (m_dot_P2 * h_pump_in + W_dot_P2) / m_dot_P2

                # Call Outputs
                self.m_dot_P2.v = m_dot_P2       # Mass flow rate out of pump 2
                self.W_dot_P2.v = W_dot_P2        # Power input to pump 2
                self.Eta_P2.v = Eta_P2             # Efficiency of Pump 2
                self.P2_point_1x.v = P2_point_1x  # Previous Flow Rate
                self.P2_point_1y.v = P2_point_1y  # Error associated with previous flow rate
                self.P2_point_2x.v = P2_point_2x  # Flow Rate sent out during this iteration

            else:
                P2_point_1x = self.P2_point_1x.v
                P2_point_1y = self.P2_point_1y.v
                P2_point_2x = self.P2_point_2x.v

                # Compute head loss in system
                # Previous pump head used in last iteration
                P2_pump_head = P_pump_prev / rho_tank_f / 9.81 - P_tank / rho_tank_f / 9.81 - L_tank - self.Length_tank2pump.v
                System_headloss = (P_pump_prev - self.P_piping_sys.v) / rho_tank_f / 9.81
                # Previous error
                P2_point_2y = P2_pump_head - self.P_SD.v / rho_tank_f / 9.81 - System_headloss + P_tank / rho_tank_f / 9.81 + L_tank + self.Length_tank2pump.v

                # Find maximum flow rate allowed for input pump speed
                Q_max_P2 = max(
                    (-B + math.sqrt(B ** 2.0 - 4.0 * A * C)) / (2.0 * A),
                    (-B - math.sqrt(B ** 2.0 - 4.0 * A * C)) / (2.0 * A)
                )

                # Solve for new flow rate leaving pump 1
                if P2_point_2x == P2_point_1x:
                    m = 0.0
                else:
                    m = (P2_point_2y - P2_point_1y) / (P2_point_2x - P2_point_1x)
                y_int = P2_point_1y - m * P2_point_1x

                if abs(P2_point_2y) <= tol:
                    Q_new_P2 = P2_point_2x   # Keep flow rate the same
                elif m != 0.0:
                    if P2_point_2y >= 0.0:   # If Error is Positive, pump head guess is too high so flow must increase
                        Q_new_P2 = max(-y_int / m, P2_point_2x)   # Ensure new Q_dot does not decrease
                        Q_new_P2 = min(Q_new_P2, Q_max_P2)         # Ensure Q_dot_new is not more than pump can provide based on curve
                        delta_Q = abs(Q_new_P2 - P2_point_2x)
                        Q_new_P2 = P2_point_2x + delta_Q * LR      # Increase Q_dot_new based on the learning rate
                    else:
                        Q_new_P2 = min(-y_int / m, P2_point_2x)   # Ensure new Q_dot does not increase
                        Q_new_P2 = max(Q_new_P2, 0.000001)          # Ensure new Q_dot is not negative
                        delta_Q = abs(Q_new_P2 - P2_point_2x)
                        Q_new_P2 = P2_point_2x - delta_Q * LR      # Increase Q_dot_new based on the learning rate
                else:   # If Error is Negative, pump head guess is too low so flow must decrease
                    if P2_point_2y >= 0.0:   # Error is Positive and slope is equal to zero, pump head is too high so flow must increase
                        Q_new_P2 = min(P2_point_2x + 0.0001, Q_max_P2)
                    else:   # Error is Negative and slope is equal to zero, pump head is too low so flow must decrease
                        Q_new_P2 = max(P2_point_2x - 0.0001, 0.000001)

                # Pressure into the pump is the pressure at the bottom of the condenser tank
                P_pump_in = P_tank + rho_tank_f * 9.81 * (L_tank + self.Length_tank2pump.v)
                # CONVERTED-NEEDS UNITS CHECK: removed /1000 Pa->kPa; eeslib uses Pa; removed *1000 kJ/kg->J/kg; eeslib returns J/kg
                h_pump_in = fp.enthalpy("water", P=P_pump_in, T=T_tank)   # J/kg
                s_pump_in = fp.entropy("water", P=P_pump_in, T=T_tank)    # J/(kg·K)
                m_dot_P2 = max(Q_new_P2 * rho_tank_f, 0.0001)
                P_P2_out = (A * Q_new_P2 ** 2.0 + B * Q_new_P2 + C) * rho_tank_f * 9.81 + P_pump_in

                # Solve for pump efficiency, power and enthalpy leaving
                Eta_P2 = max(
                    (self.P2_Eta_A.v * Q_new_P2 ** 4.0
                     + self.P2_Eta_B.v * Q_new_P2 ** 3.0
                     + self.P2_Eta_C.v * Q_new_P2 ** 2.0
                     + self.P2_Eta_D.v * Q_new_P2) / 100.0,
                    0.01
                )
                s_pump_out_s = s_pump_in
                # CONVERTED-NEEDS UNITS CHECK: removed /1000 Pa->kPa; eeslib uses Pa; removed *1000 kJ/kg->J/kg; eeslib returns J/kg
                h_pump_out_s = fp.enthalpy("water", P=P_P2_out, S=s_pump_out_s)   # J/kg
                # Isentropic energy balance to find ideal work required
                W_dot_pump_s = m_dot_P2 * (h_pump_out_s - h_pump_in)
                W_dot_P2 = W_dot_pump_s / Eta_P2
                # Actual Energy balance to find actual pump enthalpy out
                h_P2_out = (m_dot_P2 * h_pump_in + W_dot_P2) / m_dot_P2

                P2_point_1x = P2_point_2x   # Saving 2nd point as 1st point for next iteration
                P2_point_1y = P2_point_2y   # Saving 2nd point as 1st point for next iteration
                P2_point_2x = Q_new_P2       # Saving new flow as 2nd point

                # Call Outputs
                self.m_dot_P2.v = m_dot_P2       # Mass flow rate out of pump 1
                self.W_dot_P2.v = W_dot_P2        # Power input to pump 1
                self.Eta_P2.v = Eta_P2             # Efficiency of Pump 1
                self.P2_point_1x.v = P2_point_1x  # Previous Flow Rate
                self.P2_point_1y.v = P2_point_1y  # Error associated with previous flow rate
                self.P2_point_2x.v = P2_point_2x  # Flow Rate sent out during this iteration

        else:
            m_dot_P2 = 0.0
            W_dot_P2 = 0.0
            Eta_P2 = 0.0
            h_P2_out = 0.0
            self.m_dot_P2.v = m_dot_P2   # Mass flow rate out of pump 1
            self.W_dot_P2.v = W_dot_P2   # Power input to pump 1
            self.Eta_P2.v = Eta_P2        # Efficiency of Pump 1

        # -----------------------------------------------------------------------
        # Step 3: Solve for mass flow out of Constant Speed Pump 3

        if self.Power_Pump3.v == 1:
            # Adjust Pump Curve Coefficients based on Pump Speed
            A = self.P3_Coef_A.v
            B = self.P3_Coef_B.v
            C = self.P3_Coef_C.v

            if self.model.iteration == 0:   # First Iteration in the timestep
                # Previous flow out of pump 2
                m_dot_P3 = self.m_dot_P3.v

                # Find error for previous flow rate
                P3_point_1x = m_dot_P3 / rho_tank_f   # Previous flow rate
                # Previous pump head used in last iteration
                P3_pump_head = P_pump_prev / rho_tank_f / 9.81 - P_tank / rho_tank_f / 9.81 - L_tank - self.Length_tank2pump.v
                System_headloss = (P_pump_prev - self.P_piping_sys.v) / rho_tank_f / 9.81
                # Previous error
                P3_point_1y = P3_pump_head - self.P_SD.v / rho_tank_f / 9.81 - System_headloss + P_tank / rho_tank_f / 9.81 + L_tank + self.Length_tank2pump.v

                # Maximum flow allowed from pump
                Q_max_P3 = max(
                    (-B + math.sqrt(B ** 2.0 - 4.0 * A * C)) / (2.0 * A),
                    (-B - math.sqrt(B ** 2.0 - 4.0 * A * C)) / (2.0 * A)
                )

                # Solve for new flow leaving pump 2
                if abs(P3_point_1y) <= tol:
                    Q_new_P3 = P3_point_1x
                elif P3_point_1y > 0.0:   # increase flow rate
                    Q_new_P3 = min(P3_point_1x + 0.001, Q_max_P3)   # Increase flow but do not go over maximum flow
                else:   # decrease flow rate
                    Q_new_P3 = max(P3_point_1x - 0.001, 0.000001)   # Move to minimum flow pump can provide

                P3_point_2x = Q_new_P3
                P_pump_in = P_tank + rho_tank_f * 9.81 * (L_tank + self.Length_tank2pump.v)
                # CONVERTED-NEEDS UNITS CHECK: removed /1000 Pa->kPa; eeslib uses Pa; removed *1000 kJ/kg->J/kg; eeslib returns J/kg
                h_pump_in = fp.enthalpy("water", P=P_pump_in, T=T_tank)   # J/kg
                s_pump_in = fp.entropy("water", P=P_pump_in, T=T_tank)    # J/(kg·K)
                # Mass flow rate leaving pump 2
                m_dot_P3 = max(Q_new_P3 * rho_tank_f, 0.000001)
                P_P3_out = (A * Q_new_P3 ** 2.0 + B * Q_new_P3 + C) * rho_tank_f * 9.81 + P_pump_in

                # Solve for pump efficiency, power and enthalpy leaving
                Eta_P3 = max(
                    (self.P3_Eta_A.v * Q_new_P3 ** 4.0
                     + self.P3_Eta_B.v * Q_new_P3 ** 3.0
                     + self.P3_Eta_C.v * Q_new_P3 ** 2.0
                     + self.P3_Eta_D.v * Q_new_P3) / 100.0,
                    0.01
                )
                s_pump_out_s = s_pump_in
                # CONVERTED-NEEDS UNITS CHECK: removed /1000 Pa->kPa; eeslib uses Pa; removed *1000 kJ/kg->J/kg; eeslib returns J/kg
                h_pump_out_s = fp.enthalpy("water", P=P_P3_out, S=s_pump_out_s)   # J/kg
                # Isentropic energy balance to find ideal work required
                W_dot_pump_s = m_dot_P3 * (h_pump_out_s - h_pump_in)
                W_dot_P3 = W_dot_pump_s / Eta_P3
                # Actual Energy balance to find actual pump enthalpy out
                h_P3_out = (m_dot_P3 * h_pump_in + W_dot_P3) / m_dot_P3

                # Call Outputs
                self.m_dot_P3.v = m_dot_P3       # Mass flow rate out of pump 2
                self.W_dot_P3.v = W_dot_P3        # Power input to pump 2
                self.Eta_P3.v = Eta_P3             # Efficiency of Pump 2
                self.P3_point_1x.v = P3_point_1x  # Previous Flow Rate
                self.P3_point_1y.v = P3_point_1y  # Error associated with previous flow rate
                self.P3_point_2x.v = P3_point_2x  # Flow Rate sent out during this iteration

            else:
                P3_point_1x = self.P3_point_1x.v
                P3_point_1y = self.P3_point_1y.v
                P3_point_2x = self.P3_point_2x.v

                # Compute head loss in system
                # Previous pump head used in last iteration
                P3_pump_head = P_pump_prev / rho_tank_f / 9.81 - P_tank / rho_tank_f / 9.81 - L_tank - self.Length_tank2pump.v
                System_headloss = (P_pump_prev - self.P_piping_sys.v) / rho_tank_f / 9.81
                # Previous error
                P3_point_2y = P3_pump_head - self.P_SD.v / rho_tank_f / 9.81 - System_headloss + P_tank / rho_tank_f / 9.81 + L_tank + self.Length_tank2pump.v

                # Find maximum flow rate allowed for input pump speed
                Q_max_P3 = max(
                    (-B + math.sqrt(B ** 2.0 - 4.0 * A * C)) / (2.0 * A),
                    (-B - math.sqrt(B ** 2.0 - 4.0 * A * C)) / (2.0 * A)
                )

                # Solve for new flow rate leaving pump 1
                if P3_point_2x == P3_point_1x:
                    m = 0.0
                else:
                    m = (P3_point_2y - P3_point_1y) / (P3_point_2x - P3_point_1x)
                y_int = P3_point_1y - m * P3_point_1x

                if abs(P3_point_2y) <= tol:
                    Q_new_P3 = P3_point_2x   # Keep flow rate the same
                elif m != 0.0:
                    if P3_point_2y >= 0.0:   # If Error is Positive, pump head guess is too high so flow must increase
                        Q_new_P3 = max(-y_int / m, P3_point_2x)   # Ensure new Q_dot does not decrease
                        Q_new_P3 = min(Q_new_P3, Q_max_P3)         # Ensure Q_dot_new is not more than pump can provide based on curve
                        delta_Q = abs(Q_new_P3 - P3_point_2x)
                        Q_new_P3 = P3_point_2x + delta_Q * LR      # Increase Q_dot_new based on the learning rate
                    else:
                        Q_new_P3 = min(-y_int / m, P3_point_2x)   # Ensure new Q_dot does not increase
                        Q_new_P3 = max(Q_new_P3, 0.000001)          # Ensure new Q_dot is not negative
                        delta_Q = abs(Q_new_P3 - P3_point_2x)
                        Q_new_P3 = P3_point_2x - delta_Q * LR      # Increase Q_dot_new based on the learning rate
                else:   # If Error is Negative, pump head guess is too low so flow must decrease
                    if P3_point_2y >= 0.0:   # Error is Positive and slope is equal to zero, pump head is too high so flow must increase
                        Q_new_P3 = min(P3_point_2x + 0.0001, Q_max_P3)
                    else:   # Error is Negative and slope is equal to zero, pump head is too low so flow must decrease
                        Q_new_P3 = max(P3_point_2x - 0.0001, 0.000001)

                # Pressure into the pump is the pressure at the bottom of the condenser tank
                P_pump_in = P_tank + rho_tank_f * 9.81 * (L_tank + self.Length_tank2pump.v)
                # CONVERTED-NEEDS UNITS CHECK: removed /1000 Pa->kPa; eeslib uses Pa; removed *1000 kJ/kg->J/kg; eeslib returns J/kg
                h_pump_in = fp.enthalpy("water", P=P_pump_in, T=T_tank)   # J/kg
                s_pump_in = fp.entropy("water", P=P_pump_in, T=T_tank)    # J/(kg·K)
                m_dot_P3 = max(Q_new_P3 * rho_tank_f, 0.0001)
                P_P3_out = (A * Q_new_P3 ** 2.0 + B * Q_new_P3 + C) * rho_tank_f * 9.81 + P_pump_in

                # Solve for pump efficiency, power and enthalpy leaving
                Eta_P3 = max(
                    self.P3_Eta_A.v * Q_new_P3 ** 4.0
                    + self.P3_Eta_B.v * Q_new_P3 ** 3.0
                    + self.P3_Eta_C.v * Q_new_P3 ** 2.0
                    + self.P3_Eta_D.v * Q_new_P3,
                    0.01
                )
                s_pump_out_s = s_pump_in
                # CONVERTED-NEEDS UNITS CHECK: removed /1000 Pa->kPa; eeslib uses Pa; removed *1000 kJ/kg->J/kg; eeslib returns J/kg
                h_pump_out_s = fp.enthalpy("water", P=P_P3_out, S=s_pump_out_s)   # J/kg
                # Isentropic energy balance to find ideal work required
                W_dot_pump_s = m_dot_P3 * (h_pump_out_s - h_pump_in)
                W_dot_P3 = W_dot_pump_s / Eta_P3
                # Actual Energy balance to find actual pump enthalpy out
                h_P3_out = (m_dot_P3 * h_pump_in + W_dot_P3) / m_dot_P3

                P3_point_1x = P3_point_2x   # Saving 2nd point as 1st point for next iteration
                P3_point_1y = P3_point_2y   # Saving 2nd point as 1st point for next iteration
                P3_point_2x = Q_new_P3       # Saving new flow as 2nd point

                # Call Outputs
                self.m_dot_P3.v = m_dot_P3       # Mass flow rate out of pump 1
                self.W_dot_P3.v = W_dot_P3        # Power input to pump 1
                self.Eta_P3.v = Eta_P3             # Efficiency of Pump 1
                self.P3_point_1x.v = P3_point_1x  # Previous Flow Rate
                self.P3_point_1y.v = P3_point_1y  # Error associated with previous flow rate
                self.P3_point_2x.v = P3_point_2x  # Flow Rate sent out during this iteration

        else:
            m_dot_P3 = 0.0
            W_dot_P3 = 0.0
            Eta_P3 = 0.0
            h_P3_out = 0.0
            self.m_dot_P3.v = m_dot_P3   # Mass flow rate out of pump 1
            self.W_dot_P3.v = W_dot_P3   # Power input to pump 1
            self.Eta_P3.v = Eta_P3        # Efficiency of Pump 1

        # -----------------------------------------------------------------------
        # Step 4: Combine Pump Flows

        # total mass flow leaving all the pumps
        m_dot_pump = self.m_dot_P1.v + self.m_dot_P2.v + self.m_dot_P3.v
        if m_dot_pump != 0.0:
            Vol_dot_pump = m_dot_pump / 1000.0
            h_pump_out_combined = (self.m_dot_P1.v * h_P1_out + self.m_dot_P2.v * h_P2_out + self.m_dot_P3.v * h_P3_out) / m_dot_pump
            P_pump_out_combined = (self.m_dot_P1.v * P_P1_out + self.m_dot_P2.v * P_P2_out + self.m_dot_P3.v * P_P3_out) / m_dot_pump   # weighted average
            # CONVERTED-NEEDS UNITS CHECK: removed /1000 Pa->kPa; eeslib uses Pa; removed /1000 J/kg->kJ/kg; eeslib uses J/kg
            T_pump_out_combined = fp.temperature("water", P=P_pump_out_combined, H=h_pump_out_combined)
        else:
            m_dot_pump = 0.00000001
            Vol_dot_pump = 0.0
            h_pump_out_combined = self.h_pump_out.v
            P_pump_out_combined = self.P_pump_out.v
            T_pump_out_combined = self.T_pump_out.v

        W_dot_total = W_dot_P1 + W_dot_P2 + W_dot_P3

        self.m_dot_pump_out.v = m_dot_pump               # mass flow rate leaving all 3 pumps
        self.Vol_dot_pump.v = Vol_dot_pump                # Volumetric Flow Rate leaving the pumps
        self.P_pump_out.v = P_pump_out_combined           # Pressure leaving the pumps
        self.h_pump_out.v = h_pump_out_combined           # enthalpy leaving the pumps
        self.T_pump_out.v = T_pump_out_combined           # Temperature leaving the pumps
        self.W_dot_total.v = W_dot_total                  # Total Pump Power needed

        # -----------------------------------------------------------------------
        # Step 5: Solve for Bleed Entering From LP Turbine

        # mass of the tank at the beginning of the timestep
        m_tank = self.m_tank.v
        # pressure of the tank at the beginning of the timestep
        P_tank = self.P_tank.v
        # enthalpy of the tank at the beginning of the timestep
        h_tank = self.h_tank.v

        if self.Turbine_ON.v == 1.0:
            # previous bleed requested
            m_dot_LPB1_prev = self.m_dot_LPB1.v

            # enthalpy of saturated liquid in the tank for current pressure
            # CONVERTED-NEEDS UNITS CHECK: removed /1000 Pa->kPa; eeslib uses Pa; removed *1000 kJ/kg->J/kg; eeslib returns J/kg
            h_in = fp.enthalpy("water", P=P_tank, Q=0.05)   # J/kg

            # Find the mass of bleed needed to get the water entering to a saturated state
            m_dot_LPB1 = (self.m_dot_fw_in.v * (self.h_fw_in.v - h_in) + self.m_dot_TB.v * (self.h_TB.v - h_in) + self.m_dot_HPFWH.v * (self.h_HPFWH.v - h_in)) / (h_in - self.h_LPB1.v)

            if m_dot_LPB1 > self.m_dot_LPB1_max.v:
                m_dot_LPB1 = self.m_dot_LPB1_max.v   # make sure mass flow entering is not higher than maximum
            elif m_dot_LPB1 < 0.0:
                m_dot_LPB1 = 0.0   # make sure that it cannot request a negative mass

            delta_m_LPB1 = abs(m_dot_LPB1 - m_dot_LPB1_prev)   # change in extraction flow
            if self.model.time != self.model.timestep:
                if delta_m_LPB1 > self.extraction_tol.v:   # change in flow rate is greater than extraction tolerance
                    if delta_m_LPB1 > 0.0:   # bleed wants to increase
                        m_dot_LPB1 = min(m_dot_LPB1_prev + delta_m_LPB1 * 0.1, m_dot_LPB1_prev + self.DA_ss_LPB1.v)
                        m_dot_LPB1 = min(m_dot_LPB1, self.m_dot_LPB1_max.v)
                    else:
                        m_dot_LPB1 = max(m_dot_LPB1_prev - delta_m_LPB1 * 0.1, m_dot_LPB1_prev - self.DA_ss_LPB1.v)
                        m_dot_LPB1 = max(m_dot_LPB1, 0.0)
                else:
                    m_dot_LPB1 = m_dot_LPB1_prev   # keep extraction the same to help with iteration
            else:
                if self.model.iteration != 0:
                    if delta_m_LPB1 > self.extraction_tol.v:   # change in flow rate is greater than extraction tolerance
                        if delta_m_LPB1 > 0.0:   # bleed wants to increase
                            m_dot_LPB1 = min(m_dot_LPB1_prev + delta_m_LPB1 * 0.1, m_dot_LPB1_prev + self.DA_ss_LPB1.v)
                            m_dot_LPB1 = min(m_dot_LPB1, self.m_dot_LPB1_max.v)
                        else:
                            m_dot_LPB1 = max(m_dot_LPB1_prev - delta_m_LPB1 * 0.1, m_dot_LPB1_prev - self.DA_ss_LPB1.v)
                            m_dot_LPB1 = max(m_dot_LPB1, 0.0)
                    else:
                        m_dot_LPB1 = m_dot_LPB1_prev   # keep extraction the same to help with iteration

        else:   # Turbine has not been turned on, no extractions are leaving the turbine
            m_dot_LPB1 = 0.0

        # -----------------------------------------------------------------------
        # Step 6: Solve for Steam Venting Out of Deaerator

        # Initial mixing calculations
        # mass balance
        # total mass entering the tank
        m_dot_in = self.m_dot_fw_in.v + self.m_dot_TB.v + self.m_dot_HPFWH.v + m_dot_LPB1
        # energy balance
        h_in = (self.m_dot_fw_in.v * self.h_fw_in.v + self.m_dot_TB.v * self.h_TB.v + self.m_dot_HPFWH.v * self.h_HPFWH.v + m_dot_LPB1 * self.h_LPB1.v) / m_dot_in
        # CONVERTED-NEEDS UNITS CHECK: removed /1000 Pa->kPa; eeslib uses Pa; removed /1000 J/kg->kJ/kg; eeslib uses J/kg
        T_in = fp.temperature("water", P=P_tank, H=h_in)

        # amount of steam leaving the deaerator
        # CONVERTED-NEEDS UNITS CHECK: removed /1000 Pa->kPa; eeslib uses Pa; removed *1000 kJ/kg->J/kg; eeslib returns J/kg
        h_vent = fp.enthalpy("water", P=P_tank, Q=1.0)   # J/kg

        # -----------------------------------------------------------------------
        # Step 7: Tank Calculations Finding Enthalpy and Pressure Next Timestep

        # converting timestep from hr to s
        ts = self.model.timestep * 3600.0
        dh = 1000.0
        dP = 1000.0
        R_tank = self.D_tank.v / 2.0

        # aa calculations
        drhodhcp_a = drhodhcp(P_tank=P_tank, h_tank=h_tank, dh=dh)
        drhodpch_a = drhodpch(P_tank=P_tank, h_tank=h_tank, dP=dP)
        dudhcp_a = dudhcp(P_tank=P_tank, h_tank=h_tank, dh=dh)
        dudpch_a = dudpch(P_tank=P_tank, h_tank=h_tank, dP=dP)
        rho_tank = m_tank / Vol_tank
        # CONVERTED-NEEDS UNITS CHECK: removed /1000 Pa->kPa; eeslib uses Pa; removed *1000 kJ/kg->J/kg; eeslib returns J/kg
        u_tank = fp.internalenergy("water", P=P_tank, H=h_tank)   # J/kg
        # CONVERTED-NEEDS UNITS CHECK: removed /1000 Pa->kPa; eeslib uses Pa; removed *1000 kJ/kg->J/kg; eeslib returns J/kg
        h_sat_f = fp.enthalpy("water", P=P_tank, Q=0.0)   # J/kg
        denominator = dudhcp_a * drhodpch_a - dudpch_a * drhodhcp_a
        denominator = math.copysign(max(abs(denominator), 0.000025), denominator)
        drhodhcp_a = math.copysign(max(abs(drhodhcp_a), 0.0001), drhodhcp_a)
        tracker_1 = denominator
        tracker_5 = drhodhcp_a
        dPdt_aa = (((u_tank - h_in) * m_dot_in + (-u_tank + h_sat_f) * m_dot_pump - self.m_dot_vent_frac.v * (u_tank - h_vent)) * drhodhcp_a + rho_tank * dudhcp_a * (m_dot_in - m_dot_pump - self.m_dot_vent_frac.v)) / (m_tank * (dudhcp_a * drhodpch_a - dudpch_a * drhodhcp_a))
        dhdt_aa = ((m_dot_in - m_dot_pump - self.m_dot_vent_frac.v) / Vol_tank - drhodpch_a * dPdt_aa) / drhodhcp_a
        P_aa = P_tank + dPdt_aa * ts / 2.0
        h_aa = h_tank + dhdt_aa * ts / 2.0

        # bb calculations
        drhodhcp_a = drhodhcp(P_tank=P_aa, h_tank=h_aa, dh=dh)
        drhodpch_a = drhodpch(P_tank=P_aa, h_tank=h_aa, dP=dP)
        dudhcp_a = dudhcp(P_tank=P_aa, h_tank=h_aa, dh=dh)
        dudpch_a = dudpch(P_tank=P_aa, h_tank=h_aa, dP=dP)
        rho_tank = m_tank / Vol_tank
        # CONVERTED-NEEDS UNITS CHECK: removed /1000 Pa->kPa; eeslib uses Pa; removed *1000 kJ/kg->J/kg; eeslib returns J/kg
        u_tank = fp.internalenergy("water", P=P_aa, H=h_aa)   # J/kg
        T_tank = fp.temperature("water", P=P_aa, H=h_aa)
        # CONVERTED-NEEDS UNITS CHECK: removed /1000 Pa->kPa; eeslib uses Pa; removed *1000 kJ/kg->J/kg; eeslib returns J/kg
        h_sat_f = fp.enthalpy("water", P=P_aa, Q=0.0)   # J/kg
        P_pump_in = P_aa + L_tank * rho_tank_f * 9.81
        denominator = dudhcp_a * drhodpch_a - dudpch_a * drhodhcp_a
        denominator = math.copysign(max(abs(denominator), 0.000025), denominator)
        drhodhcp_a = math.copysign(max(abs(drhodhcp_a), 0.0001), drhodhcp_a)
        tracker_2 = denominator
        tracker_6 = drhodhcp_a
        dPdt_bb = (((u_tank - h_in) * m_dot_in + (-u_tank + h_sat_f) * m_dot_pump - self.m_dot_vent_frac.v * (u_tank - h_vent)) * drhodhcp_a + rho_tank * dudhcp_a * (m_dot_in - m_dot_pump - self.m_dot_vent_frac.v)) / (m_tank * (dudhcp_a * drhodpch_a - dudpch_a * drhodhcp_a))
        dhdt_bb = ((m_dot_in - m_dot_pump - self.m_dot_vent_frac.v) / Vol_tank - drhodpch_a * dPdt_bb) / drhodhcp_a
        P_bb = P_tank + dPdt_bb * ts / 2.0
        h_bb = h_tank + dhdt_bb * ts / 2.0

        # cc calculations
        drhodhcp_a = drhodhcp(P_tank=P_bb, h_tank=h_bb, dh=dh)
        drhodpch_a = drhodpch(P_tank=P_bb, h_tank=h_bb, dP=dP)
        dudhcp_a = dudhcp(P_tank=P_bb, h_tank=h_bb, dh=dh)
        dudpch_a = dudpch(P_tank=P_bb, h_tank=h_bb, dP=dP)
        rho_tank = m_tank / Vol_tank
        # CONVERTED-NEEDS UNITS CHECK: removed /1000 Pa->kPa; eeslib uses Pa; removed *1000 kJ/kg->J/kg; eeslib returns J/kg
        u_tank = fp.internalenergy("water", P=P_bb, H=h_bb)   # J/kg
        T_tank = fp.temperature("water", P=P_bb, H=h_bb)
        # CONVERTED-NEEDS UNITS CHECK: removed /1000 Pa->kPa; eeslib uses Pa; removed *1000 kJ/kg->J/kg; eeslib returns J/kg
        h_sat_f = fp.enthalpy("water", P=P_bb, Q=0.0)   # J/kg
        denominator = dudhcp_a * drhodpch_a - dudpch_a * drhodhcp_a
        denominator = math.copysign(max(abs(denominator), 0.000025), denominator)
        drhodhcp_a = math.copysign(max(abs(drhodhcp_a), 0.0001), drhodhcp_a)
        tracker_3 = denominator
        tracker_7 = drhodhcp_a
        dPdt_cc = (((u_tank - h_in) * m_dot_in + (-u_tank + h_sat_f) * m_dot_pump - self.m_dot_vent_frac.v * (u_tank - h_vent)) * drhodhcp_a + rho_tank * dudhcp_a * (m_dot_in - m_dot_pump - self.m_dot_vent_frac.v)) / (m_tank * (dudhcp_a * drhodpch_a - dudpch_a * drhodhcp_a))
        dhdt_cc = ((m_dot_in - m_dot_pump - self.m_dot_vent_frac.v) / Vol_tank - drhodpch_a * dPdt_cc) / drhodhcp_a
        P_cc = P_tank + dPdt_cc * ts
        h_cc = h_tank + dhdt_cc * ts

        # dd calculations
        drhodhcp_a = drhodhcp(P_tank=P_cc, h_tank=h_cc, dh=dh)
        drhodpch_a = drhodpch(P_tank=P_cc, h_tank=h_cc, dP=dP)
        dudhcp_a = dudhcp(P_tank=P_cc, h_tank=h_cc, dh=dh)
        dudpch_a = dudpch(P_tank=P_cc, h_tank=h_cc, dP=dP)
        rho_tank = m_tank / Vol_tank
        # CONVERTED-NEEDS UNITS CHECK: removed /1000 Pa->kPa; eeslib uses Pa; removed *1000 kJ/kg->J/kg; eeslib returns J/kg
        u_tank = fp.internalenergy("water", P=P_cc, H=h_cc)   # J/kg
        T_tank = fp.temperature("water", P=P_cc, H=h_cc)
        # CONVERTED-NEEDS UNITS CHECK: removed /1000 Pa->kPa; eeslib uses Pa; removed *1000 kJ/kg->J/kg; eeslib returns J/kg
        h_sat_f = fp.enthalpy("water", P=P_cc, Q=0.0)   # J/kg
        denominator = dudhcp_a * drhodpch_a - dudpch_a * drhodhcp_a
        denominator = math.copysign(max(abs(denominator), 0.000025), denominator)
        drhodhcp_a = math.copysign(max(abs(drhodhcp_a), 0.0001), drhodhcp_a)
        tracker_4 = denominator
        tracker_8 = drhodhcp_a
        dPdt_dd = (((u_tank - h_in) * m_dot_in + (-u_tank + h_sat_f) * m_dot_pump - self.m_dot_vent_frac.v * (u_tank - h_vent)) * drhodhcp_a + rho_tank * dudhcp_a * (m_dot_in - m_dot_pump - self.m_dot_vent_frac.v)) / (m_tank * (dudhcp_a * drhodpch_a - dudpch_a * drhodhcp_a))
        dhdt_dd = ((m_dot_in - m_dot_pump - self.m_dot_vent_frac.v) / Vol_tank - drhodpch_a * dPdt_dd) / drhodhcp_a

        # End of timestep Pressure and Enthalpy
        P_tank_new = P_tank + (dPdt_aa + 2.0 * dPdt_bb + 2.0 * dPdt_cc + dPdt_dd) * ts / 6.0
        h_tank_new = h_tank + (dhdt_aa + 2.0 * dhdt_bb + 2.0 * dhdt_cc + dhdt_dd) * ts / 6.0
        # CONVERTED-NEEDS UNITS CHECK: removed /1000 Pa->kPa; eeslib uses Pa; removed /1000 J/kg->kJ/kg; eeslib uses J/kg
        x_tank_new = fp.quality("water", P=P_tank_new, H=h_tank_new)
        rho_tank_new = fp.density("water", P=P_tank_new, H=h_tank_new)
        T_tank_new = fp.temperature("water", P=P_tank_new, H=h_tank_new)

        # -----------------------------------------------------------------------
        # Step 8: Finding new tank level

        L_tank_prev = self.L_tank.v
        dmdt = m_dot_in - m_dot_pump - self.m_dot_vent_frac.v
        m_tank_new = m_tank + dmdt * ts
        m_tank_g_new = m_tank_new * x_tank_new
        m_tank_f_new = m_tank_new - max(m_tank_g_new, 0.0)
        # CONVERTED-NEEDS UNITS CHECK: removed /1000 Pa->kPa; eeslib uses Pa
        rho_f_new = fp.density("water", P=P_tank_new, Q=0.0)
        Vol_liquid = m_tank_f_new / rho_f_new
        level_tol = 0.01   # [m]
        L_tank_new = tank_level(Vol_liquid, self.D_tank.v, self.Length_tank.v, L_tank_prev, level_tol)

        # Set the Output Values
        self.m_tank_new.v = m_tank_new      # total mass in the tank for the next timestep
        self.T_tank_new.v = T_tank_new      # Temperature of the tank for the next timestep
        self.L_tank_new.v = L_tank_new      # Level of the tank for the next timestep
        self.P_tank_new.v = P_tank_new      # Pressure of the tank for the next timestep
        self.h_tank_new.v = h_tank_new      # Enthalpy of the tank for the next timestep
        self.m_dot_vent_out.v = self.m_dot_vent_frac.v   # mass of steam exiting through deaerator vents
        self.h_vent.v = h_vent              # enthalpy of steam exiting through deaerator vents
        self.m_dot_LPB1.v = m_dot_LPB1     # mass flow requested from LP extraction 1
