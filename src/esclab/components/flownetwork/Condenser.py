"""Condenser component model (Type 6007)."""

import math

import numpy as np
from eeslib import fluid_properties as fp

from esclab.simulate import Component
from esclab.components.flownetwork.simple_pipe import FricFactor_IC
from esclab.components.flownetwork.esol6015_helpers import f_cp_water, drhodhcp, drhodpch, dudhcp, dudpch


class Condenser(Component):
    """
    TRNSYS Type 6007: ESOL6007-Condenser.

    Models a condenser hotwell tank with a condensate pump, a cooling-water
    heat exchanger, and a reservoir mass-flow balance.  Three sequential
    calculation steps are performed each iteration:

    1. Pump flow solving (secant-like iteration between system-curve and pump
       curve).
    2. Condenser HX heat-transfer calculation.
    3. Hotwell tank energy/mass balance via RK4 integration with an inner
       iteration to keep the tank level constant via a reservoir flow.

    Parameters
    ----------
    Height_tank : float
        Height of the condenser tank [m].
    Area_tank : float
        Cross-sectional area of the condenser tank [m^2].
    L_tank_ini : float
        Initial height of the water level in the condenser tank [m].
    P_tank_ini : float
        Initial pressure of the condenser tank [Pa].
    Coef_a : float
        Condenser pump curve coefficient (Coef_a*Q^2 + Coef_b*Q + Coef_c = pump head).
    Coef_b : float
        Condenser pump curve coefficient.
    Coef_c : float
        Condenser pump curve coefficient.
    D_pump_out : float
        Diameter of pump outlet [m].
    Eta_A : float
        Condenser pump efficiency curve coefficient (Eta_A*Q^4 + Eta_B*Q^3 + Eta_C*Q^2 + Eta_D*Q).
    Eta_B : float
        Condenser pump efficiency curve coefficient.
    Eta_C : float
        Condenser pump efficiency curve coefficient.
    Eta_D : float
        Condenser pump efficiency curve coefficient.
    TTD : float
        Terminal temperature difference of steam entering the HX [K].
    No_tubes : float
        Number of tubes travelling through condenser HX.
    Length_Tubes : float
        Length of tubes travelling through condenser HX [m].
    ID : float
        Inner diameter of condenser tubes [m].
    Th : float
        Tube thickness [m].

    Inputs
    ------
    m_dot_in : float
        Turbine exhaust entering the condenser [kg/s].
    h_in : float
        Enthalpy of turbine exhaust entering the condenser [J/kg].
    P_in : float
        Pressure of turbine exhaust entering the condenser [Pa].
    m_dot_cool_in : float
        Cooling flow entering the condenser HX [kg/s].
    h_cool_in : float
        Enthalpy of the cooling flow entering the condenser HX [J/kg].
    P_cool_in : float
        Pressure of the cooling flow entering the condenser HX [Pa].
    P_DA : float
        Pressure of the deaerator used to solve for point on pump curve [Pa].
    Sys_losses : float
        System losses between the tanks, used for pump iteration process [Pa].
    pump_speed : float
        Pump speed the condenser pump is running at [MAX == 1].

    Outputs
    -------
    m_dot_pump : float
        Pump mass out [kg/s].
    Vol_pump_out : float
        Volumetric flow rate leaving the condenser pump [m^3/s].
    T_pump_out : float
        Pump temperature out [K].
    h_pump_out : float
        Pump enthalpy out [J/kg].
    P_pump_out : float
        Pump pressure out [Pa].
    m_dot_cool_out : float
        Mass flow of cooling flow exiting condenser [kg/s].
    Vol_cool_out : float
        Volumetric flow rate of cooling flow exiting the condenser [m^3/s].
    T_cool_out : float
        Cooling water temperature exiting the condenser [K].
    P_cool_out : float
        Cooling water pressure exiting the condenser [Pa].
    h_cool_out : float
        Cooling water enthalpy exiting the condenser [J/kg].
    L_tank_next : float
        Tank level for the next timestep [m].
    P_tank_next : float
        Tank pressure for the next timestep [Pa].
    h_tank_next : float
        Tank enthalpy for the next timestep [J/kg].
    T_tank_next : float
        Tank temperature for the next timestep [K].
    m_tank_next : float
        Mass in the tank for the next timestep [kg].
    L_tank_cur : float
        Current tank level for this timestep [m].
    P_tank_cur : float
        Current tank pressure for this timestep [Pa].
    h_tank_cur : float
        Current tank enthalpy for this timestep [J/kg].
    m_tank_cur : float
        Current mass in the tank for this timestep [kg].
    prev_error : float
        Previous pump error for last guessed pump flow rate (Point 1y) [m].
    Q_new_out : float
        Current pump flow rate guessed for the next iteration (Point_2x).
    prev_flow : float
        Previous pump flow rate guessed for the last iteration (Point_1x).
    Eta_pump_out : float
        Pump efficiency [1 = MAX].
    W_dot_pump_out : float
        Pump power required [W].
    m_dot_res_out : float
        Mass flow from reservoir tank [kg/s].
    h_guess_out : float
        Guess exit enthalpy for condenser heat exchanger function.
    Q_dot_hx : float
        Total heat removed by condenser HX [W].
    ff_out : float
        Guess friction factor.
    """

    # *** Model Parameters ***
    # Height of the condenser tank in [m]
    Height_tank = Component.Parameter()
    # Area of the condenser tank in [m^2]
    Area_tank = Component.Parameter()
    # Initial height of the water level in the condenser tank [m]
    L_tank_ini = Component.Parameter()
    # Initial Pressure of the Condenser Tank [Pa]
    P_tank_ini = Component.Parameter()
    # Condenser Pump Curve Coefficient ... Coef_a*Flow^2 + Coef_b*Flow + Coef_c = Pump_head
    Coef_a = Component.Parameter()
    # Condenser Pump Curve Coefficient
    Coef_b = Component.Parameter()
    # Condenser Pump Curve Coefficient
    Coef_c = Component.Parameter()
    # Diameter of pump outlet
    D_pump_out = Component.Parameter()
    # Condenser Pump Efficiency Curve ... Eta_A*Flow^4 + Eta_B*Flow^3 + Eta_C*Flow^2 + Eta_D*Flow
    Eta_A = Component.Parameter()
    # Condenser Pump Efficiency Curve
    Eta_B = Component.Parameter()
    # Condenser Pump Efficiency Curve
    Eta_C = Component.Parameter()
    # Condenser Pump Efficiency Curve
    Eta_D = Component.Parameter()
    # Terminal Temperature Difference of steam entering the HX
    TTD = Component.Parameter()
    # Number of tubes traveling through condenser HX
    No_tubes = Component.Parameter()
    # Length of tubes traveling through condenser HX
    Length_Tubes = Component.Parameter()
    # Inner Diameter of Tubes
    ID = Component.Parameter()
    # Tube Thickness
    Th = Component.Parameter()

    # *** Model Inputs ***
    # Turbine exhaust entering the condenser [kg/s]
    m_dot_in = Component.Input()
    # Enthalpy of turbine exhaust entering the condenser [J/kg]
    h_in = Component.Input()
    # Pressure of turbine exhaust entering the condenser [Pa]
    P_in = Component.Input()
    # Cooling Flow entering the condenser HX [kg/s]
    m_dot_cool_in = Component.Input()
    # Enthalpy of the Cooling Flow entering the condenser HX [J/kg]
    h_cool_in = Component.Input()
    # Pressure of the Cooling Flow entering the condenser HX [Pa]
    P_cool_in = Component.Input()
    # Pressure of the Deaerator used to solve for point on pump curve [Pa]
    # !!!!!Ensure its the pressure for that timestep!!!!!
    P_DA = Component.Input()
    # System Losses between the tanks - used for pump iteration process
    Sys_losses = Component.Input()
    # Pump Speed the condenser pump is running at [MAX == 1]
    pump_speed = Component.Input()

    # *** Model Outputs ***
    # Pump Mass Out [kg/s]
    m_dot_pump = Component.Output()
    # Volumetric Flow Rate of Pump Flow [m^3/s]
    Vol_pump_out = Component.Output()
    # Pump temperature out [K]
    T_pump_out = Component.Output()
    # Pump enthalpy out [J/kg]
    h_pump_out = Component.Output()
    # Pump Pressure Out [Pa]
    P_pump_out = Component.Output()
    # Mass flow of cooling flow exiting condenser [kg/s]
    m_dot_cool_out = Component.Output()
    # Volumetric flow rate of cooling flow exiting the condenser [m^3/s]
    Vol_cool_out = Component.Output()
    # Cooling Water Temperature exiting the condenser [K]
    T_cool_out = Component.Output()
    # Cooling Water Pressure exiting the condenser [Pa]
    P_cool_out = Component.Output()
    # Cooling Water Enthalpy exiting the condenser [J/kg]
    h_cool_out = Component.Output()
    # Tank Level for the next timestep [m]
    L_tank_next = Component.Output()
    # Tank Pressure for the next timestep [Pa]
    P_tank_next = Component.Output()
    # Tank Enthalpy for the next timestep [J/kg]
    h_tank_next = Component.Output()
    # Tank Temperature for the next timestep [K]
    T_tank_next = Component.Output()
    # Mass in the tank for the next timestep [kg]
    m_tank_next = Component.Output()
    # Current Tank Level for this timestep [m]
    L_tank_cur = Component.Output()
    # Current Tank Pressure for this timestep [Pa]
    P_tank_cur = Component.Output()
    # Current Tank Enthalpy for this timestep [J/kg]
    h_tank_cur = Component.Output()
    # Current Mass in the tank for this timestep [kg]
    m_tank_cur = Component.Output()
    # Previous Pump Error for last guessed pump flow rate (Point 1y) [m]
    prev_error = Component.Output()
    # Current Pump Flow rate guessed for the next iteration (Point_2x)
    Q_new_out = Component.Output()
    # Previous Pump flow rate guessed for the last iteration (Point 1x)
    prev_flow = Component.Output()
    # Pump Efficiency [1 = MAX]
    Eta_pump_out = Component.Output()
    # Pump Power [W]
    W_dot_pump_out = Component.Output()
    # Mass flow from reservoir tank
    m_dot_res_out = Component.Output()
    # Guess exit enthalpy for condenser heat exchanger function
    h_guess_out = Component.Output()
    # Total Heat removed by condenser HX [W]
    Q_dot_hx = Component.Output()
    # Guess Friction Factor
    ff_out = Component.Output()

    def calculate(self):
        # Perform Any "After Convergence" Manipulations That May Be Required at the End of Each Timestep
        if self.model.is_converged:
            # Set the new tank level, pressure and enthalpy as start of next timestep values
            L_tank_new = self.L_tank_next.v
            P_tank_new = self.P_tank_next.v
            H_tank_new = self.h_tank_next.v
            m_tank_new = self.m_tank_next.v
            self.L_tank_cur.v = L_tank_new
            self.P_tank_cur.v = P_tank_new
            self.h_tank_cur.v = H_tank_new
            self.m_tank_cur.v = m_tank_new
            return

        # Do All of the First Timestep Manipulations Here - There Are No Iterations at the Initial Time
        if self.model.is_first_step:
            # Parameters
            # Height of the condenser tank in [m]
            # Area of the condenser tank in [m^2]
            # Initial height of the water level in the condenser tank [m]
            # Initial Pressure of the Condenser Tank [Pa]
            P_tank = self.P_tank_ini.v
            # Condenser Pump Curve Coefficient ... Coef_a*Flow^2 + Coef_b*Flow + Coef_c = Pump_head
            # Condenser Pump Efficiency Curve ... Eta_A*Flow^4 + Eta_B*Flow^3 + Eta_C*Flow^2 + Eta_D*Flow
            # Terminal Temperature Difference of steam entering the HX
            # Number of tubes traveling through condenser HX
            # Length of tubes traveling through condenser HX
            # Inner Diameter of Tubes
            # Tube Thickness

            # Inputs
            # Turbine exhaust entering the condenser [kg/s]
            # Enthalpy of turbine exhaust entering the condenser [J/kg]
            # Pressure of turbine exhaust entering the condenser [Pa]
            # Cooling Flow entering the condenser HX [kg/s]
            # Temperature of the Cooling Flow entering the condenser HX [K]
            # Pressure of the Cooling Flow entering the condenser HX [Pa]
            # Pressure of the Deaerator used to solve for point on pump curve [Pa]
            # !!!!!Ensure its the pressure for that timestep!!!!!
            # Pump Speed the condenser pump is running at [MAX == 1]

            # !!!!!!!SOLVE FOR INITIAL TANK ENTHALPY AND MASS!!!!!!!!!!
            # TODO-NEEDS UNITS CHECK: FIT_PQ passes P in kPa (P/1000); eeslib fp may use Pa
            T_tank = fp.temperature("water", P=P_tank / 1000.0, x=1.0)
            rho_tank_g = fp.density("water", P=P_tank / 1000.0, x=1.0)
            # TODO-NEEDS UNITS CHECK: FIT_PQ passes P in kPa; enthalpy returned in kJ/kg then *1000 to J/kg
            h_pump_in = fp.enthalpy("water", P=P_tank / 1000.0, x=0.0)
            rho_f = 1000.0  # assumed incompressible
            h_pump_in = h_pump_in * 1000.0  # converting from kJ/kg to J/kg
            m_tank_f = self.Area_tank.v * self.L_tank_ini.v * rho_f
            m_tank_g = (self.Height_tank.v - self.L_tank_ini.v) * self.Area_tank.v * rho_tank_g
            m_tank = m_tank_f + m_tank_g
            x_tank = m_tank_g / m_tank
            # TODO-NEEDS UNITS CHECK: FIT_PQ passes P in kPa; enthalpy returned in kJ/kg then *1000 to J/kg
            h_tank = fp.enthalpy("water", P=P_tank / 1000.0, x=x_tank)
            h_tank = h_tank * 1000.0  # converting from KJ/kg to J/kg

            # !!!!!!!!!!!!!!SEND OUT PUMP FLOW !!!!!!!!!!!!!!!!!!!!!
            # Pump pressure assuming no flow
            P_pump_out = self.Coef_c.v * self.pump_speed.v ** 2 * rho_f * 9.81 + P_tank + self.L_tank_ini.v * rho_f * 9.81
            T_pump_out = T_tank
            # TODO-NEEDS UNITS CHECK: FIT_TP passes P in kPa; enthalpy returned in kJ/kg then *1000 to J/kg
            h_pump_out = fp.enthalpy("water", P=P_pump_out / 1000.0, T=T_tank)
            h_pump_out = h_pump_out * 1000.0  # converting from kJ/kg to J/kg

            # Calculating Volumetric flows
            Vol_pump_out = 0.0
            Vol_cool_out = self.m_dot_cool_in.v / 1000.0

            # Set the Initial Values of the Outputs (#,Value)
            self.m_dot_pump.v = 0.000001      # Pump Mass Out [kg/s]
            self.Vol_pump_out.v = Vol_pump_out  # Volumetric Flow Rate of Pump Flow [m^3/s]
            self.T_pump_out.v = T_pump_out      # Pump temperature out [K]
            self.h_pump_out.v = h_pump_out      # Pump enthalpy out [J/kg]
            self.P_pump_out.v = P_pump_out      # Pump Pressure Out [Pa]
            self.m_dot_cool_out.v = self.m_dot_cool_in.v  # Mass flow of cooling flow exiting condenser [kg/s]
            self.Vol_cool_out.v = Vol_cool_out  # Volumetric flow rate of cooling flow exiting the condenser [m^3/s]
            self.T_cool_out.v = 300.0           # Cooling Water Temperature exiting the condenser [K]
            self.P_cool_out.v = self.P_cool_in.v  # Cooling Water Pressure exiting the condenser [Pa]
            self.h_cool_out.v = self.h_cool_in.v  # Cooling Water Enthalpy exiting the condenser [J/kg]
            self.L_tank_next.v = self.L_tank_ini.v  # Tank Level for the next timestep [m]
            self.P_tank_next.v = P_tank           # Tank Pressure for the next timestep [Pa]
            self.h_tank_next.v = h_tank           # Tank Enthalpy for the next timestep [J/kg]
            self.T_tank_next.v = T_tank           # Tank Temperature for the next timestep [K]
            self.m_tank_next.v = m_tank           # Mass in the tank for the next timestep [kg]
            self.L_tank_cur.v = self.L_tank_ini.v  # Current Tank Level for this timestep [m]
            self.P_tank_cur.v = P_tank            # Current Tank Pressure for this timestep [Pa]
            self.h_tank_cur.v = h_tank            # Current Tank Enthalpy for this timestep [J/kg]
            self.m_tank_cur.v = m_tank            # Current Mass in the tank for this timestep [kg]
            self.prev_error.v = 0.0               # Previous Pump Error for last guessed pump flow rate (Point 1y) [m]
            self.Q_new_out.v = 0.0                # Current Pump Flow rate guessed for the next iteration (Point_2x)
            self.prev_flow.v = 0.0                # Previous Pump flow rate guessed for the last iteration (Point 1x)
            self.Eta_pump_out.v = 0.0             # Pump Efficiency [1 = MAX]
            self.W_dot_pump_out.v = 0.0           # Pump Power [W]
            self.m_dot_res_out.v = 0.0            # Mass flow from reservoir tank
            self.h_guess_out.v = 0.0              # Guess exit enthalpy for condenser heat exchanger function
            self.Q_dot_hx.v = 0.0                 # Total Heat removed by condenser HX [W]
            self.ff_out.v = 0.1                   # Guess Friction Factor
            return

        # Read the Parameters
        # Height of the condenser tank in [m]
        # Area of the condenser tank in [m^2]
        # Initial height of the water level in the condenser tank [m]
        # Initial Pressure of the Condenser Tank [Pa]
        # Condenser Pump Curve Coefficient ... Coef_a*Flow^2 + Coef_b*Flow + Coef_c = Pump_head
        # Condenser Pump Curve Coefficient
        # Condenser Pump Curve Coefficient
        # Diameter of pump outlet
        # Condenser Pump Efficiency Curve ... Eta_A*Flow^4 + Eta_B*Flow^3 + Eta_C*Flow^2 + Eta_D*Flow
        # Condenser Pump Efficiency Curve
        # Condenser Pump Efficiency Curve
        # Condenser Pump Efficiency Curve
        # Terminal Temperature Difference of steam entering the HX
        # Number of tubes traveling through condenser HX
        # Length of tubes traveling through condenser HX
        # Inner Diameter of Tubes
        # Tube Thickness

        # Inputs
        # Turbine exhaust entering the condenser [kg/s]
        # Enthalpy of turbine exhaust entering the condenser [J/kg]
        # Pressure of turbine exhaust entering the condenser [Pa]
        # Cooling Flow entering the condenser HX [kg/s]
        # Enthalpy of the Cooling Flow entering the condenser HX [J/kg]
        # Pressure of the Cooling Flow entering the condenser HX [Pa]
        # Pressure of the Deaerator used to solve for point on pump curve [Pa]
        # !!!!!Ensure its the pressure for that timestep!!!!!
        # System Losses between the tanks - used for pump iteration process
        # Pump Speed the condenser pump is running at [MAX == 1]

        # !!!!!!!!!!!!!!! STEP ONE SOLVE FOR MASS FLOW LEAVING THE PUMPS BASED ON SYSTEM LOSSES AND PUMP CURVE!!!!!!!!!!!!!!!!!!!!!
        tol = 1.0  # Error within 1 kPa is accepted
        A = self.Coef_a.v
        B = self.Coef_b.v * self.pump_speed.v
        C = self.Coef_c.v * (self.pump_speed.v ** 2.0)
        rho_f = 1000.0

        # !!!!!!!!!!!!!!!! STEP 1: Solving for flow leaving the Condenser Pumps !!!!!!!!!!!!!!!!!!!!!!!

        if self.model.timestep_iteration == 0:
            # First Iteration of Timestep, adjust pump flow slightly if needed

            m_dot_pump = self.m_dot_pump.v  # Previous flow leaving the pump
            P_tank = self.P_tank_cur.v      # Pressure of tank at this timestep
            L_tank = self.L_tank_cur.v      # Level of tank at this timestep
            P_pump_prev = self.h_pump_out.v  # Pump outlet pressure computed at last iteration
            # NOTE: Fortran getOutputValue(4) is h_pump_out (output 4), but comment says "Pump outlet pressure at last
            # iteration". The variable P_pump_prev here is loaded from output 4 (h_pump_out). This appears to be a
            # comment/variable name mismatch in the original Fortran.
            # TODO-NEEDS CONVERSION REVIEW: getOutputValue(4) is h_pump_out slot but used as P_pump_prev

            # Find error for previous flow rate
            # TODO-NEEDS UNITS CHECK: FIT_PQ passes P in kPa
            T_tank = fp.temperature("water", P=P_tank / 1000.0, x=0.0)

            Point_1x = m_dot_pump / rho_f
            Pump_head = P_pump_prev / rho_f / 9.81 - P_tank / rho_f / 9.81 - L_tank
            System_headloss = (P_pump_prev - self.Sys_losses.v) / rho_f / 9.81
            Point_1y = Pump_head - self.P_DA.v / rho_f / 9.81 - System_headloss + P_tank / rho_f / 9.81 + L_tank

            # Maximum flow allowed for given pump speed input
            Q_max = max(
                (-B + math.sqrt(B ** 2.0 - 4.0 * A * C)) / (2.0 * A),
                (-B - math.sqrt(B ** 2.0 - 4.0 * A * C)) / (2.0 * A),
            )

            if abs(Point_1y) <= tol:
                Q_new = Point_1x
            elif Point_1y > 0.0:  # increase flow rate
                Q_new = min(Point_1x + 0.01, Q_max)  # Move to maximum flow pump can provide at given speed
            elif Point_1y < 0.0:  # decrease flow rate
                Q_new = max(Point_1x - 0.01, 0.000000001)  # Move to minimum flow pump can provide
            else:
                Q_new = Point_1x  # Keep Q at same point

            P_pump_in = P_tank + rho_f * 9.81 * L_tank
            # TODO-NEEDS UNITS CHECK: FIT_TP passes P in kPa; enthalpy in kJ/kg then *1000; entropy in kJ/kg-K then *1000
            h_pump_in = fp.enthalpy("water", P=P_pump_in / 1000.0, T=T_tank)
            s_pump_in = fp.entropy("water", P=P_pump_in / 1000.0, T=T_tank)
            h_pump_in = h_pump_in * 1000.0  # converting from kJ/kg to J/kg
            s_pump_in = s_pump_in * 1000.0  # converting from KJ/kg-K to J/kg-K
            m_dot_pump = Q_new * rho_f
            P_pump_out = (A * Q_new ** 2.0 + B * Q_new + C) * rho_f * 9.81 + P_pump_in

            pump_speed_safe = max(self.pump_speed.v, 0.00000001)
            Eta_pump = max(
                self.Eta_A.v * (Q_new / pump_speed_safe) ** 4.0
                + self.Eta_B.v * (Q_new / pump_speed_safe) ** 3.0
                + self.Eta_C.v * (Q_new / pump_speed_safe) ** 2.0
                + self.Eta_D.v * Q_new / pump_speed_safe,
                0.01,
            )

            s_pump_out_s = s_pump_in
            # TODO-NEEDS UNITS CHECK: FIT_PS passes P in kPa, S in kJ/kg-K; enthalpy returned in kJ/kg then *1000
            h_pump_out_s = fp.enthalpy("water", P=P_pump_out / 1000.0, s=s_pump_out_s / 1000.0)
            h_pump_out_s = h_pump_out_s * 1000.0  # converting from kJ/kg to J/kg
            W_dot_pump_s = m_dot_pump * (h_pump_out_s - h_pump_in)  # Isentropic energy balance to find ideal work required
            W_dot_pump = W_dot_pump_s / Eta_pump
            h_pump_out = (m_dot_pump * h_pump_in + W_dot_pump) / m_dot_pump  # Actual Energy balance to find actual pump enthalpy out
            # TODO-NEEDS UNITS CHECK: FIT_PH passes P in kPa, h in kJ/kg; specific volume returned
            T_pump_out = fp.temperature("water", P=P_pump_out / 1000.0, h=h_pump_out / 1000.0)
            v_pump_out = fp.volume("water", P=P_pump_out / 1000.0, h=h_pump_out / 1000.0)

            Vol_pump_out = m_dot_pump * v_pump_out

            self.m_dot_pump.v = m_dot_pump          # Pump Mass Out [kg/s]
            self.Vol_pump_out.v = Vol_pump_out       # Volumetric Flow Rate leaving the condenser pump [m^3/s]
            self.T_pump_out.v = T_pump_out           # Pump temperature out [K]
            self.h_pump_out.v = h_pump_out           # Pump enthalpy out [J/kg]
            self.P_pump_out.v = P_pump_out           # Pump Pressure Out [Pa]
            self.prev_error.v = Point_1y             # prev_error (Point 1y)
            self.Q_new_out.v = Q_new                 # New flow rate (point 2x)
            self.prev_flow.v = Point_1x              # Prev flow rate (point 1x)
            self.Eta_pump_out.v = Eta_pump           # Pump Efficiency
            self.W_dot_pump_out.v = W_dot_pump       # Pump Power Required

        else:
            # Previous iteration/timestep information needed
            Point_1y = self.prev_error.v
            Point_2x = self.Q_new_out.v
            Point_1x = self.prev_flow.v
            P_tank = self.P_tank_cur.v          # Pressure of tank at beginning of timestep
            L_tank = self.L_tank_cur.v          # Level of tank at beginning of timestep
            P_pump_prev = self.P_pump_out.v     # Pump outlet pressure computed at last iteration

            # Compute head loss in system
            # TODO-NEEDS UNITS CHECK: FIT_PQ passes P in kPa
            T_tank = fp.temperature("water", P=P_tank / 1000.0, x=0.0)
            Pump_head = P_pump_prev / rho_f / 9.81 - P_tank / rho_f / 9.81 - L_tank
            System_headloss = (P_pump_prev - self.Sys_losses.v) / rho_f / 9.81
            Point_2y = Pump_head - self.P_DA.v / rho_f / 9.81 - System_headloss + P_tank / rho_f / 9.81 + L_tank

            # Find maximum flow rate allowed for input pump speed
            Q_max = max(
                (-B + math.sqrt(B ** 2.0 - 4.0 * A * C)) / (2.0 * A),
                (-B - math.sqrt(B ** 2.0 - 4.0 * A * C)) / (2.0 * A),
            )

            # Create line between points
            if Point_2x == Point_1x:
                m = 0.0
            else:
                m = (Point_2y - Point_1y) / (Point_2x - Point_1x)
            y_int = Point_1y - m * Point_1x

            if abs(Point_2y) <= tol:
                Q_new = Point_2x
            elif m != 0.0:
                if Point_2y >= 0.0:
                    Q_new = max(-y_int / m, Point_2x)   # Ensure new Q_dot does not decrease
                    Q_new = min(Q_new, Q_max)            # Ensure Q_dot_new is not more than pump can provide based on curve
                else:
                    Q_new = min(-y_int / m, Point_2x)   # Ensure new Q_dot does not increase
                    Q_new = max(Q_new, 0.000000001)      # Ensure new Q_dot is not negative
            else:
                if Point_2y >= 0:  # Check new point slightly to the right
                    Q_new = min(Point_2x + 0.001, Q_max)
                else:              # Check new point slightly to the left
                    Q_new = max(Point_2x - 0.001, 0.000000001)

            P_pump_in = P_tank + rho_f * 9.81 * L_tank  # Pressure into the pump is the pressure at the bottom of the condenser tank
            # TODO-NEEDS UNITS CHECK: FIT_TP passes P in kPa; enthalpy in kJ/kg then *1000; entropy in kJ/kg-K then *1000
            h_pump_in = fp.enthalpy("water", P=P_pump_in / 1000.0, T=T_tank)
            s_pump_in = fp.entropy("water", P=P_pump_in / 1000.0, T=T_tank)
            h_pump_in = h_pump_in * 1000.0  # converting from kJ/kg to J/kg
            s_pump_in = s_pump_in * 1000.0  # converting from KJ/kg-K to J/kg-K
            m_dot_pump = Q_new * rho_f
            P_pump_out = (A * Q_new ** 2.0 + B * Q_new + C) * rho_f * 9.81 + P_pump_in

            pump_speed_safe = max(self.pump_speed.v, 0.00000001)
            Eta_pump = max(
                self.Eta_A.v * (Q_new / pump_speed_safe) ** 4.0
                + self.Eta_B.v * (Q_new / pump_speed_safe) ** 3.0
                + self.Eta_C.v * (Q_new / pump_speed_safe) ** 2.0
                + self.Eta_D.v * Q_new / pump_speed_safe,
                0.01,
            )

            s_pump_out_s = s_pump_in
            # TODO-NEEDS UNITS CHECK: FIT_PS passes P in kPa, S in kJ/kg-K; enthalpy returned in kJ/kg then *1000
            h_pump_out_s = fp.enthalpy("water", P=P_pump_out / 1000.0, s=s_pump_out_s / 1000.0)
            h_pump_out_s = h_pump_out_s * 1000.0  # converting from kJ/kg to J/kg
            W_dot_pump_s = m_dot_pump * (h_pump_out_s - h_pump_in)  # Isentropic energy balance to find ideal work required
            W_dot_pump = W_dot_pump_s / Eta_pump
            h_pump_out = (m_dot_pump * h_pump_in + W_dot_pump) / m_dot_pump  # Actual Energy balance to find actual pump enthalpy out
            # TODO-NEEDS UNITS CHECK: FIT_PH passes P in kPa, h in kJ/kg; specific volume returned
            T_pump_out = fp.temperature("water", P=P_pump_out / 1000.0, h=h_pump_out / 1000.0)
            v_pump_out = fp.volume("water", P=P_pump_out / 1000.0, h=h_pump_out / 1000.0)

            Point_1x = Point_2x  # Saving 2nd point as 1st point for next iteration
            Point_1y = Point_2y  # Saving 2nd point as 1st point for next iteration
            Point_2x = Q_new     # Saving new flow as 2nd point

            Vol_pump_out = m_dot_pump * v_pump_out

            # Set the Outputs from this Model (#,Value)
            self.m_dot_pump.v = m_dot_pump          # Pump Mass Out [kg/s]
            self.Vol_pump_out.v = Vol_pump_out       # Volumetric Flow Rate leaving the condenser pump [m^3/s]
            self.T_pump_out.v = T_pump_out           # Pump temperature out [K]
            self.h_pump_out.v = h_pump_out           # Pump enthalpy out [J/kg]
            self.P_pump_out.v = P_pump_out           # Pump Pressure Out [Pa]
            self.prev_error.v = Point_1y             # prev_error (Point 1y)
            self.Q_new_out.v = Point_2x              # New flow rate (point 2x)
            self.prev_flow.v = Point_1x              # Prev flow rate (point 1x)
            self.Eta_pump_out.v = Eta_pump           # Pump Efficiency [MAX == 1]
            self.W_dot_pump_out.v = W_dot_pump       # Pump Power Required [W]

        # After either branch the following locals must be available; read back from outputs for the subsequent steps
        m_dot_pump = self.m_dot_pump.v
        P_pump_out = self.P_pump_out.v
        h_pump_out = self.h_pump_out.v
        rho_f = 1000.0

        # !!!!!!!!!!! STEP 2: CONDENSER HX CALCULATIONS !!!!!!!!!!!!!!!!!!!!!!!!!!!

        P_tank_prev = self.P_tank_cur.v
        h_tank_prev = self.h_tank_cur.v
        ff = self.ff_out.v
        OD = self.ID.v + 2.0 * self.Th.v
        A_s_inner = 3.14 * self.ID.v * self.Length_Tubes.v * self.No_tubes.v
        A_s_outer = 3.14 * OD * self.Length_Tubes.v * self.No_tubes.v

        if self.m_dot_cool_in.v > 1.0:
            # Calculate the temperature of the cooling water out based on the TTD
            # TODO-NEEDS UNITS CHECK: FIT_PH passes P in kPa, h in kJ/kg
            T_cool_in = fp.temperature("water", P=self.P_cool_in.v / 1000.0, h=self.h_cool_in.v / 1000.0)
            # TODO-NEEDS UNITS CHECK: FIT_PQ passes P in kPa
            T_sat = fp.temperature("water", P=P_tank_prev / 1000.0, x=1.0)
            T_cool_out = max(T_sat - self.TTD.v, T_cool_in)
            T_cool_out = min(T_cool_out, 342.0)  # Cooling water temperature can not be higher than 70 [C] otherwise Cooling Tower Type will Fail

            # Calculate the Resistance on inside of the condenser tubes from cooling water
            vel = self.m_dot_cool_in.v / self.No_tubes.v / (3.14 / 4.0 * self.ID.v ** 2.0) / 1000.0
            # TODO-NEEDS UNITS CHECK: FIT_PH passes P in kPa, h in kJ/kg; viscosity in microPa-s, conductivity returned
            mu_water = fp.viscosity("water", P=self.P_cool_in.v / 1000.0, h=self.h_cool_in.v / 1000.0)
            k_water = fp.conductivity("water", P=self.P_cool_in.v / 1000.0, h=self.h_cool_in.v / 1000.0)
            mu_water = mu_water / 1000000.0  # converting from microPa-s to Pa-s
            if mu_water == 0.0:
                mu_water = 0.001  # set as default
            Re_CW = 1000.0 * vel * self.ID.v / mu_water

            # Find heat transfer coefficient on inside of tubes based on cooling water
            if Re_CW > 2300.0:  # turbulent flow -> Gnielinski Correlation
                ff = FricFactor_IC(0.0, Re_CW, ff)
                Pr = mu_water * 4200.0 / k_water
                Nu = ((ff / 8.0) * (Re_CW - 1000.0) * Pr) / (1.0 + 12.7 * (ff / 8.0) ** (0.5) * (Pr ** (2.0 / 3.0) - 1.0))
                h_bar_cool = Nu * k_water / self.ID.v
            else:  # laminar flow -> laminar Nusselt number for uniform temperature on a circular tube
                Nu = 3.67
                h_bar_cool = Nu * k_water / self.ID.v

            if h_bar_cool > 0.0:
                R_cool = 1.0 / (h_bar_cool * A_s_inner)
            else:
                R_cool = 10000000.0

            # Find condensation coefficient on the outside of the tubes
            # -> Incropera and DeWitt Condensation on horizontal pipes
            T_cw = (T_cool_out + T_cool_in) / 2.0
            # TODO-NEEDS UNITS CHECK: FIT_PQ passes P in kPa; densities returned; enthalpy returned in kJ/kg then *1000
            rho_sat_g = fp.density("water", P=P_tank_prev / 1000.0, x=1.0)
            h_sat_g = fp.enthalpy("water", P=P_tank_prev / 1000.0, x=1.0)
            rho_sat_f = fp.density("water", P=P_tank_prev / 1000.0, x=0.0)
            k_L = fp.conductivity("water", P=P_tank_prev / 1000.0, x=0.0)
            h_sat_f = fp.enthalpy("water", P=P_tank_prev / 1000.0, x=0.0)
            # TODO-NEEDS UNITS CHECK: viscosity returned in microPa-s then /1000000 to Pa-s
            mu_L = fp.viscosity("water", P=P_tank_prev / 1000.0, x=0.0)
            mu_L = mu_L / 1000000.0  # converting from microPa-s to Pa-s
            h_sat_g = h_sat_g * 1000.0
            h_sat_f = h_sat_f * 1000.0
            cp_L = f_cp_water(P_tank_prev, T_sat - 2.0)
            # Modified latent heat as recommended by Rosenhow (near EQ.27 of Incropera and DeWitt 2002)
            h_fg_mod = (h_sat_g - h_sat_f) + 0.68 * cp_L * (T_sat - T_cw)
            # Eq. 10.45
            h_bar_cond = 0.729 * (
                (9.81 * rho_sat_f * (rho_sat_f - rho_sat_g) * k_L ** 3.0 * h_fg_mod)
                / (mu_L * (T_sat - T_cw) * OD)
            ) ** (0.25)

            if h_bar_cond > 0.0:
                R_cond = 1.0 / (h_bar_cond * A_s_outer)
            else:
                R_cond = 10000000.0

            R_tot = R_cool + R_cond

            # Calculate the heat transfer from the cooling water
            # TODO-NEEDS UNITS CHECK: FIT_TP passes P in kPa; enthalpy returned in kJ/kg then *1000
            h_cool_out = fp.enthalpy("water", P=self.P_cool_in.v / 1000.0, T=T_cool_out)
            h_cool_out = h_cool_out * 1000.0
            Q_dot_R = (T_sat - T_cw) / (R_tot)  # heat transfer based on the resistance network
            Q_dot_cw = self.m_dot_cool_in.v * (h_cool_out - self.h_cool_in.v)  # heat transfer based on the amount of cooling flow traveling through the type
            Q_dot = min(Q_dot_R, Q_dot_cw)  # limiting heat transfer amount

            # Calculating volumetric cooling flow
            m_dot_cool_out = self.m_dot_cool_in.v
            h_cool_out = (self.m_dot_cool_in.v * self.h_cool_in.v + Q_dot) / self.m_dot_cool_in.v
            # TODO-NEEDS UNITS CHECK: FIT_PH passes P in kPa, h in kJ/kg
            T_cool_out = fp.temperature("water", P=self.P_cool_in.v / 1000.0, h=h_cool_out / 1000.0)
            P_cool_out = self.P_cool_in.v
            Vol_cool_out = m_dot_cool_out / 1000.0
        else:  # No Cooling Flow is entering the Condenser
            m_dot_cool_out = self.m_dot_cool_in.v
            Vol_cool_out = m_dot_cool_out / 1000.0
            h_cool_out = self.h_cool_in.v
            P_cool_out = self.P_cool_in.v
            Q_dot = 0.0
            # TODO-NEEDS UNITS CHECK: FIT_PH passes P in kPa, h in kJ/kg
            T_cool_out = fp.temperature("water", P=self.P_cool_in.v / 1000.0, h=h_cool_out / 1000.0)

        # -------------------Calling OUTPUT Values----------------------------------------------------------------!
        self.m_dot_cool_out.v = m_dot_cool_out  # Cooling mass flow leaving the condenser HX [kg/s]
        self.Vol_cool_out.v = Vol_cool_out       # Volumetric flow rate leaving the condenser HX [m^3/s]
        self.T_cool_out.v = T_cool_out           # Temperature of the Cooling Flow leaving the condenser HX [K]
        self.P_cool_out.v = P_cool_out           # Pressure of the Cooling Flow leaving the condenser HX [Pa]
        self.h_cool_out.v = h_cool_out           # Enthalpy of Cooling Flow leaving the condenser HX [J/kg]
        self.Q_dot_hx.v = Q_dot                  # Total Heat Removed by the condenser HX [W]
        self.ff_out.v = ff                        # friction factor guess
        # --------------------------------------------------------------------------------------------------------!

        # !!!!!!!!!!!!!! STEP 3: CONDENSER HOTWELL TANK EQUATIONS + FINDING RESERVOIR FLOW IN/OUT TO KEEP LEVEL THE SAME!!!!!!!!!!!!!!
        ts = self.model.timestep * 3600.0
        t_crit = 0.1  # critical timestep to solve for the change in pressure in a tank
        Vol_tank = self.Height_tank.v * self.Area_tank.v
        dh = 1000.0
        dP = 1000.0
        m_tank_prev = self.m_tank_cur.v
        L_tank_prev = self.L_tank_cur.v
        rho_f = 1000.0

        # Reservoir conditions (same as tank conditions to not have an effect that will confuse operators)
        m_dot_res = self.m_dot_res_out.v

        if ts <= t_crit:  # Continue
            res_min = -(m_tank_prev - m_dot_pump * ts) + 0.01  # lowest amount of reservoir flow that can be removed based on the mass in the tank and the amount leaving from the pump
            res_max = Vol_tank * rho_f - m_tank_prev
            res_tol = 0.001
            res_error2 = res_tol + 1.0
            res_error_min = 10.0  # placeholder to keep lowest reservoir flow rate that has smallest error if iterations does not converge
            alpha = 0.2  # learning rate for calculating new reservoir flow rate based on previous 2 flows and their errors
            whileiterations = 0.0
            # Variables to track previous iteration for secant method
            m_dot_res1 = m_dot_res
            res_error1 = 0.0

            while abs(res_error2) > res_tol:
                whileiterations = whileiterations + 1.0

                # RK4 Step Forward in Time to solve for new pressure and enthalpy at the end of the timestep
                # AA tank calculations
                # TODO-NEEDS UNITS CHECK: FIT_PQ passes P in kPa; enthalpy returned in kJ/kg then *1000
                T_tank = fp.temperature("water", P=P_tank_prev / 1000.0, x=0.0)
                h_res = fp.enthalpy("water", P=P_tank_prev / 1000.0, x=0.0)
                h_res = h_res * 1000.0  # Converting from KJ/kg to J/kg
                P_pump_in = P_tank_prev + self.L_tank_ini.v * rho_f * 9.81
                # TODO-NEEDS UNITS CHECK: FIT_TP passes P in kPa; enthalpy returned in kJ/kg then *1000
                h_pump = fp.enthalpy("water", P=P_pump_in / 1000.0, T=T_tank)
                h_pump = h_pump * 1000.0  # Converting from KJ/kg to J/kg

                drhodhcp_a = drhodhcp(P_tank=P_tank_prev, h_tank=h_tank_prev, dh=dh)
                drhodpch_a = drhodpch(P_tank=P_tank_prev, h_tank=h_tank_prev, dP=dP)
                dudhcp_a = dudhcp(P_tank=P_tank_prev, h_tank=h_tank_prev, dh=dh)
                dudpch_a = dudpch(P_tank=P_tank_prev, h_tank=h_tank_prev, dP=dP)
                rho_tank = m_tank_prev / Vol_tank
                u_tank = fp.intenergy("water", P=P_tank_prev / 1000.0, h=h_tank_prev / 1000.0)
                u_tank = u_tank * 1000.0  # converting kJ/kg to J/kg

                dPdt_aa = (
                    ((-h_pump + u_tank) * m_dot_pump + (h_res - u_tank) * m_dot_res + (self.h_in.v - u_tank) * self.m_dot_in.v - Q_dot) * drhodhcp_a
                    - rho_tank * dudhcp_a * (self.m_dot_in.v + m_dot_res - m_dot_pump)
                ) / (Vol_tank * rho_tank * (drhodhcp_a * dudpch_a - drhodpch_a * dudhcp_a))
                dhdt_aa = -(Vol_tank * dPdt_aa * drhodpch_a - self.m_dot_in.v + m_dot_pump - m_dot_res) / (Vol_tank * drhodhcp_a)

                P_aa = P_tank_prev + dPdt_aa * ts / 2.0
                h_aa = h_tank_prev + dhdt_aa * ts / 2.0

                # BB tank calculations
                # TODO-NEEDS UNITS CHECK: FIT_PQ passes P in kPa; enthalpy returned in kJ/kg then *1000
                T_tank = fp.temperature("water", P=P_aa / 1000.0, x=0.0)
                h_res = fp.enthalpy("water", P=P_aa / 1000.0, x=0.0)
                h_res = h_res * 1000.0
                P_pump_in = P_aa + self.L_tank_ini.v * rho_f * 9.81
                # TODO-NEEDS UNITS CHECK: FIT_TP passes P in kPa; enthalpy returned in kJ/kg then *1000
                h_pump = fp.enthalpy("water", P=P_pump_in / 1000.0, T=T_tank)
                h_pump = h_pump * 1000.0

                drhodhcp_a = drhodhcp(P_tank=P_aa, h_tank=h_aa, dh=dh)
                drhodpch_a = drhodpch(P_tank=P_aa, h_tank=h_aa, dP=dP)
                dudhcp_a = dudhcp(P_tank=P_aa, h_tank=h_aa, dh=dh)
                dudpch_a = dudpch(P_tank=P_aa, h_tank=h_aa, dP=dP)
                u_tank = fp.intenergy("water", P=P_aa / 1000.0, h=h_aa / 1000.0)
                u_tank = u_tank * 1000.0  # converting kJ/kg to J/kg

                dPdt_bb = (
                    ((-h_pump + u_tank) * m_dot_pump + (h_res - u_tank) * m_dot_res + (self.h_in.v - u_tank) * self.m_dot_in.v - Q_dot) * drhodhcp_a
                    - rho_tank * dudhcp_a * (self.m_dot_in.v + m_dot_res - m_dot_pump)
                ) / (Vol_tank * rho_tank * (drhodhcp_a * dudpch_a - drhodpch_a * dudhcp_a))
                dhdt_bb = -(Vol_tank * dPdt_bb * drhodpch_a - self.m_dot_in.v + m_dot_pump - m_dot_res) / (Vol_tank * drhodhcp_a)

                P_bb = P_tank_prev + dPdt_bb * ts / 2.0
                h_bb = h_tank_prev + dhdt_bb * ts / 2.0

                # CC Tank Calculations
                # TODO-NEEDS UNITS CHECK: FIT_PQ passes P in kPa; enthalpy returned in kJ/kg then *1000
                T_tank = fp.temperature("water", P=P_bb / 1000.0, x=0.0)
                h_res = fp.enthalpy("water", P=P_bb / 1000.0, x=0.0)
                h_res = h_res * 1000.0
                P_pump_in = P_bb + self.L_tank_ini.v + rho_f * 9.81
                # TODO-NEEDS UNITS CHECK: FIT_TP passes P in kPa; enthalpy returned in kJ/kg then *1000
                h_pump = fp.enthalpy("water", P=P_pump_in / 1000.0, T=T_tank)
                h_pump = h_pump * 1000.0

                drhodhcp_a = drhodhcp(P_tank=P_bb, h_tank=h_bb, dh=dh)
                drhodpch_a = drhodpch(P_tank=P_bb, h_tank=h_bb, dP=dP)
                dudhcp_a = dudhcp(P_tank=P_bb, h_tank=h_bb, dh=dh)
                dudpch_a = dudpch(P_tank=P_bb, h_tank=h_bb, dP=dP)
                rho_tank = m_tank_prev / Vol_tank
                u_tank = fp.intenergy("water", P=P_bb / 1000.0, h=h_bb / 1000.0)
                u_tank = u_tank * 1000.0  # converting kJ/kg to J/kg

                dPdt_cc = (
                    ((-h_pump + u_tank) * m_dot_pump + (h_res - u_tank) * m_dot_res + (self.h_in.v - u_tank) * self.m_dot_in.v - Q_dot) * drhodhcp_a
                    - rho_tank * dudhcp_a * (self.m_dot_in.v + m_dot_res - m_dot_pump)
                ) / (Vol_tank * rho_tank * (drhodhcp_a * dudpch_a - drhodpch_a * dudhcp_a))
                dhdt_cc = -(Vol_tank * dPdt_cc * drhodpch_a - self.m_dot_in.v + m_dot_pump - m_dot_res) / (Vol_tank * drhodhcp_a)

                P_cc = P_tank_prev + dPdt_cc * ts
                h_cc = h_tank_prev + dhdt_cc * ts

                # DD Tank Calculations
                # TODO-NEEDS UNITS CHECK: FIT_PQ passes P in kPa; enthalpy returned in kJ/kg then *1000
                T_tank = fp.temperature("water", P=P_cc / 1000.0, x=0.0)
                h_res = fp.enthalpy("water", P=P_cc / 1000.0, x=0.0)
                h_res = h_res * 1000.0
                P_pump_in = P_cc + self.L_tank_ini.v + rho_f * 9.81
                # TODO-NEEDS UNITS CHECK: FIT_TP passes P in kPa; enthalpy returned in kJ/kg then *1000
                h_pump = fp.enthalpy("water", P=P_pump_in / 1000.0, T=T_tank)
                h_pump = h_pump * 1000.0

                drhodhcp_a = drhodhcp(P_tank=P_cc, h_tank=h_cc, dh=dh)
                drhodpch_a = drhodpch(P_tank=P_cc, h_tank=h_cc, dP=dP)
                dudhcp_a = dudhcp(P_tank=P_cc, h_tank=h_cc, dh=dh)
                dudpch_a = dudpch(P_tank=P_cc, h_tank=h_cc, dP=dP)
                # TODO-NEEDS UNITS CHECK: FIT_PH passes P in kPa, h in kJ/kg; internal energy returned in kJ/kg then *1000
                u_tank = fp.intenergy("water", P=P_cc / 1000.0, h=h_cc / 1000.0)
                u_tank = u_tank * 1000.0  # converting kJ/kg to J/kg

                dPdt_dd = (
                    ((-h_pump + u_tank) * m_dot_pump + (h_res - u_tank) * m_dot_res + (self.h_in.v - u_tank) * self.m_dot_in.v - Q_dot) * drhodhcp_a
                    - rho_tank * dudhcp_a * (self.m_dot_in.v + m_dot_res - m_dot_pump)
                ) / (Vol_tank * rho_tank * (drhodhcp_a * dudpch_a - drhodpch_a * dudhcp_a))
                dhdt_dd = -(Vol_tank * dPdt_dd * drhodpch_a - self.m_dot_in.v + m_dot_pump - m_dot_res) / (Vol_tank * drhodhcp_a)

                # Solving for pressure at the end of the timestep
                P_tank_new = P_tank_prev + (dPdt_aa + 2.0 * dPdt_bb + 2.0 * dPdt_cc + dPdt_dd) * ts / 6.0
                h_tank_new = h_tank_prev + (dhdt_aa + 2.0 * dhdt_bb + 2.0 * dhdt_cc + dhdt_dd) * ts / 6.0
                # TODO-NEEDS UNITS CHECK: FIT_PH passes P in kPa, h in kJ/kg
                x_tank_new = fp.quality("water", P=P_tank_new / 1000.0, h=h_tank_new / 1000.0)
                rho_tank_new = fp.density("water", P=P_tank_new / 1000.0, h=h_tank_new / 1000.0)
                T_tank_new = fp.temperature("water", P=P_tank_new / 1000.0, h=h_tank_new / 1000.0)

                dmdt = self.m_dot_in.v - m_dot_pump + m_dot_res
                m_tank_new = m_tank_prev + dmdt * ts
                m_tank_g_new = max(m_tank_new * x_tank_new, 0.0)
                m_tank_f_new = max(m_tank_new - m_tank_g_new, 0.0)

                L_tank_new = m_tank_f_new / self.Area_tank.v / (rho_f)

                # Keep level constant by adjusting flow coming/going from reservoir
                res_error2 = L_tank_new - self.L_tank_ini.v

                if abs(res_error2) < abs(res_error_min):
                    res_error_min = res_error2
                    m_dot_res_min = m_dot_res

                m_dot_res2 = m_dot_res  # save m_dot_res as new second point
                # Calculating new reservoir mass to get Level change = 0

                if abs(res_error2) < res_tol:
                    break
                elif whileiterations == 100.0:
                    m_dot_res = m_dot_res_min
                elif whileiterations > 100.0:
                    break
                elif whileiterations > 1.0:
                    if abs(m_dot_res2 - m_dot_res1) < 0.0001:  # computing slope between points will be undefined
                        if res_error2 > 0.0:  # New tank level is higher than desired, decrease m_dot_res
                            m_dot_res3 = max(m_dot_res2 - 10.0, res_min)
                        else:               # new tank level is lower than desired, increase m_dot_res
                            m_dot_res3 = min(m_dot_res2 + 10.0, res_max)
                    else:
                        m_res = (res_error2 - res_error1) / (m_dot_res2 - m_dot_res1)
                        b_res = res_error2 - m_res * m_dot_res2
                        if m_res != 0.0:
                            m_dot_res3p = max(min(-b_res / m_res, res_max), res_min)  # calculate m_dot_res 3 prime which is where it wants to go next
                            delta_res = m_dot_res3p - m_dot_res2
                            if res_error2 > 0.0:  # need to remove water from tank, m_dot_res should decrease
                                if delta_res < 0.0:  # all good
                                    m_dot_res3 = m_dot_res2 + alpha * (delta_res)
                                else:
                                    m_dot_res3 = max(m_dot_res2 - 10.0, res_min)
                            else:
                                if delta_res > 0.0:  # all good continue
                                    m_dot_res3 = m_dot_res2 + alpha * (delta_res)
                                else:
                                    m_dot_res3 = min(m_dot_res2 + 10.0, res_max)
                        else:
                            if res_error2 > 0.0:  # New tank level is higher than desired, decrease m_dot_res
                                m_dot_res3 = max(m_dot_res2 - 10.0, res_min)
                            else:               # new tank level is lower than desired, increase m_dot_res
                                m_dot_res3 = min(m_dot_res2 + 10.0, res_max)
                else:
                    if res_error2 > 0.0:  # New tank level is higher than desired, decrease m_dot_res
                        m_dot_res3 = max(m_dot_res2 - 10.0, res_min)
                    else:               # new tank level is lower than desired, increase m_dot_res
                        m_dot_res3 = min(m_dot_res2 + 10.0, res_max)

                m_dot_res1 = m_dot_res2
                res_error1 = res_error2

                m_dot_res = m_dot_res3

        else:  # timestep is larger than critical timestep for condenser, break up timesteps into subtimesteps
            ts_sub_f = ts / t_crit
            ts_sub_i = math.ceil(ts_sub_f)
            ts_sub_f = ts / float(ts_sub_i)
            res_tol = 0.001
            res_error2 = res_tol + 1.0
            res_error_min = 10.0  # placeholder to keep lowest reservoir flow rate that has smallest error if iterations does not converge
            alpha = 0.01  # learning rate for calculating new reservoir flow rate based on previous 2 flows and their errors
            # Initialize variables for sub-step loop
            P_tank_new = P_tank_prev
            h_tank_new = h_tank_prev
            m_tank_new = m_tank_prev
            T_tank_new = fp.temperature("water", P=P_tank_prev / 1000.0, x=0.0)  # TODO-NEEDS UNITS CHECK

            for i in range(1, ts_sub_i + 1):
                if i > 1:
                    P_tank_prev = P_tank_new
                    h_tank_prev = h_tank_new
                    m_tank_prev = m_tank_new
                res_min = -(m_tank_prev - m_dot_pump * ts_sub_f) + 0.01  # lowest amount of reservoir flow that can be removed based on the mass in the tank and the amount leaving from the pump
                res_max = Vol_tank * rho_f - m_tank_prev
                whileiterations = 0.0
                m_dot_res1 = m_dot_res
                res_error1 = 0.0

                while abs(res_error2) > res_tol:
                    whileiterations = whileiterations + 1.0

                    # RK4 Step Forward in Time to solve for new pressure and enthalpy at the end of the timestep
                    # AA tank calculations
                    # TODO-NEEDS UNITS CHECK: FIT_PQ passes P in kPa; enthalpy returned in kJ/kg then *1000
                    T_tank = fp.temperature("water", P=P_tank_prev / 1000.0, x=0.0)
                    h_res = fp.enthalpy("water", P=P_tank_prev / 1000.0, x=0.0)
                    h_res = h_res * 1000.0  # Converting from KJ/kg to J/kg
                    P_pump_in = P_tank_prev + self.L_tank_ini.v * rho_f * 9.81
                    # TODO-NEEDS UNITS CHECK: FIT_TP passes P in kPa; enthalpy returned in kJ/kg then *1000
                    h_pump = fp.enthalpy("water", P=P_pump_in / 1000.0, T=T_tank)
                    h_pump = h_pump * 1000.0  # Converting from KJ/kg to J/kg

                    drhodhcp_a = drhodhcp(P_tank=P_tank_prev, h_tank=h_tank_prev, dh=dh)
                    drhodpch_a = drhodpch(P_tank=P_tank_prev, h_tank=h_tank_prev, dP=dP)
                    dudhcp_a = dudhcp(P_tank=P_tank_prev, h_tank=h_tank_prev, dh=dh)
                    dudpch_a = dudpch(P_tank=P_tank_prev, h_tank=h_tank_prev, dP=dP)
                    rho_tank = m_tank_prev / Vol_tank
                    u_tank = fp.intenergy("water", P=P_tank_prev / 1000.0, h=h_tank_prev / 1000.0)
                    u_tank = u_tank * 1000.0  # converting kJ/kg to J/kg

                    dPdt_aa = (
                        ((-h_pump + u_tank) * m_dot_pump + (h_res - u_tank) * m_dot_res + (self.h_in.v - u_tank) * self.m_dot_in.v - Q_dot) * drhodhcp_a
                        - rho_tank * dudhcp_a * (self.m_dot_in.v + m_dot_res - m_dot_pump)
                    ) / (Vol_tank * rho_tank * (drhodhcp_a * dudpch_a - drhodpch_a * dudhcp_a))
                    dhdt_aa = -(Vol_tank * dPdt_aa * drhodpch_a - self.m_dot_in.v + m_dot_pump - m_dot_res) / (Vol_tank * drhodhcp_a)

                    P_aa = P_tank_prev + dPdt_aa * ts_sub_f / 2.0
                    h_aa = h_tank_prev + dhdt_aa * ts_sub_f / 2.0

                    # BB tank calculations
                    # TODO-NEEDS UNITS CHECK: FIT_PQ passes P in kPa; enthalpy returned in kJ/kg then *1000
                    T_tank = fp.temperature("water", P=P_aa / 1000.0, x=0.0)
                    h_res = fp.enthalpy("water", P=P_aa / 1000.0, x=0.0)
                    h_res = h_res * 1000.0
                    P_pump_in = P_aa + self.L_tank_ini.v * rho_f * 9.81
                    # TODO-NEEDS UNITS CHECK: FIT_TP passes P in kPa; enthalpy returned in kJ/kg then *1000
                    h_pump = fp.enthalpy("water", P=P_pump_in / 1000.0, T=T_tank)
                    h_pump = h_pump * 1000.0

                    drhodhcp_a = drhodhcp(P_tank=P_aa, h_tank=h_aa, dh=dh)
                    drhodpch_a = drhodpch(P_tank=P_aa, h_tank=h_aa, dP=dP)
                    dudhcp_a = dudhcp(P_tank=P_aa, h_tank=h_aa, dh=dh)
                    dudpch_a = dudpch(P_tank=P_aa, h_tank=h_aa, dP=dP)
                    u_tank = fp.intenergy("water", P=P_aa / 1000.0, h=h_aa / 1000.0)
                    u_tank = u_tank * 1000.0  # converting kJ/kg to J/kg

                    dPdt_bb = (
                        ((-h_pump + u_tank) * m_dot_pump + (h_res - u_tank) * m_dot_res + (self.h_in.v - u_tank) * self.m_dot_in.v - Q_dot) * drhodhcp_a
                        - rho_tank * dudhcp_a * (self.m_dot_in.v + m_dot_res - m_dot_pump)
                    ) / (Vol_tank * rho_tank * (drhodhcp_a * dudpch_a - drhodpch_a * dudhcp_a))
                    dhdt_bb = -(Vol_tank * dPdt_bb * drhodpch_a - self.m_dot_in.v + m_dot_pump - m_dot_res) / (Vol_tank * drhodhcp_a)

                    P_bb = P_tank_prev + dPdt_bb * ts_sub_f / 2.0
                    h_bb = h_tank_prev + dhdt_bb * ts_sub_f / 2.0

                    # CC Tank Calculations
                    # TODO-NEEDS UNITS CHECK: FIT_PQ passes P in kPa; enthalpy returned in kJ/kg then *1000
                    T_tank = fp.temperature("water", P=P_bb / 1000.0, x=0.0)
                    h_res = fp.enthalpy("water", P=P_bb / 1000.0, x=0.0)
                    h_res = h_res * 1000.0
                    P_pump_in = P_bb + self.L_tank_ini.v + rho_f * 9.81
                    # TODO-NEEDS UNITS CHECK: FIT_TP passes P in kPa; enthalpy returned in kJ/kg then *1000
                    h_pump = fp.enthalpy("water", P=P_pump_in / 1000.0, T=T_tank)
                    h_pump = h_pump * 1000.0

                    drhodhcp_a = drhodhcp(P_tank=P_bb, h_tank=h_bb, dh=dh)
                    drhodpch_a = drhodpch(P_tank=P_bb, h_tank=h_bb, dP=dP)
                    dudhcp_a = dudhcp(P_tank=P_bb, h_tank=h_bb, dh=dh)
                    dudpch_a = dudpch(P_tank=P_bb, h_tank=h_bb, dP=dP)
                    rho_tank = m_tank_prev / Vol_tank
                    u_tank = fp.intenergy("water", P=P_bb / 1000.0, h=h_bb / 1000.0)
                    u_tank = u_tank * 1000.0  # converting kJ/kg to J/kg

                    dPdt_cc = (
                        ((-h_pump + u_tank) * m_dot_pump + (h_res - u_tank) * m_dot_res + (self.h_in.v - u_tank) * self.m_dot_in.v - Q_dot) * drhodhcp_a
                        - rho_tank * dudhcp_a * (self.m_dot_in.v + m_dot_res - m_dot_pump)
                    ) / (Vol_tank * rho_tank * (drhodhcp_a * dudpch_a - drhodpch_a * dudhcp_a))
                    dhdt_cc = -(Vol_tank * dPdt_cc * drhodpch_a - self.m_dot_in.v + m_dot_pump - m_dot_res) / (Vol_tank * drhodhcp_a)

                    P_cc = P_tank_prev + dPdt_cc * ts_sub_f
                    h_cc = h_tank_prev + dhdt_cc * ts_sub_f

                    # DD Tank Calculations
                    # TODO-NEEDS UNITS CHECK: FIT_PQ passes P in kPa; density returned; enthalpy in kJ/kg then *1000
                    T_tank = fp.temperature("water", P=P_cc / 1000.0, x=0.0)
                    h_res = fp.enthalpy("water", P=P_cc / 1000.0, x=0.0)
                    rho_f = fp.density("water", P=P_cc / 1000.0, x=0.0)
                    h_res = h_res * 1000.0
                    P_pump_in = P_cc + self.L_tank_ini.v + rho_f * 9.81
                    # TODO-NEEDS UNITS CHECK: FIT_TP passes P in kPa; enthalpy returned in kJ/kg then *1000
                    h_pump = fp.enthalpy("water", P=P_pump_in / 1000.0, T=T_tank)
                    h_pump = h_pump * 1000.0

                    drhodhcp_a = drhodhcp(P_tank=P_cc, h_tank=h_cc, dh=dh)
                    drhodpch_a = drhodpch(P_tank=P_cc, h_tank=h_cc, dP=dP)
                    dudhcp_a = dudhcp(P_tank=P_cc, h_tank=h_cc, dh=dh)
                    dudpch_a = dudpch(P_tank=P_cc, h_tank=h_cc, dP=dP)
                    u_tank = fp.intenergy("water", P=P_cc / 1000.0, h=h_cc / 1000.0)
                    u_tank = u_tank * 1000.0  # converting kJ/kg to J/kg

                    dPdt_dd = (
                        ((-h_pump + u_tank) * m_dot_pump + (h_res - u_tank) * m_dot_res + (self.h_in.v - u_tank) * self.m_dot_in.v - Q_dot) * drhodhcp_a
                        - rho_tank * dudhcp_a * (self.m_dot_in.v + m_dot_res - m_dot_pump)
                    ) / (Vol_tank * rho_tank * (drhodhcp_a * dudpch_a - drhodpch_a * dudhcp_a))
                    dhdt_dd = -(Vol_tank * dPdt_dd * drhodpch_a - self.m_dot_in.v + m_dot_pump - m_dot_res) / (Vol_tank * drhodhcp_a)

                    # Solving for pressure at the end of the timestep
                    P_tank_new = P_tank_prev + (dPdt_aa + 2.0 * dPdt_bb + 2.0 * dPdt_cc + dPdt_dd) * ts_sub_f / 6.0
                    h_tank_new = h_tank_prev + (dhdt_aa + 2.0 * dhdt_bb + 2.0 * dhdt_cc + dhdt_dd) * ts_sub_f / 6.0
                    # TODO-NEEDS UNITS CHECK: FIT_PH passes P in kPa, h in kJ/kg
                    x_tank_new = fp.quality("water", P=P_tank_new / 1000.0, h=h_tank_new / 1000.0)
                    rho_tank_new = fp.density("water", P=P_tank_new / 1000.0, h=h_tank_new / 1000.0)
                    T_tank_new = fp.temperature("water", P=P_tank_new / 1000.0, h=h_tank_new / 1000.0)

                    dmdt = self.m_dot_in.v - m_dot_pump + m_dot_res
                    m_tank_new = m_tank_prev + dmdt * ts_sub_f
                    m_tank_g_new = max(m_tank_new * x_tank_new, 0.0)
                    m_tank_f_new = max(m_tank_new - m_tank_g_new, 0.0)

                    L_tank_new = m_tank_f_new / self.Area_tank.v / rho_f

                    # Keep level constant by adjusting flow coming/going from reservoir
                    res_error2 = L_tank_new - self.L_tank_ini.v

                    if abs(res_error2) < abs(res_error_min):
                        res_error_min = res_error2
                        m_dot_res_min = m_dot_res

                    m_dot_res2 = m_dot_res  # save m_dot_res as new second point
                    # Calculating new reservoir mass to get Level change = 0

                    if abs(res_error2) < res_tol:
                        break
                    elif whileiterations == 200.0:
                        m_dot_res3 = m_dot_res_min
                    elif whileiterations > 200.0:
                        break
                    elif whileiterations > 1.0:
                        if abs(m_dot_res2 - m_dot_res1) < 0.0001:  # computing slope between points will be undefined
                            if res_error2 > 0.0:  # New tank level is higher than desired, decrease m_dot_res
                                m_dot_res3 = max(m_dot_res2 - 10.0, res_min)
                            else:               # new tank level is lower than desired, increase m_dot_res
                                m_dot_res3 = min(m_dot_res2 + 10.0, res_max)
                        else:
                            m_res = (res_error2 - res_error1) / (m_dot_res2 - m_dot_res1)
                            b_res = res_error2 - m_res * m_dot_res2
                            if m_res != 0.0:
                                m_dot_res3p = max(min(-b_res / m_res, res_max), res_min)  # calculate m_dot_res 3 prime which is where it wants to go next
                                delta_res = m_dot_res3p - m_dot_res2
                                if res_error2 > 0.0:  # need to remove water from tank, m_dot_res should decrease
                                    if delta_res < 0.0:  # all good
                                        m_dot_res3 = m_dot_res2 + alpha * (delta_res)
                                    else:
                                        m_dot_res3 = max(m_dot_res2 - 10.0, res_min)
                                else:
                                    if delta_res > 0.0:  # all good continue
                                        m_dot_res3 = m_dot_res2 + alpha * (delta_res)
                                    else:
                                        m_dot_res3 = min(m_dot_res2 + 10.0, res_max)
                            else:
                                if res_error2 > 0.0:  # New tank level is higher than desired, decrease m_dot_res
                                    m_dot_res3 = max(m_dot_res2 - 10.0, res_min)
                                else:               # new tank level is lower than desired, increase m_dot_res
                                    m_dot_res3 = min(m_dot_res2 + 10.0, res_max)
                    else:
                        if res_error2 > 0.0:  # New tank level is higher than desired, decrease m_dot_res
                            m_dot_res3 = max(m_dot_res2 - 10.0, res_min)
                        else:               # new tank level is lower than desired, increase m_dot_res
                            m_dot_res3 = min(m_dot_res2 + 10.0, res_max)

                    m_dot_res1 = m_dot_res2
                    res_error1 = res_error2

                    m_dot_res = m_dot_res3

        # ---------------------CALLING OUTPUT VALUES------------------------------------------------------------------------!
        self.L_tank_next.v = L_tank_new     # Level of the tank for the next timestep [m]
        self.P_tank_next.v = P_tank_new     # Pressure of the tank for the next timestep [Pa]
        self.h_tank_next.v = h_tank_new     # Enthalpy of the tank for the next timestep [J/kg]
        self.T_tank_next.v = T_tank_new     # Temperature of the tank for the next timestep [K]
        self.m_tank_next.v = m_tank_new     # Mass in the tank for the next timestep [kg]
        self.m_dot_res_out.v = m_dot_res    # Mass need from reservoir to make level change 0
