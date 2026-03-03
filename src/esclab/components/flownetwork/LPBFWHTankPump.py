"""LP Boiler Feedwater Heater Tank-Pump component model (Type 6014)."""

import eeslib.fluid_properties as fp

from esclab.simulate import Component
from esclab.components.flownetwork.esol6015_helpers import PB_CV_data


class LPBFWHTankPump(Component):
    """
    Object: LPBFWH Tank-Pump
    Simulation Studio Model: ESOL6014-LPBFWH-Tank-Pump

    Author:
    Editor:
    Date:     June 01, 2023
    last modified: June 01, 2023
    """

    # *** Model Parameters ***
    #   H_tank        - Height of the vertical cylindrical receiver tank
    #   D_tank        - Diameter of the vertical cylindrical receiver tank
    #   perc_tank_ini - Initial percentage that the tank is filled with water
    #   T_tank_ini    - Initial temperature of the tank
    #   Valve_speed   - Maximum speed [deg/s] that the valve position can change at
    #   D_valve       - Diameter of the valve
    #   Valve_type    - Type of valve used (1 = concentric butterfly valve, 2 = triple offset butterfly valve)
    #   PC_Coef_A     - Coefficient A used to determine pump head - Coef_A * Flow^2 + Coef_B * Flow + Coef_C = Pump_Head
    #   PC_Coef_B     - Coefficient B used to determine pump head - Coef_A * Flow^2 + Coef_B * Flow + Coef_C = Pump_Head
    #   PC_Coef_C     - Coefficient C used to determine pump head - Coef_A * Flow^2 + Coef_B * Flow + Coef_C = Pump_Head
    #   Eta_Coef_A    - Coefficient A used to determine pump efficiency - Coef_A*Flow^4 + Coef_B*Flow^3 + Coef_C*Flow^2 + Coef_D*Flow = Eta_pump
    #   Eta_Coef_B    - Coefficient B used to determine pump efficiency - Coef_A*Flow^4 + Coef_B*Flow^3 + Coef_C*Flow^2 + Coef_D*Flow = Eta_pump
    #   Eta_Coef_C    - Coefficient C used to determine pump efficiency - Coef_A*Flow^4 + Coef_B*Flow^3 + Coef_C*Flow^2 + Coef_D*Flow = Eta_pump
    #   Eta_Coef_D    - Coefficient D used to determine pump efficiency - Coef_A*Flow^4 + Coef_B*Flow^3 + Coef_C*Flow^2 + Coef_D*Flow = Eta_pump
    H_tank = Component.Parameter()         # Height of the vertical cylindrical receiver tank
    D_tank = Component.Parameter()         # Diameter of the vertical cylindrical receiver tank
    perc_tank_ini = Component.Parameter()  # Initial percentage that the tank is filled with water
    T_tank_ini = Component.Parameter()     # Initial temperature of the tank
    Valve_speed = Component.Parameter()    # Maximum speed [deg/s] that the valve position can change at
    D_valve = Component.Parameter()        # Diameter of the valve
    Valve_type = Component.Parameter()     # Type of valve used (1 = concentric butterfly valve, 2 = triple offset butterfly valve)
    PC_Coef_A = Component.Parameter()     # Coefficient A used to determine pump head - Coef_A * Flow^2 + Coef_B * Flow + Coef_C = Pump_Head
    PC_Coef_B = Component.Parameter()     # Coefficient B used to determine pump head - Coef_A * Flow^2 + Coef_B * Flow + Coef_C = Pump_Head
    PC_Coef_C = Component.Parameter()     # Coefficient C used to determine pump head - Coef_A * Flow^2 + Coef_B * Flow + Coef_C = Pump_Head
    Eta_Coef_A = Component.Parameter()    # Coefficient A used to determine pump efficiency - Coef_A*Flow^4 + Coef_B*Flow^3 + Coef_C*Flow^2 + Coef_D*Flow = Eta_pump
    Eta_Coef_B = Component.Parameter()    # Coefficient B used to determine pump efficiency - Coef_A*Flow^4 + Coef_B*Flow^3 + Coef_C*Flow^2 + Coef_D*Flow = Eta_pump
    Eta_Coef_C = Component.Parameter()    # Coefficient C used to determine pump efficiency - Coef_A*Flow^4 + Coef_B*Flow^3 + Coef_C*Flow^2 + Coef_D*Flow = Eta_pump
    Eta_Coef_D = Component.Parameter()    # Coefficient D used to determine pump efficiency - Coef_A*Flow^4 + Coef_B*Flow^3 + Coef_C*Flow^2 + Coef_D*Flow = Eta_pump

    # *** Model Inputs ***
    #   auto_control  - manual or automatic (1 if automatic otherwise manual mode)
    #   pump_power    - Pump power state (1 if on, 0 if off)
    #   VP_input      - Valve position input signal must be between 0 and 1
    #   Turbine_ON    - Turbine on signal (1 indicated turbine is on, anything else means turbine is off)
    #   m_dot_fw_in   - mass flow rate of feedwater that bfwh extraction will be added to
    #   h_fw_in       - enthalpy of feedwater that bfwh extraction will mix with
    #   P_fw_in       - Pressure of feedwater that bfwh extraction will be added to
    #   m_dot_BFWH   - mass flow rate of water entering receiver tank from Boiler feedwater heaters
    #   h_BFWH        - enthalpy of water entering receiver tank from boiler feedwater heater
    #   P_BFWH        - Pressure of water entering receiver tank from boiler feedwater heater
    #   PID_signal    - Valve Position Signal Determined by PID Controllers
    auto_control = Component.Input()   # manual or automatic (1 if automatic otherwise manual mode)
    pump_power = Component.Input()     # Pump power state (1 if on, 0 if off)
    VP_input = Component.Input()       # Valve position input signal must be between 0 and 1
    Turbine_ON = Component.Input()     # Turbine on signal (1 indicated turbine is on, anything else means turbine is off)
    m_dot_fw_in = Component.Input()    # mass flow rate of feedwater that bfwh extraction will be added to
    h_fw_in = Component.Input()        # enthalpy of feedwater that bfwh extraction will mix with
    P_fw_in = Component.Input()        # Pressure of feedwater that bfwh extraction will be added to
    m_dot_BFWH = Component.Input()    # mass flow rate of water entering receiver tank from Boiler feedwater heaters
    h_BFWH = Component.Input()         # enthalpy of water entering receiver tank from boiler feedwater heater
    P_BFWH = Component.Input()         # Pressure of water entering receiver tank from boiler feedwater heater
    PID_signal = Component.Input()     # Valve Position Signal Determined by PID Controllers

    # *** Model Outputs ***
    #   VP_output      - Valve position output [0-1]
    #   m_dot_fw_out   - Feedwater exiting through main tee [kg/s]
    #   h_fw_out       - Enthalpy of feedwater exiting [J/kg]
    #   P_fw_out       - Pressure of feedwater exiting [Pa]
    #   m_dot_pump     - Mass flow rate through pump [kg/s]
    #   P_pump_out     - Pressure at pump outlet [Pa]
    #   W_dot_pump     - Pump power [W]
    #   L_tank_start   - Level of the tank at the start of a timestep [m]  (start-of-ts state, updated at convergence)
    #   L_tank_end     - Level of the tank at the end of a timestep [m]
    #   T_tank_start   - Temperature of the tank at the start of a timestep [K]  (start-of-ts state, updated at convergence)
    #   T_tank_end     - Temperature of the tank at the end of a timestep [K]
    #   m_tank_start   - Mass in the tank at the start of the timestep [kg]  (start-of-ts state, updated at convergence)
    #   m_tank_end     - Mass in the tank at the end of the timestep [kg]
    #   h_tank_start   - Enthalpy of the tank at the start of the timestep [J/kg]  (start-of-ts state, updated at convergence)
    #   h_tank_end     - Enthalpy of the tank at the end of the timestep [J/kg]
    VP_output = Component.Output()      # output 1
    m_dot_fw_out = Component.Output()   # output 2
    h_fw_out = Component.Output()       # output 3
    P_fw_out = Component.Output()       # output 4
    m_dot_pump = Component.Output()     # output 5
    P_pump_out = Component.Output()     # output 6
    W_dot_pump = Component.Output()     # output 7
    L_tank_start = Component.Output()   # output 8  - level at start of timestep (updated at convergence)
    L_tank_end = Component.Output()     # output 9  - level at end of timestep
    T_tank_start = Component.Output()   # output 10 - temperature at start of timestep (updated at convergence)
    T_tank_end = Component.Output()     # output 11 - temperature at end of timestep
    m_tank_start = Component.Output()   # output 12 - mass at start of timestep (updated at convergence)
    m_tank_end = Component.Output()     # output 13 - mass at end of timestep
    h_tank_start = Component.Output()   # output 14 - enthalpy at start of timestep (updated at convergence)
    h_tank_end = Component.Output()     # output 15 - enthalpy at end of timestep

    # -----------------------------------------------------------------------------------------------------------------------
    def presim_setup(self, **kwargs):
        """Initialize outputs for first timestep."""
        P_tank = 101325.0
        # Make sure initial tank temperature is not saturated
        T_tank_ini_val = min(self.T_tank_ini.v, 373.15)
        # Determine initial tank enthalpy
        # CONVERTED-NEEDS UNITS CHECK: removed unit scale factors; eeslib uses SI
        h_tank_start_val = fp.enthalpy("water", T=T_tank_ini_val, P=P_tank)

        # Determine the amount of mass in the tank at the start of the timestep
        rho_water = 1000.0  # density of water
        L_tank_start_val = self.perc_tank_ini.v * self.H_tank.v
        m_tank_start_val = L_tank_start_val * 3.14 * (self.D_tank.v / 2.0) ** 2 * rho_water

        self.VP_output.v = self.VP_input.v
        self.L_tank_start.v = L_tank_start_val    # Level of the tank at the start of a timestep
        self.L_tank_end.v = L_tank_start_val
        self.T_tank_start.v = T_tank_ini_val      # Temperature of the tank at the start of a timestep
        self.T_tank_end.v = T_tank_ini_val
        self.m_tank_start.v = m_tank_start_val    # Mass in the tank at the start of the timestep
        self.m_tank_end.v = m_tank_start_val
        self.h_tank_start.v = h_tank_start_val    # Enthalpy of the tank at the start of the timestep
        self.h_tank_end.v = h_tank_start_val

    # -----------------------------------------------------------------------------------------------------------------------
    def calculate(self):

        # -----------------------------------------------------------------------------------------------------------------------
        # Perform Any "After Convergence" Manipulations That May Be Required at the End of Each Timestep
        if self.model.is_converged:
            # UPDATE TANK CONDITIONS
            # Set tank level at the end of this timestep as the tank level at the start of the next timestep
            self.L_tank_start.v = self.L_tank_end.v
            # Set tank temperature at the end of this timestep as the tank temperature at the start of the next timestep
            self.T_tank_start.v = self.T_tank_end.v
            # Set tank mass at the end of this timestep as the tank mass at the start of the next timestep
            self.m_tank_start.v = self.m_tank_end.v
            # Set tank enthalpy at the end of this timestep as the tank enthalpy at the start of the next timestep
            self.h_tank_start.v = self.h_tank_end.v
            return

        P_tank = 101325.0
        rho_water = 1000.0
        ts = self.model.timestep * 3600.0  # Converting timestep from hr to s

        # max pump head and max pump flow based on pump curve
        tol = 1000.0  # tolerance for while loops
        LR = 0.4  # Learning rate for while loops
        Coef_C_adj = max(self.PC_Coef_C.v - (self.P_fw_in.v - P_tank) / rho_water / 9.81, 1.0)
        m_dot_pump_max = (-self.PC_Coef_B.v - (self.PC_Coef_B.v ** 2.0 - 4 * self.PC_Coef_A.v * Coef_C_adj) ** 0.5) / (2 * self.PC_Coef_A.v) * rho_water
        m_dot_pump_min = 0.00000001

        ##### Update valve position based on PID controller or manual control if First Iteration in Timestep #####
        if self.auto_control.v == 1.0:  # Automatic Mode
            # Ensure pump power is on
            self.pump_power.v = 1.0
            if self.model.iteration == 0:  # First Iteration in Timestep, adjust VP Output
                if self.model.time != self.model.timestep:
                    # Valve Position Requested, converting from percent open to degrees
                    VP_input_d = self.VP_input.v * 90.0
                    # Last Timesteps Valve Position, converting from percent open to degrees
                    VP_output_d = self.VP_output.v * 90.0
                    # Convert timestep from hr to s
                    ts_val = self.model.timestep * 3600.0
                    if VP_input_d == VP_output_d:  # If input requested matched last timesteps input then do nothing
                        self.VP_output.v = self.VP_input.v
                    elif VP_input_d > VP_output_d:  # If input requested is greater than last timesteps value, open valve position more
                        self.VP_output.v = min(VP_output_d + (self.Valve_speed.v * ts_val), VP_input_d) / 90.0
                    else:  # If input requested is lower than last timesteps value, close valve more
                        self.VP_output.v = max(VP_output_d - (self.Valve_speed.v * ts_val), VP_input_d) / 90.0
                else:  # First Timestep do not change valve position
                    self.VP_output.v = self.PID_signal.v
                # Make sure valve input is not at zero because it will cause divide by zero errors
                if self.VP_output.v <= 0.00001:
                    self.VP_output.v = 0.00001
        else:  # Manual Control
            if self.model.iteration == 0:  # First Iteration in Timestep, adjust VP Output
                if self.model.time != self.model.timestep:
                    # Valve Position Requested, converting from percent open to degrees
                    VP_input_d = self.VP_input.v * 90.0
                    # Last Timesteps Valve Position, converting from percent open to degrees
                    VP_output_d = self.VP_output.v * 90.0
                    # Convert timestep from hr to s
                    ts_val = self.model.timestep * 3600.0
                    if VP_input_d == VP_output_d:  # If input requested matched last timesteps input then do nothing
                        self.VP_output.v = self.VP_input.v
                    elif VP_input_d > VP_output_d:  # If input requested is greater than last timesteps value, open valve position more
                        self.VP_output.v = min(VP_output_d + (self.Valve_speed.v * ts_val), VP_input_d) / 90.0
                    else:  # If input requested is lower than last timesteps value, close valve more
                        self.VP_output.v = max(VP_output_d - (self.Valve_speed.v * ts_val), VP_input_d) / 90.0
                else:  # First Timestep do not change valve position
                    self.VP_output.v = self.VP_input.v
                # Make sure valve input is not at zero because it will cause divide by zero errors
                if self.VP_output.v <= 0.00001:
                    self.VP_output.v = 0.00001

        if self.pump_power.v != 1.0:  # Pump is OFF
            m_dot_pump = 0.0
            P_pump_out = self.P_fw_in.v
            h_pump_out = self.h_fw_in.v
            W_dot_pump = 0.0
        else:  # Pump is ON
            # Find the amount of flow leaving the pump based on the valve position
            CV = PB_CV_data(int(self.Valve_type.v), self.D_valve.v, self.VP_output.v)  # retrieve CV value based on current valve position
            m_dot_pump = self.m_dot_pump.v
            if m_dot_pump == 0.0:
                m_dot_pump = m_dot_pump_min
            error = tol + 1.0
            whileiterations = 0.0
            m_dot_pump_prev = m_dot_pump
            error_prev = 0.0
            do_iterate = True
            while abs(error) > tol:
                whileiterations = whileiterations + 1.0
                Vol_in = m_dot_pump / rho_water                     # Volumetric flow entering the valve [m^3/s]
                Vol_in_gpm = Vol_in * 15850.3                       # Volumetric flow entering valve [GPM]
                DELTA_P_psi = 1.0 / (CV / Vol_in_gpm) ** 2.0       # pressure drop in psi across valve
                DELTA_P = DELTA_P_psi * 6894.76                     # pressure drop in Pa across valve
                P_pump_guess = self.P_fw_in.v + DELTA_P - P_tank - self.L_tank_start.v * 1000.0 * 9.81  # Necessary pump head to overcome valve losses
                P_pump_out = (self.PC_Coef_A.v * Vol_in ** 2.0 + self.PC_Coef_B.v * Vol_in + self.PC_Coef_C.v) * 9.81 * rho_water + P_tank
                error = P_pump_guess - P_pump_out
                if m_dot_pump < 0.001:
                    if error > 0.0:
                        m_dot_pump = m_dot_pump_min
                        error = tol / 2.0
                if abs(m_dot_pump - m_dot_pump_max) < 0.001:
                    if error < 0.0:
                        m_dot_pump = m_dot_pump_max
                        error = tol / 2.0

                if whileiterations == 1.0:  # only one point of information
                    if abs(error) > tol:
                        m_dot_pump_prev = m_dot_pump
                        error_prev = error
                        if error < 0.0:
                            m_dot_pump = min(m_dot_pump + 1.0, m_dot_pump_max)
                        else:
                            m_dot_pump = max(m_dot_pump - 1.0, m_dot_pump_min)
                else:  # two points of information
                    if abs(error) > tol:
                        if m_dot_pump_prev != m_dot_pump:
                            m = (error_prev - error) / (m_dot_pump_prev - m_dot_pump)
                            y_int = error_prev - m * m_dot_pump_prev
                            if m != 0.0:
                                m_dot_new = -y_int / m
                                m_dot_new = min(max(m_dot_new, m_dot_pump_min), m_dot_pump_max)
                                m_dot_pump_prev = m_dot_pump
                                error_prev = error
                                m_dot_pump = m_dot_pump + (m_dot_new - m_dot_pump) * LR
                            else:
                                m_dot_pump_prev = m_dot_pump
                                error_prev = error
                                if error < 0.0:
                                    m_dot_pump = min(m_dot_pump + 1.0, m_dot_pump_max)
                                else:
                                    m_dot_pump = max(m_dot_pump - 1.0, m_dot_pump_min)
                        else:
                            m_dot_pump_prev = m_dot_pump
                            error_prev = error
                            if error < 0.0:
                                m_dot_pump = min(m_dot_pump + 1.0, m_dot_pump_max)
                            else:
                                m_dot_pump = max(m_dot_pump - 1.0, m_dot_pump_min)

        #### COMPLETE PUMP CALCS ####
        if m_dot_pump != m_dot_pump_min:
            h_pump_in = self.h_tank_start.v
            P_pump_in = P_tank
            flow = m_dot_pump / rho_water
            Eta_pump = max(
                self.Eta_Coef_A.v * flow ** 4.0 +
                self.Eta_Coef_B.v * flow ** 3.0 +
                self.Eta_Coef_C.v * flow ** 2.0 +
                self.Eta_Coef_D.v * flow,
                0.2
            )
            # CONVERTED-NEEDS UNITS CHECK: removed unit scale factors; eeslib uses SI
            s_pump_in = fp.entropy("water", P=P_pump_in, h=h_pump_in)
            s_pump_out_s = s_pump_in
            # CONVERTED-NEEDS UNITS CHECK: removed unit scale factors; eeslib uses SI
            h_pump_out_s = fp.enthalpy("water", P=P_pump_out, s=s_pump_out_s)
            W_dot_pump_s = m_dot_pump * (h_pump_out_s - h_pump_in)
            W_dot_pump = W_dot_pump_s / Eta_pump
            if m_dot_pump == 0.0:
                m_dot_pump = 0.000001
            h_pump_out = (h_pump_in * m_dot_pump + W_dot_pump) / m_dot_pump
        else:
            P_pump_out = self.P_fw_in.v
            h_pump_out = self.h_fw_in.v
            W_dot_pump = 0.0
            Eta_pump = 0.0

        #### Complete mixing tank calcs ####
        m_tank_end = self.m_tank_start.v + self.m_dot_BFWH.v * ts - m_dot_pump * ts
        L_tank_end = m_tank_end / rho_water / (3.14 * (self.D_tank.v / 2) ** 2.0)
        h_tank_end = (self.m_tank_start.v * self.h_tank_start.v + self.m_dot_BFWH.v * ts * self.h_BFWH.v - m_dot_pump * ts * self.h_tank_start.v) / m_tank_end
        # CONVERTED-NEEDS UNITS CHECK: removed unit scale factors; eeslib uses SI
        T_tank_end = fp.temperature("water", P=P_tank, h=h_tank_end)

        #### COMPLETE FEEDWATER MIXING CALCS ####
        m_dot_fw_out = self.m_dot_fw_in.v + m_dot_pump
        P_fw_out = self.P_fw_in.v
        h_fw_out = (self.m_dot_fw_in.v * self.h_fw_in.v + m_dot_pump * h_pump_out) / m_dot_fw_out

        self.VP_output.v = self.VP_output.v
        self.m_dot_fw_out.v = m_dot_fw_out
        self.h_fw_out.v = h_fw_out
        self.P_fw_out.v = P_fw_out
        self.m_dot_pump.v = m_dot_pump
        self.P_pump_out.v = P_pump_out
        self.W_dot_pump.v = W_dot_pump
        self.L_tank_end.v = L_tank_end
        self.T_tank_end.v = T_tank_end
        self.m_tank_end.v = m_tank_end
        self.h_tank_end.v = h_tank_end
