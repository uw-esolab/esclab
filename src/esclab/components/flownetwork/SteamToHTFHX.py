"""Steam-to-HTF shell-and-tube heat exchanger component model (Type 6017).

Object: ESOL6017-STHX
Simulation Studio Model: ESOL6017-STHX
"""

import math

import numpy as np
from eeslib import fluid_properties as fp

from esclab.simulate import Component
from esclab.components.esol_properties import Incompressible as Inc
from esclab.components.flownetwork.esol6015_helpers import f_cp_water


class SteamToHTFHX(Component):
    """
    TRNSYS Type 6017: ESOL6017-STHX.

    Shell-and-tube heat exchanger transferring heat from an HTF (shell side)
    to feedwater / steam (tube side).  Uses the NTU-effectiveness method with
    an off-design UA scaling law to predict heat transfer at off-rated
    conditions.  Alarm and trip logic for high/low temperatures, flow rates,
    pressures, and one-minute heating rates are evaluated at the end of each
    converged timestep.

    Parameters
    ----------
    heat_transfer_rated : float
        Rated overall heat-transfer coefficient [W/(m^2·K)].
    m_dot_fw_rated : float
        Rated feedwater mass flow rate [kg/s].
    m_dot_htf_rated : float
        Rated HTF mass flow rate [kg/s].
    rated_exp : float
        Exponent for off-design UA scaling law [-].
    No_shell_passes : float
        Number of shell passes [-].
    No_tube_passes : float
        Number of tube passes [-].
    Length_HX : float
        Length of heat exchanger [m].
    Tube_OD : float
        Tube outer diameter [m].
    Tube_th : float
        Tube wall thickness [m].
    No_tubes : float
        Number of tubes [-].
    Fluid_ID : str
        HTF fluid identifier passed to esol_properties.Incompressible.
    Alarm_HTF_high_temp_in : float
        High HTF Temperature entering HX Alarm setpoint [K].
    Trip_HTF_high_temp_in : float
        High HTF Temperature entering HX Trip setpoint [K].
    Alarm_HTF_low_temp_in : float
        Low HTF Temperature entering HX Alarm setpoint [K].
    Trip_HTF_low_temp_in : float
        Low HTF Temperature entering HX Trip setpoint [K].
    Alarm_HTF_low_temp_out : float
        Low HTF Temperature leaving HX Alarm setpoint [K].
    Trip_HTF_low_temp_out : float
        Low HTF Temperature leaving HX Trip setpoint [K].
    Alarm_HTF_high_flow_in : float
        High HTF Flow entering HX Alarm setpoint [kg/s].
    Trip_HTF_high_flow_in : float
        High HTF Flow entering HX Trip setpoint [kg/s].
    Alarm_HTF_high_press_in : float
        High HTF Pressure entering HX Alarm setpoint [Pa].
    Trip_HTF_high_press_in : float
        High HTF Pressure entering HX Trip setpoint [Pa].
    Alarm_FW_low_temp_in : float
        Low FW Temperature entering HX Alarm setpoint [K].
    Trip_FW_low_temp_in : float
        Low FW Temperature entering HX Trip setpoint [K].
    Alarm_HR_HTF_in : float
        Alarm Condition for High Heating Rate of HTF Inlet [K/min].
    Trip_HR_HTF_in : float
        Trip Condition for High Heating Rate of HTF Inlet [K/min].
    Alarm_HR_fw_out : float
        Alarm Condition for High Heating Rate of FW Outlet Temp [K/min].
    Trip_HR_fw_out : float
        Trip Condition for High Heating Rate of FW Outlet Temp [K/min].

    Inputs
    ------
    m_dot_fw : float
        Mass flow rate of feedwater entering the tube side [kg/s].
    P_fw : float
        Pressure of feedwater entering the tube side [Pa].
    h_fw : float
        Specific enthalpy of feedwater entering the tube side [J/kg].
    m_dot_HTF : float
        Mass flow rate of HTF entering the shell side [kg/s].
    P_htf : float
        Pressure of HTF entering the shell side [Pa].
    T_htf : float
        Temperature of HTF entering the shell side [K].

    Outputs
    -------
    m_dot_fw_out : float
        Mass flow rate of feedwater leaving the tube side [kg/s].
    Vol_dot_fw : float
        Volumetric flow rate of feedwater leaving the tube side [m^3/s].
    P_fw_out : float
        Pressure of feedwater flow leaving the tube side [Pa].
    h_fw_out : float
        Enthalpy of feedwater leaving the tube side [J/kg].
    T_fw_out : float
        Temperature of feedwater leaving the tube side [K].
    m_dot_HTF_out : float
        Mass flow rate of HTF leaving the shell side [kg/s].
    Vol_dot_HTF_out : float
        Volumetric flow rate of HTF leaving the shell side [m^3/s].
    P_HTF_out : float
        Pressure of HTF leaving the shell side [Pa].
    T_HTF_out : float
        Temperature of HTF leaving the shell side [K].
    Q_dot_hx : float
        Heat transfer between feedwater and HTF [W].
    eta_OD : float
        Effectiveness of heat exchanger [-].
    HR_HTF_in : float
        Heating rate of HTF inlet over 1-minute window [K/min].
    HR_fw_out_val : float
        Heating rate of FW outlet temp over 1-minute window [K/min].
    Out14_Alarm_HTF_high_temp_in : float
        High HTF Temp in Alarm [0/1].
    Out15_Trip_HTF_high_temp_in : float
        High HTF Temp in Trip [0/1].
    Out16_Alarm_HTF_low_temp_in : float
        Low HTF Temp in Alarm [0/1].
    Out17_Trip_HTF_low_temp_in : float
        Low HTF Temp in Trip [0/1].
    Out18_Alarm_HTF_low_temp_out : float
        Low HTF Temp out Alarm [0/1].
    Out19_Trip_HTF_low_temp_out : float
        Low HTF Temp out Trip [0/1].
    Out20_Alarm_HTF_high_flow_in : float
        High HTF Flow entering HX Alarm [0/1].
    Out21_Trip_HTF_high_flow_in : float
        High HTF Flow entering HX Trip [0/1].
    Out22_Alarm_HTF_high_press_in : float
        High HTF Pressure entering HX Alarm [0/1].
    Out23_Trip_HTF_high_press_in : float
        High HTF Pressure entering HX Trip [0/1].
    Out24_Alarm_FW_low_temp_in : float
        Low FW Temperature entering HX Alarm [0/1].
    Out25_Trip_FW_low_temp_in : float
        Low FW Temperature entering HX Trip [0/1].
    Out26_Alarm_HR_HTF_in : float
        High Heating Rate Alarm for HTF Inlet over 1-minute time [0/1].
    Out27_Trip_HR_HTF_in : float
        High Heating Rate Trip for HTF Inlet over 1-minute time [0/1].
    Out28_Alarm_HR_fw_out : float
        High Heating Rate Alarm for FW Outlet Temp over 1-minute time [0/1].
    Out29_Trip_HR_fw_out : float
        High Heating Rate Trip for FW Outlet Temp over 1-minute time [0/1].
    """

    # *** Model Parameters ***
    #    PARAMETERS
    heat_transfer_rated = Component.Parameter()   # rated overall heat transfer coefficient [W/(m^2·K)]
    m_dot_fw_rated = Component.Parameter()        # rated feedwater mass flow rate [kg/s]
    m_dot_htf_rated = Component.Parameter()       # rated HTF mass flow rate [kg/s]
    rated_exp = Component.Parameter()             # exponent for off-design UA scaling [-]
    No_shell_passes = Component.Parameter()       # number of shell passes [-]
    No_tube_passes = Component.Parameter()        # number of tube passes [-]
    Length_HX = Component.Parameter()             # length of heat exchanger [m]
    Tube_OD = Component.Parameter()               # tube outer diameter [m]
    Tube_th = Component.Parameter()               # tube wall thickness [m]
    No_tubes = Component.Parameter()              # number of tubes [-]
    Fluid_ID = Component.Parameter()              # HTF fluid identifier string
    # Alarm/trip parameters (params 12-27)
    Alarm_HTF_high_temp_in = Component.Parameter()   # High HTF Temperature entering HX Alarm [K]
    Trip_HTF_high_temp_in = Component.Parameter()    # High HTF Temperature entering HX Trip [K]
    Alarm_HTF_low_temp_in = Component.Parameter()    # Low HTF Temperature entering HX Alarm [K]
    Trip_HTF_low_temp_in = Component.Parameter()     # Low HTF Temperature entering HX Trip [K]
    Alarm_HTF_low_temp_out = Component.Parameter()   # Low HTF Temperature leaving HX Alarm [K]
    Trip_HTF_low_temp_out = Component.Parameter()    # Low HTF Temperature leaving HX Trip [K]
    Alarm_HTF_high_flow_in = Component.Parameter()   # High HTF Flow entering HX Alarm [kg/s]
    Trip_HTF_high_flow_in = Component.Parameter()    # High HTF Flow entering HX Trip [kg/s]
    Alarm_HTF_high_press_in = Component.Parameter()  # High HTF Pressure entering HX Alarm [Pa]
    Trip_HTF_high_press_in = Component.Parameter()   # High HTF Pressure entering HX Trip [Pa]
    Alarm_FW_low_temp_in = Component.Parameter()     # Low FW Temperature entering HX Alarm [K]
    Trip_FW_low_temp_in = Component.Parameter()      # Low FW Temperature entering HX Trip [K]
    Alarm_HR_HTF_in = Component.Parameter()          # Alarm Condition for High Heating Rate of HTF Inlet [K/min]
    Trip_HR_HTF_in = Component.Parameter()           # Trip Condition for High Heating Rate of HTF Inlet [K/min]
    Alarm_HR_fw_out = Component.Parameter()          # Alarm Condition for High Heating Rate of FW Outlet [K/min]
    Trip_HR_fw_out = Component.Parameter()           # Trip Condition for High Heating Rate of FW Outlet [K/min]

    # *** Model Inputs ***
    m_dot_fw = Component.Input()   # mass flow rate of feedwater [kg/s]
    P_fw = Component.Input()       # pressure of feedwater [Pa]
    h_fw = Component.Input()       # specific enthalpy of feedwater [J/kg]
    m_dot_HTF = Component.Input()  # mass flow rate of HTF [kg/s]
    P_htf = Component.Input()      # pressure of HTF [Pa]
    T_htf = Component.Input()      # temperature of HTF [K]

    # *** Model Outputs ***
    m_dot_fw_out = Component.Output()              # output  1: mass flow rate of feedwater leaving the tube side [kg/s]
    Vol_dot_fw = Component.Output()                # output  2: volumetric flow rate of feedwater leaving the tube side [m^3/s]
    P_fw_out = Component.Output()                  # output  3: Pressure of feedwater flow leaving the tube side [Pa]
    h_fw_out = Component.Output()                  # output  4: Enthalpy of feedwater leaving the tube side [J/kg]
    T_fw_out = Component.Output()                  # output  5: Temperature of feedwater leaving the tube side [K]
    m_dot_HTF_out = Component.Output()             # output  6: mass flow rate of htf leaving the shell side [kg/s]
    Vol_dot_HTF_out = Component.Output()           # output  7: volumetric flow rate of htf leaving the shell side [m^3/s]
    P_HTF_out = Component.Output()                 # output  8: Pressure of htf leaving the shell side [Pa]
    T_HTF_out = Component.Output()                 # output  9: Temperature of htf leaving the shell side [K]
    Q_dot_hx = Component.Output()                  # output 10: Heat transfer between feedwater and htf [W]
    eta_OD = Component.Output()                    # output 11: effectiveness of heat exchanger [-]
    HR_HTF_in = Component.Output()                 # output 12: heating rate of HTF inlet [K/min]
    HR_fw_out_val = Component.Output()             # output 13: heating rate of FW outlet temp [K/min]
    Out14_Alarm_HTF_high_temp_in = Component.Output()   # output 14: High HTF Temp in Alarm
    Out15_Trip_HTF_high_temp_in = Component.Output()    # output 15: High HTF Temp in Trip
    Out16_Alarm_HTF_low_temp_in = Component.Output()    # output 16: Low HTF Temp in Alarm
    Out17_Trip_HTF_low_temp_in = Component.Output()     # output 17: Low HTF Temp in Trip
    Out18_Alarm_HTF_low_temp_out = Component.Output()   # output 18: Low HTF Temp out Alarm
    Out19_Trip_HTF_low_temp_out = Component.Output()    # output 19: Low HTF Temp out Trip
    Out20_Alarm_HTF_high_flow_in = Component.Output()   # output 20: Low HTF Flow entering HX Alarm
    Out21_Trip_HTF_high_flow_in = Component.Output()    # output 21: Low HTF Flow entering HX Trip
    Out22_Alarm_HTF_high_press_in = Component.Output()  # output 22: High HTF Pressure entering HX Alarm
    Out23_Trip_HTF_high_press_in = Component.Output()   # output 23: High HTF Pressure entering HX Trip
    Out24_Alarm_FW_low_temp_in = Component.Output()     # output 24: Low FW Temperature entering HX Alarm
    Out25_Trip_FW_low_temp_in = Component.Output()      # output 25: Low FW Temperature entering HX Trip
    Out26_Alarm_HR_HTF_in = Component.Output()          # output 26: High Heating Rate Alarm for HTF Inlet over 1 minute time
    Out27_Trip_HR_HTF_in = Component.Output()           # output 27: High Heating Rate Trip for HTF Inlet over 1 minute time
    Out28_Alarm_HR_fw_out = Component.Output()          # output 28: High Heating Rate Alarm for FW Outlet Temp over 1 minute time
    Out29_Trip_HR_fw_out = Component.Output()           # output 29: High Heating Rate Trip for FW Outlet Temp over 1 minute time

    # -----------------------------------------------------------------------------------------------------------------------
    def presim_setup(self, **kwargs):
        """First call of simulation: configure array sizes.
        Equivalent to Fortran getIsFirstCallofSimulation() block."""

        # The number of static variables sized at first call is 2*N_int
        # N_int = ceiling(60 / (Timestep [hr] * 3600)) = ceiling(steps per minute)
        # TODO-NEEDS CONVERSION REVIEW: verify model.settings.timestep units (hours vs seconds);
        # Fortran uses Timestep in hours: N = 60 / (Timestep * 3600)
        N = 60.0 / (self.model.settings.timestep * 3600.0)
        N_int = math.ceil(N)
        self._N_int = N_int

        # Heating Rate for HTF Inlet (indices 0..N_int-1) in static variables array
        # Heating Rate for FW Outlet (indices N_int..2*N_int-1) in static variables array
        self._static_T = np.zeros(2 * N_int)

        # Flag to track whether the static arrays have been initialized with
        # meaningful temperatures (done at first end-of-timestep).
        self._static_initialized = False

    # -----------------------------------------------------------------------------------------------------------------------
    def calculate(self):
        """Main calculation called every iteration."""

        # ------------------------------------------------------------------
        # First timestep: set initial output values; no heat transfer computed.
        # Equivalent to Fortran getIsStartTime() block.
        # ------------------------------------------------------------------
        if self.model.is_first_step:
            # Set the Initial Values of the Outputs (#,Value)
            self.m_dot_fw_out.v = self.m_dot_fw.v     # mass flow rate of feedwater leaving the tube side
            self.Vol_dot_fw.v = 0.0                    # volumetric flow rate of feedwater leaving the tube side
            self.P_fw_out.v = self.P_fw.v              # Pressure of feedwater flow leaving the tube side
            self.h_fw_out.v = self.h_fw.v              # Enthalpy of feedwater leaving the tube side
            self.T_fw_out.v = 0.0                      # Temperature of feedwater leaving the tube side
            self.m_dot_HTF_out.v = self.m_dot_HTF.v   # mass flow rate of htf leaving the shell side
            self.Vol_dot_HTF_out.v = 0.0               # volumetric flow rate of htf leaving the shell side
            self.P_HTF_out.v = self.P_htf.v            # Pressure of htf leaving the shell side
            self.T_HTF_out.v = self.T_htf.v            # Temperature of htf leaving the shell side
            self.Q_dot_hx.v = 0.0                      # Heat transfer between feedwater and htf
            self.eta_OD.v = 0.0                        # effectiveness of heat exchanger
            return

        # ------------------------------------------------------------------
        # Main iterative calculation (all subsequent timesteps and iterations)
        # ------------------------------------------------------------------

        T_htf_local = self.T_htf.v
        # Default temperature until actual temperature enters
        if T_htf_local == 0.0:
            T_htf_local = 500.0

        if self.m_dot_HTF.v > 0.01:
            if self.m_dot_fw.v > 0.01:
                if self.P_fw.v > 0.0:  # Hydraulic calculations are reasonable

                    # Saturation temperature and enthalpy of saturated liquid at feedwater pressure
                    # eeslib.fluid_properties equivalent of FIT_PQ for T_sat and h_sat_f at x=0
                    # CONVERTED-NEEDS UNITS CHECK: removed /1000 Pa->kPa and /1000 J/kg->kJ/kg and *1000 kJ/kg->J/kg; eeslib uses SI
                    T_sat = fp.t_sat("water", P=self.P_fw.v)
                    h_sat_f = fp.enthalpy("water", P=self.P_fw.v, x=0.0)

                    # Finding temperature of feedwater in based on pressure and enthalpy in
                    # eeslib.fluid_properties equivalent of FIT_PH for T(P,h)
                    # CONVERTED-NEEDS UNITS CHECK: removed /1000 Pa->kPa and /1000 J/kg->kJ/kg; eeslib uses SI
                    T_fw = fp.temperature("water", P=self.P_fw.v, h=self.h_fw.v)

                    # solve for specific heat values of fw and htf
                    # finding the specific heat of the htf entering heat exchanger
                    cp_htf = Inc.specheat(self.Fluid_ID.v, T=T_htf_local, P=self.P_htf.v)

                    if abs(T_fw - T_sat) > 1.0:
                        # specific heat of feedwater entering heat exchanger
                        cp_fw = f_cp_water(P=self.P_fw.v, T=T_fw)
                    else:
                        # T_fw is too close to saturation temperature to give accurate cp value
                        cp_fw = f_cp_water(P=self.P_fw.v, T=((T_htf_local + T_fw) / 2.0))

                    # Find effectiveness of heat exchanger based on inlet conditions of HTF and FW
                    # Surface area of the feedwater side of the heat exchanger
                    A_s = (3.14 * (self.Tube_OD.v - 2.0 * self.Tube_th.v)
                           * self.Length_HX.v * self.No_tube_passes.v * self.No_tubes.v)
                    UA_rated = self.heat_transfer_rated.v * A_s
                    # finding off design UA for heat exchanger
                    UA_OD = UA_rated * (self.m_dot_HTF.v / self.m_dot_htf_rated.v) ** self.rated_exp.v
                    Cap_fw = self.m_dot_fw.v * cp_fw
                    Cap_htf = self.m_dot_HTF.v * cp_htf
                    Cap_min = min(Cap_fw, Cap_htf)
                    CR = max(Cap_min / max(Cap_fw, Cap_htf), 0.001)
                    # Off design NTU for heat exchanger
                    NTU_OD = UA_OD / Cap_min
                    # effectiveness for one pass shell and tube heat exchanger
                    eta_1pass = (2.0 * (1.0 + CR + math.sqrt(1.0 + CR ** 2.0)
                                        * ((1.0 + math.exp(-NTU_OD * math.sqrt(1.0 + CR ** 2.0)))
                                           / (1.0 - math.exp(-NTU_OD * math.sqrt(1.0 + CR ** 2.0))))) ** (-1.0))
                    eta_od_calc = (
                        (((1.0 - eta_1pass * CR) / (1.0 - eta_1pass)) ** self.No_shell_passes.v - 1.0)
                        / (((1.0 - eta_1pass * CR) / (1.0 - eta_1pass)) ** self.No_shell_passes.v - CR)
                    )
                    # Check that HTF temp is higher than FW temp

                    if T_fw < T_htf_local:  # Heat Transfer is going the correct way

                        # eeslib.fluid_properties equivalent of FIT_TP for h(P,T)
                        # CONVERTED-NEEDS UNITS CHECK: removed /1000 Pa->kPa and *1000 kJ/kg->J/kg; eeslib uses SI
                        h_fw_out_s = fp.enthalpy("water", P=self.P_fw.v, T=T_htf_local)

                        Q_dot_hx_calc = min(
                            eta_od_calc * self.m_dot_fw.v * (h_fw_out_s - self.h_fw.v),
                            eta_od_calc * self.m_dot_HTF.v * cp_htf * (T_htf_local - T_fw)
                        )
                        h_fw_out_calc = max(
                            (self.m_dot_fw.v * self.h_fw.v + Q_dot_hx_calc) / self.m_dot_fw.v,
                            self.h_fw.v
                        )
                        # Flow is entering subcooled; make sure it is not passing T_sat
                        if self.h_fw.v < h_sat_f:
                            if h_fw_out_calc > h_sat_f:
                                h_fw_out_calc = h_sat_f
                                Q_dot_hx_calc = self.m_dot_fw.v * (h_fw_out_calc - self.h_fw.v)

                        # eeslib.fluid_properties equivalent of FIT_PH for T(P,h)
                        # CONVERTED-NEEDS UNITS CHECK: removed /1000 Pa->kPa and /1000 J/kg->kJ/kg; eeslib uses SI
                        T_fw_out_calc = fp.temperature("water", P=self.P_fw.v, h=h_fw_out_calc)

                        # energy balance on HTF side of heat exchanger
                        T_htf_out_calc = (
                            (self.m_dot_HTF.v * cp_htf * T_htf_local) - Q_dot_hx_calc
                        ) / (self.m_dot_HTF.v * cp_htf)

                    else:
                        # htf temperature is lower than feedwater temperature,
                        # heat transfer is going the wrong way
                        # eeslib.fluid_properties equivalent of FIT_TP for h(P,T)
                        # CONVERTED-NEEDS UNITS CHECK: removed /1000 Pa->kPa and *1000 kJ/kg->J/kg; eeslib uses SI
                        h_fw_out_s = fp.enthalpy("water", P=self.P_fw.v, T=T_htf_local)

                        Q_dot_hx_calc = min(
                            eta_od_calc * self.m_dot_fw.v * (self.h_fw.v - h_fw_out_s),
                            eta_od_calc * self.m_dot_HTF.v * cp_htf * (T_fw - T_htf_local)
                        )
                        if self.m_dot_fw.v < 1.0:
                            Q_dot_hx_calc = 0.0
                            h_fw_out_calc = self.h_fw.v   # enthalpy out = enthalpy in
                            T_fw_out_calc = T_fw
                        else:
                            h_fw_out_calc = (self.m_dot_fw.v * self.h_fw.v - Q_dot_hx_calc) / self.m_dot_fw.v

                            # eeslib.fluid_properties equivalent of FIT_PQ for h_sat_g at x=1
                            # CONVERTED-NEEDS UNITS CHECK: removed /1000 Pa->kPa and *1000 kJ/kg->J/kg; eeslib uses SI
                            h_sat_g = fp.enthalpy("water", P=self.P_fw.v, x=1.0)

                            # if this type is the superheater do not allow it to go below saturation
                            # - would break other types
                            if self.h_fw.v >= h_sat_g:
                                h_fw_out_calc = h_sat_g

                            # eeslib.fluid_properties equivalent of FIT_PH for T(P,h)
                            # CONVERTED-NEEDS UNITS CHECK: removed /1000 Pa->kPa and /1000 J/kg->kJ/kg; eeslib uses SI
                            T_fw_out_calc = fp.temperature("water", P=self.P_fw.v, h=self.h_fw.v)

                        # energy balance on HTF side of heat exchanger
                        T_htf_out_calc = (
                            (self.m_dot_HTF.v * cp_htf * T_htf_local) + Q_dot_hx_calc
                        ) / (self.m_dot_HTF.v * cp_htf)

                    # calculating the volumetric flow rates
                    # eeslib.fluid_properties equivalent of FIT_PH for density(P,h)
                    # CONVERTED-NEEDS UNITS CHECK: removed /1000 Pa->kPa and /1000 J/kg->kJ/kg; eeslib uses SI
                    rho_fw = fp.density("water", P=self.P_fw.v, h=h_fw_out_calc)
                    if rho_fw > 0.0:
                        Vol_dot_fw_calc = self.m_dot_fw.v / rho_fw
                    else:
                        Vol_dot_fw_calc = 0.0

                    rho_htf = Inc.density(self.Fluid_ID.v, T=T_htf_out_calc, P=self.P_htf.v)
                    Vol_dot_htf_calc = self.m_dot_HTF.v / rho_htf

                    # Set the Outputs from this Model (#,Value)
                    self.m_dot_fw_out.v = self.m_dot_fw.v     # mass flow rate of feedwater leaving the tube side
                    self.Vol_dot_fw.v = Vol_dot_fw_calc        # volumetric flow rate of feedwater leaving the tube side
                    self.P_fw_out.v = self.P_fw.v              # Pressure of feedwater flow leaving the tube side
                    self.h_fw_out.v = h_fw_out_calc            # Enthalpy of feedwater leaving the tube side
                    self.T_fw_out.v = T_fw_out_calc            # Temperature of feedwater leaving the tube side
                    self.m_dot_HTF_out.v = self.m_dot_HTF.v   # mass flow rate of htf leaving the shell side
                    self.Vol_dot_HTF_out.v = Vol_dot_htf_calc  # volumetric flow rate of htf leaving the shell side
                    self.P_HTF_out.v = self.P_htf.v            # Pressure of htf leaving the shell side
                    self.T_HTF_out.v = T_htf_out_calc          # Temperature of htf leaving the shell side
                    self.Q_dot_hx.v = Q_dot_hx_calc            # Heat transfer between feedwater and htf
                    self.eta_OD.v = eta_od_calc                 # effectiveness of heat exchanger

                else:
                    # Pressure is not possible, need to wait for next iteration to compute temperatures
                    # keep values the same
                    # eeslib.fluid_properties equivalent of FIT_PH for T(P,h)
                    # CONVERTED-NEEDS UNITS CHECK: removed /1000 Pa->kPa and /1000 J/kg->kJ/kg; eeslib uses SI
                    T_fw = fp.temperature("water", P=self.P_fw.v, h=self.h_fw.v)

                    self.m_dot_fw_out.v = self.m_dot_fw.v    # m_dot_fw
                    self.P_fw_out.v = self.P_fw.v            # P_fw
                    self.h_fw_out.v = self.h_fw.v            # h_fw
                    self.T_fw_out.v = T_fw                   # T_fw
                    self.m_dot_HTF_out.v = self.m_dot_HTF.v  # m_dot_HTF
                    self.P_HTF_out.v = self.P_htf.v          # P_HTF
                    self.T_HTF_out.v = T_htf_local           # T_HTF
                    self.Q_dot_hx.v = 0.0                    # Q_dot_HX
                    self.eta_OD.v = 0.0                      # Epsilon_HX

            else:
                # no FW Flow entering the system; set HTF inputs as the same
                # HTF Outputs
                T_htf_out_calc = T_htf_local
                rho_HTF = Inc.density(self.Fluid_ID.v, T=T_htf_local, P=self.P_htf.v)
                Vol_dot_HTF_calc = self.m_dot_HTF.v / rho_HTF

                # FW Outputs
                # eeslib.fluid_properties equivalent of FIT_PH for T(P,h)
                # CONVERTED-NEEDS UNITS CHECK: removed /1000 Pa->kPa and /1000 J/kg->kJ/kg; eeslib uses SI
                T_fw_out_calc = fp.temperature("water", P=self.P_fw.v, h=self.h_fw.v)
                Vol_dot_fw_calc = self.m_dot_fw.v
                h_fw_out_calc = self.h_fw.v

                Q_dot_hx_calc = 0.0
                eta_od_calc = 0.0

                self.m_dot_fw_out.v = self.m_dot_fw.v         # mass flow rate of feedwater leaving the tube side
                self.Vol_dot_fw.v = Vol_dot_fw_calc            # volumetric flow rate of feedwater leaving the tube side
                self.P_fw_out.v = self.P_fw.v                  # Pressure of feedwater flow leaving the tube side
                self.h_fw_out.v = h_fw_out_calc                # Enthalpy of feedwater leaving the tube side
                self.T_fw_out.v = T_fw_out_calc                # Temperature of feedwater leaving the tube side
                self.m_dot_HTF_out.v = self.m_dot_HTF.v       # mass flow rate of htf leaving the shell side
                self.Vol_dot_HTF_out.v = Vol_dot_HTF_calc      # volumetric flow rate of htf leaving the shell side
                self.P_HTF_out.v = self.P_htf.v                # Pressure of htf leaving the shell side
                self.T_HTF_out.v = T_htf_out_calc              # Temperature of htf leaving the shell side
                self.Q_dot_hx.v = Q_dot_hx_calc                # Heat transfer between feedwater and htf
                self.eta_OD.v = eta_od_calc                     # effectiveness of heat exchanger

        else:
            # no HTF Flow entering the system; set feedwater outlet the same as the inlet
            # Set the Outputs from this Model (#,Value)
            # HTF Outputs
            Vol_dot_HTF_calc = 0.0
            T_htf_out_calc = T_htf_local

            # FW Outputs
            # eeslib.fluid_properties equivalent of FIT_PH for T(P,h) and density(P,h)
            # CONVERTED-NEEDS UNITS CHECK: removed /1000 Pa->kPa and /1000 J/kg->kJ/kg; eeslib uses SI
            rho_fw = fp.density("water", P=self.P_fw.v, h=self.h_fw.v)
            T_fw_out_calc = fp.temperature("water", P=self.P_fw.v, h=self.h_fw.v)
            if rho_fw > 0.0:
                Vol_dot_fw_calc = self.m_dot_fw.v / rho_fw
            else:
                Vol_dot_fw_calc = 0.0  # Imitate as bad sensor
            h_fw_out_calc = self.h_fw.v

            Q_dot_hx_calc = 0.0
            eta_od_calc = 0.0

            self.m_dot_fw_out.v = self.m_dot_fw.v         # mass flow rate of feedwater leaving the tube side
            self.Vol_dot_fw.v = Vol_dot_fw_calc            # volumetric flow rate of feedwater leaving the tube side
            self.P_fw_out.v = self.P_fw.v                  # Pressure of feedwater flow leaving the tube side
            self.h_fw_out.v = h_fw_out_calc                # Enthalpy of feedwater leaving the tube side
            self.T_fw_out.v = T_fw_out_calc                # Temperature of feedwater leaving the tube side
            self.m_dot_HTF_out.v = self.m_dot_HTF.v       # mass flow rate of htf leaving the shell side
            self.Vol_dot_HTF_out.v = Vol_dot_HTF_calc      # volumetric flow rate of htf leaving the shell side
            self.P_HTF_out.v = self.P_htf.v                # Pressure of htf leaving the shell side
            self.T_HTF_out.v = T_htf_out_calc              # Temperature of htf leaving the shell side
            self.Q_dot_hx.v = Q_dot_hx_calc                # Heat transfer between feedwater and htf
            self.eta_OD.v = eta_od_calc                     # effectiveness of heat exchanger

    # -----------------------------------------------------------------------------------------------------------------------
    def converged(self):
        """End-of-timestep operations after convergence.
        Equivalent to Fortran getIsEndOfTimestep() block.
        Computes heating rates and evaluates alarm/trip conditions."""

        N_int = self._N_int

        # TODO-NEEDS CONVERSION REVIEW: verify model.settings.timestep units (hours vs seconds);
        # Fortran Timestep is in hours. N = 60 / (Timestep * 3600) where Timestep is hours.
        timestep_hrs = self.model.settings.timestep  # assumed hours to match Fortran convention

        # At beginning of simulation set initial temperature arrays for Heating Rates
        if not self._static_initialized:
            # Heating Rate for HTF Inlet (indices 0..N_int-1) in Static Variables
            # Heating Rate for FW Outlet (indices N_int..2*N_int-1) in Static Variables
            T_HTF_init = self.T_htf.v
            T_fw_out_init = self.T_fw_out.v
            for i in range(N_int):
                self._static_T[i] = T_HTF_init
                self._static_T[N_int + i] = T_fw_out_init
            self._static_initialized = True

        # !!!!!High or Low Temperature Alarms!!!!
        T_HTF = self.T_htf.v          # HTF Inlet Temperature
        T_HTF_out = self.T_HTF_out.v  # HTF Outlet Temperature
        m_dot_HTF_val = self.m_dot_HTF.v   # HTF Flow Rate
        P_HTF_val = self.P_htf.v           # HTF inlet Pressure

        # eeslib.fluid_properties equivalent of FIT_PH for T(P,h)
        # CONVERTED-NEEDS UNITS CHECK: removed /1000 Pa->kPa and /1000 J/kg->kJ/kg; eeslib uses SI
        T_fw = fp.temperature("water", P=self.P_fw.v, h=self.h_fw.v)

        # High HTF Temp In Check
        if self.Alarm_HTF_high_temp_in.v > T_HTF:  # HTF Inlet Temperature is less than Alarm (Alarm = 0)
            Alarm = 0.0
            Trip = 0.0
        else:
            Alarm = 1.0
            if self.Trip_HTF_high_temp_in.v > T_HTF:  # HTF Inlet Temperature is less than Trip
                Trip = 0.0
            else:
                Trip = 1.0
        self.Out14_Alarm_HTF_high_temp_in.v = Alarm   # High HTF Temp in Alarm
        self.Out15_Trip_HTF_high_temp_in.v = Trip     # High HTF Temp in Trip

        if self.Alarm_HTF_low_temp_in.v < T_HTF:  # HTF Inlet Temperature is above low temp alarm (Alarm = 0)
            Alarm = 0.0
            Trip = 0.0
        else:
            Alarm = 1.0
            if self.Trip_HTF_low_temp_in.v < T_HTF:  # HTF Inlet Temperature is lower than alarm but higher than trip
                Trip = 0.0
            else:
                Trip = 1.0
        self.Out16_Alarm_HTF_low_temp_in.v = Alarm   # Low HTF Temp in Alarm
        self.Out17_Trip_HTF_low_temp_in.v = Trip     # Low HTF Temp in Trip

        if self.Alarm_HTF_low_temp_out.v < T_HTF_out:
            Alarm = 0.0
            Trip = 0.0
        else:
            Alarm = 1.0
            if self.Trip_HTF_low_temp_out.v < T_HTF_out:
                Trip = 0.0
            else:
                Trip = 1.0
        self.Out18_Alarm_HTF_low_temp_out.v = Alarm   # Low HTF Temp out Alarm
        self.Out19_Trip_HTF_low_temp_out.v = Trip     # Low HTF Temp out Trip

        if m_dot_HTF_val < self.Alarm_HTF_high_flow_in.v:
            Alarm = 0.0
            Trip = 0.0
        else:
            Alarm = 1.0
            if m_dot_HTF_val < self.Trip_HTF_high_flow_in.v:
                Trip = 0.0
            else:
                Trip = 1.0
        self.Out20_Alarm_HTF_high_flow_in.v = Alarm   # Low HTF Flow entering HX Alarm
        self.Out21_Trip_HTF_high_flow_in.v = Trip     # Low HTF Flow entering HX Trip

        if P_HTF_val < self.Alarm_HTF_high_press_in.v:
            Alarm = 0.0
            Trip = 0.0
        else:
            Alarm = 1.0
            if P_HTF_val < self.Trip_HTF_high_press_in.v:
                Trip = 0.0
            else:
                Trip = 1.0
        self.Out22_Alarm_HTF_high_press_in.v = Alarm   # High HTF Pressure entering HX Alarm
        self.Out23_Trip_HTF_high_press_in.v = Trip     # High HTF Pressure entering HX Trip

        if T_fw > self.Alarm_FW_low_temp_in.v:
            Alarm = 0.0
            Trip = 0.0
        else:
            Alarm = 1.0
            if T_fw > self.Trip_FW_low_temp_in.v:
                Trip = 0.0
            else:
                Trip = 1.0
        self.Out24_Alarm_FW_low_temp_in.v = Alarm   # Low FW Temperature entering HX Alarm
        self.Out25_Trip_FW_low_temp_in.v = Trip     # Low FW Temperature entering HX Trip

        # !!!!!!High Heating Rates Alarms and Trips!!!!!
        # HTF Inlet Heating Rate stored as a static variable between index 0 and N_int-1

        T_arr = np.zeros(N_int)

        # --- HTF Inlet Heating Rate ---
        HR = 0.0
        for i in range(N_int):
            # Location where HTF Inlet Temperatures for past minute are stored (Fortran indices 1..N_int)
            T_arr[i] = self._static_T[i]
            if i > 0:
                HR = HR + (T_arr[i] - T_arr[i - 1]) / (timestep_hrs * 3600.0)
        HR = HR + (T_HTF - T_arr[N_int - 1]) / (timestep_hrs * 3600.0)
        # Divide by number of HR solved for and multiply by 60 to get K/min
        HR = HR / N_int * 60.0

        if abs(HR) < self.Alarm_HR_HTF_in.v:  # No Alarm Condition
            Alarm = 0.0
            Trip = 0.0
        else:  # Alarm Condition
            Alarm = 1.0
            if abs(HR) < self.Trip_HR_HTF_in.v:  # No Trip Condition
                Trip = 0.0
            else:
                Trip = 1.0
        self.HR_HTF_in.v = HR
        self.Out26_Alarm_HR_HTF_in.v = Alarm   # High Heating Rate Alarm for HTF Inlet over 1 minute time
        self.Out27_Trip_HR_HTF_in.v = Trip     # High Heating Rate Trip for HTF Inlet over 1 minute time

        # Update static temperature array for HTF inlet (shift left, append current)
        for i in range(N_int):
            if i != N_int - 1:
                self._static_T[i] = T_arr[i + 1]
            else:
                self._static_T[i] = T_HTF

        # --- Outlet FW Heating Rate ---
        T_fw_out_now = self.T_fw_out.v

        HR = 0.0
        for i in range(N_int):
            # Location where Feedwater Outlet Temperature for past minute are stored
            T_arr[i] = self._static_T[N_int + i]
            if i > 0:
                HR = HR + (T_arr[i] - T_arr[i - 1]) / timestep_hrs
        # TODO-NEEDS CONVERSION REVIEW: Fortran uses T_HTF (not T_fw_out) here – possible bug in original source
        HR = HR + (T_HTF - T_arr[N_int - 1]) / timestep_hrs
        # Divide by number of HR solved for and multiply by 60 to get K/min
        HR = HR / N_int * 60

        if abs(HR) < self.Alarm_HR_fw_out.v:  # No Alarm Condition
            Alarm = 0.0
            Trip = 0.0
        else:  # Alarm Condition
            Alarm = 1.0
            if abs(HR) < self.Trip_HR_fw_out.v:  # No Trip Condition
                Trip = 0.0
            else:
                Trip = 1.0
        self.HR_fw_out_val.v = HR
        self.Out28_Alarm_HR_fw_out.v = Alarm   # High Heating Rate Alarm for FW Outlet Temp over 1 minute time
        self.Out29_Trip_HR_fw_out.v = Trip     # High Heating Rate Trip for FW Outlet Temp over 1 minute time

        # Update Temperature variables for next end of timestep
        # TODO-NEEDS CONVERSION REVIEW: Fortran stores T_HTF (not T_fw_out_now) in FW outlet history –
        #   possible bug in original source; porting as-is.
        for i in range(N_int):
            if i != N_int - 1:
                self._static_T[N_int + i] = T_arr[i + 1]
            else:
                self._static_T[N_int + i] = T_HTF
