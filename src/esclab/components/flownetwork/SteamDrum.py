"""Steam drum with integrated evaporator component model (Type 6003)."""

import math

import numpy as np
from eeslib import fluid_properties as fp

from esclab.simulate import Component
from esclab.components.flownetwork.esol6015_helpers import drhodhcp, drhodpch, dudhcp, dudpch, tank_level


class SteamDrum(Component):
    """
    Object: Steam Drum
    Simulation Studio Model: ESOL6003-SteamDrum

    Models a steam drum with an integrated evaporator (shell-and-tube HX).
    Feedwater enters the drum and is boiled by heat transferred from an HTF
    stream.  The drum pressure and enthalpy are advanced in time using a
    4th-order Runge-Kutta integration.  Alarm and trip logic monitors drum
    level, pressure, temperature differentials, and heating rates.

    Parameters
    ----------
    D_tank : float
        Tank (drum) diameter [m].
    Length_tank : float
        Tank (drum) length [m].
    L_tank_init : float
        Initial drum water level [m].
    P_tank_init : float
        Initial drum pressure [Pa].
    Rated_heat_transfer : float
        Rated heat-transfer coefficient of the evaporator [W/m²/K].
    Rated_htf_flow : float
        Rated HTF mass flow rate [kg/s].
    Rated_exp : float
        Exponent for off-design UA correction [-].
    No_shell_passes : float
        Number of shell passes in the evaporator [-].
    No_tube_passes : float
        Number of tube passes in the evaporator [-].
    Length_HX : float
        Tube length of the evaporator [m].
    Tube_OD : float
        Tube outer diameter [m].
    Tube_th : float
        Tube wall thickness [m].
    No_tubes : float
        Number of tubes in the evaporator [-].
    Fluid_ID : str
        HTF fluid identifier for property lookups.
    LowLevel_Alarm_cond : float
        Low-level alarm setpoint [m].
    LowLevel_Trip_cond : float
        Low-level trip setpoint [m].
    HighLevel_Alarm_cond : float
        High-level alarm setpoint [m].
    HighLevel_Trip_cond : float
        High-level trip setpoint [m].
    HighPressure_Alarm_cond : float
        High-pressure alarm setpoint [Pa].
    HighPressure_Trip_cond : float
        High-pressure trip setpoint [Pa].
    HighDeltaT_Alarm_cond : float
        High |T_fw_in - T_sat| alarm setpoint [K].
    HighDeltaT_Trip_cond : float
        High |T_fw_in - T_sat| trip setpoint [K].
    HighHR_Alarm_cond : float
        High drum heating-rate alarm setpoint [K/min].
    HighHR_Trip_cond : float
        High drum heating-rate trip setpoint [K/min].
    HighTempDiffEvap_Alarm_cond : float
        High (T_HTF_in - T_sat) alarm setpoint [K].
    HighTempDiffEvap_Trip_cond : float
        High (T_HTF_in - T_sat) trip setpoint [K].
    HighHTF_HR_Alarm_cond : float
        High HTF inlet heating-rate alarm setpoint [K/min].
    HighHTF_HR_Trip_cond : float
        High HTF inlet heating-rate trip setpoint [K/min].

    Inputs
    ------
    m_dot_in : float
        Feedwater mass flow rate entering the drum [kg/s].
    h_in : float
        Feedwater enthalpy entering the drum [J/kg].
    P_in : float
        Feedwater pressure entering the drum [Pa].
    HTF_mass_in : float
        HTF mass flow rate entering the evaporator [kg/s].
    HTF_Temp_in : float
        HTF temperature entering the evaporator [K].
    HTF_P_in : float
        HTF pressure entering the evaporator [Pa].
    m_dot_superheat : float
        Steam mass flow rate leaving the drum (feedback/iteration input) [kg/s].

    Outputs
    -------
    m_dot_steam : float
        Mass flow rate of steam leaving the steam drum [kg/s].
    Vol_dot_fw : float
        Volumetric flow rate of steam leaving the steam drum [m³/s].
    P_steam_out : float
        Pressure of steam leaving the steam drum [Pa].
    T_steam_out : float
        Temperature of steam leaving the steam drum [K].
    h_steam_out : float
        Enthalpy of steam leaving the steam drum [J/kg].
    HTF_mass_out : float
        HTF mass flow rate leaving the evaporator [kg/s].
    Vol_dot_HTF : float
        HTF volumetric flow rate leaving the evaporator [m³/s].
    HTF_Temp_out : float
        HTF temperature leaving the evaporator [K].
    HTF_P_out : float
        HTF pressure leaving the evaporator [Pa].
    m_tank_new : float
        Total mass in the drum at the end of the current timestep [kg].
    L_tank_new : float
        Drum water level at the end of the current timestep [m].
    P_tank_new : float
        Drum pressure at the end of the current timestep [Pa].
    T_tank_new : float
        Drum saturation temperature at the end of the current timestep [K].
    h_tank_new : float
        Drum mixture enthalpy at the end of the current timestep [J/kg].
    m_tank_prev : float
        Total mass in the drum at the beginning of the current timestep [kg].
    L_tank_prev : float
        Drum water level at the beginning of the current timestep [m].
    P_tank_prev : float
        Drum pressure at the beginning of the current timestep [Pa].
    h_tank_prev : float
        Drum mixture enthalpy at the beginning of the current timestep [J/kg].
    eta_OD : float
        Off-design effectiveness of the evaporator [-].
    m_dot_evap : float
        Rate of water evaporation in the drum during the timestep [kg/s].
    HR_drum : float
        Drum temperature heating rate [K/min].
    HR_HTF : float
        HTF inlet heating rate [K/min].
    LowLevel_Alarm : float
        Low-level alarm signal (0 = no alarm, 1 = alarm).
    LowLevel_Trip : float
        Low-level trip signal (0 = no trip, 1 = trip).
    HighLevel_Alarm : float
        High-level alarm signal.
    HighLevel_Trip : float
        High-level trip signal.
    HighPressure_Alarm : float
        High-pressure alarm signal.
    HighPressure_Trip : float
        High-pressure trip signal.
    HighDeltaT_Alarm : float
        High |T_fw_in - T_sat| alarm signal.
    HighDeltaT_Trip : float
        High |T_fw_in - T_sat| trip signal.
    HighHR_Alarm : float
        High drum heating-rate alarm signal.
    HighHR_Trip : float
        High drum heating-rate trip signal.
    HighTempDiffEvap_Alarm : float
        High (T_HTF_in - T_sat) alarm signal.
    HighTempDiffEvap_Trip : float
        High (T_HTF_in - T_sat) trip signal.
    HighHTF_HR_Alarm : float
        High HTF inlet heating-rate alarm signal.
    HighHTF_HR_Trip : float
        High HTF inlet heating-rate trip signal.
    """

    # *** Model Parameters ***
    D_tank = Component.Parameter()                  # Tank Diameter [m]
    Length_tank = Component.Parameter()             # Tank Length [m]
    L_tank_init = Component.Parameter()             # Initial Level of the tank [m]
    P_tank_init = Component.Parameter()             # Initial Pressure in the tank [Pa]
    Rated_heat_transfer = Component.Parameter()     # Rated heat transfer coefficient [W/m²/K]
    Rated_htf_flow = Component.Parameter()          # Rated HTF flow rate [kg/s]
    Rated_exp = Component.Parameter()               # Off-design UA exponent [-]
    No_shell_passes = Component.Parameter()         # Number of shell passes [-]
    No_tube_passes = Component.Parameter()          # Number of tube passes [-]
    Length_HX = Component.Parameter()               # Heat exchanger tube length [m]
    Tube_OD = Component.Parameter()                 # Tube outer diameter [m]
    Tube_th = Component.Parameter()                 # Tube wall thickness [m]
    No_tubes = Component.Parameter()                # Number of tubes [-]
    Fluid_ID = Component.Parameter()                # HTF fluid identifier [-]
    LowLevel_Alarm_cond = Component.Parameter()     # Low level alarm condition [m]
    LowLevel_Trip_cond = Component.Parameter()      # Low level trip condition [m]
    HighLevel_Alarm_cond = Component.Parameter()    # High level alarm condition [m]
    HighLevel_Trip_cond = Component.Parameter()     # High level trip condition [m]
    HighPressure_Alarm_cond = Component.Parameter() # High pressure alarm condition [Pa]
    HighPressure_Trip_cond = Component.Parameter()  # High pressure trip condition [Pa]
    HighDeltaT_Alarm_cond = Component.Parameter()   # High delta T alarm condition [K]
    HighDeltaT_Trip_cond = Component.Parameter()    # High delta T trip condition [K]
    HighHR_Alarm_cond = Component.Parameter()       # High HR drum alarm condition [K/min]
    HighHR_Trip_cond = Component.Parameter()        # High HR drum trip condition [K/min]
    HighTempDiffEvap_Alarm_cond = Component.Parameter()  # High temp diff evap alarm condition [K]
    HighTempDiffEvap_Trip_cond = Component.Parameter()   # High temp diff evap trip condition [K]
    HighHTF_HR_Alarm_cond = Component.Parameter()   # High HTF HR alarm condition [K/min]
    HighHTF_HR_Trip_cond = Component.Parameter()    # High HTF HR trip condition [K/min]

    # *** Model Inputs ***
    m_dot_in = Component.Input()        # feedwater mass flow rate entering the drum [kg/s]
    h_in = Component.Input()             # feedwater enthalpy entering the drum [J/kg]
    P_in = Component.Input()             # feedwater pressure entering the drum [Pa]
    HTF_mass_in = Component.Input()      # HTF mass flow rate entering the evaporator [kg/s]
    HTF_Temp_in = Component.Input()      # HTF temperature entering the evaporator [K]
    HTF_P_in = Component.Input()         # HTF pressure entering the evaporator [Pa]
    m_dot_superheat = Component.Input()  # steam mass flow rate leaving the drum (feedback) [kg/s]

    # *** Model Outputs ***
    m_dot_steam = Component.Output()    # mass flow rate of steam leaving the steam drum [kg/s]
    Vol_dot_fw = Component.Output()     # Volumetric flow rate of steam leaving the steam drum [m³/s]
    P_steam_out = Component.Output()    # Pressure of steam leaving the steam drum [Pa]
    T_steam_out = Component.Output()    # Temperature of steam leaving the steam drum [K]
    h_steam_out = Component.Output()    # Enthalpy of steam leaving the steam drum [J/kg]
    HTF_mass_out = Component.Output()   # HTF mass flow rate leaving the evaporator [kg/s]
    Vol_dot_HTF = Component.Output()    # HTF volumetric flow rate leaving the evaporator [m³/s]
    HTF_Temp_out = Component.Output()   # HTF temperature leaving the evaporator [K]
    HTF_P_out = Component.Output()      # HTF pressure leaving the evaporator [Pa]
    m_tank_new = Component.Output()     # total mass in the drum at end of timestep [kg]
    L_tank_new = Component.Output()     # drum water level at end of timestep [m]
    P_tank_new = Component.Output()     # drum pressure at end of timestep [Pa]
    T_tank_new = Component.Output()     # drum saturation temperature at end of timestep [K]
    h_tank_new = Component.Output()     # drum mixture enthalpy at end of timestep [J/kg]
    m_tank_prev = Component.Output()    # total mass in the drum at beginning of timestep [kg]
    L_tank_prev = Component.Output()    # drum water level at beginning of timestep [m]
    P_tank_prev = Component.Output()    # drum pressure at beginning of timestep [Pa]
    h_tank_prev = Component.Output()    # drum mixture enthalpy at beginning of timestep [J/kg]
    eta_OD = Component.Output()         # off-design effectiveness of the evaporator [-]
    m_dot_evap = Component.Output()     # rate of water evaporation during the timestep [kg/s]
    HR_drum = Component.Output()        # Heating Rate of drum [K/min]
    HR_HTF = Component.Output()         # Heating Rate of HTF at evaporator inlet [K/min]
    LowLevel_Alarm = Component.Output()             # Low Level Alarm [-]
    LowLevel_Trip = Component.Output()              # Low Level Trip [-]
    HighLevel_Alarm = Component.Output()            # High Level Alarm [-]
    HighLevel_Trip = Component.Output()             # High Level Trip [-]
    HighPressure_Alarm = Component.Output()         # High Pressure Alarm [-]
    HighPressure_Trip = Component.Output()          # High Pressure Trip [-]
    HighDeltaT_Alarm = Component.Output()           # High Temp Diff Alarm (|T_fw_in - T_sat|) [-]
    HighDeltaT_Trip = Component.Output()            # High Temp Diff Trip (|T_fw_in - T_sat|) [-]
    HighHR_Alarm = Component.Output()               # Alarm for High Heating rate in the steam drum [-]
    HighHR_Trip = Component.Output()                # Trip for High Heating Rate in the Steam Drum [-]
    HighTempDiffEvap_Alarm = Component.Output()     # Alarm for High Temp Diff between evap HTF in and T_sat [-]
    HighTempDiffEvap_Trip = Component.Output()      # Trip for High Temp Diff between evap HTF in and T_sat [-]
    HighHTF_HR_Alarm = Component.Output()           # Alarm for High Heating rate at evaporator inlet [-]
    HighHTF_HR_Trip = Component.Output()            # Trip for High Heating Rate at evaporator inlet [-]

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        # State arrays for drum temperature and HTF inlet temperature history
        # used for heating-rate alarm calculations (replacing TRNSYS static arrays).
        # Size is 2*N_int; initialized in the first-step block when N_int is known.
        self._N_int: int = 1
        self._T_static: np.ndarray = np.zeros(2)  # [drum temps | HTF temps], length 2*N_int

    def calculate(self):

        # ------------------------------------------------------------------
        # Perform Any "After Convergence" Manipulations at the End of Each Timestep
        # ------------------------------------------------------------------
        if self.model.is_converged:

            N = 60.0 / (self.model.timestep * 3600.0)
            N_int = math.ceil(N)

            L_tank_new_val = self.L_tank_new.v
            P_tank_new_val = self.P_tank_new.v
            # Save new tank level, pressure and enthalpy as previous values for next timestep
            self.m_tank_prev.v = self.m_tank_new.v          # Mass in Tank
            self.L_tank_prev.v = L_tank_new_val             # Tank Level at beginning of timestep
            self.P_tank_prev.v = P_tank_new_val             # Tank Pressure at beginning of timestep
            self.h_tank_prev.v = self.h_tank_new.v          # Enthalpy of Tank at beginning of timestep

            # Set Initial Temperature arrays to calculate HR over time,
            # Drum Heating Rate (indices 0..N_int-1), HTF Inlet HR (indices N_int..2*N_int-1)
            if self.model.is_first_step:
                T_tank_val = self.T_tank_new.v
                HTF_Temp_in_val = self.HTF_Temp_in.v
                for i in range(1, N_int + 1):
                    self._T_static[i - 1] = T_tank_val           # drum temperature history
                    self._T_static[N_int + i - 1] = HTF_Temp_in_val  # HTF inlet temperature history

            # !!!! Find Alarms and Trips !!!!
            # !!!! Steam Drum Alarms !!!!
            # Low Level Alarm
            if L_tank_new_val > self.LowLevel_Alarm_cond.v:  # No Alarm
                Alarm = 0.0
                Trip = 0.0
            else:
                Alarm = 1.0
                if L_tank_new_val > self.LowLevel_Trip_cond.v:  # No Trip
                    Trip = 0.0
                else:
                    Trip = 1.0
            self.LowLevel_Alarm.v = Alarm  # Low Level Alarm
            self.LowLevel_Trip.v = Trip    # Low Level Trip

            # High Level Alarm
            if L_tank_new_val < self.HighLevel_Alarm_cond.v:
                Alarm = 0.0
                Trip = 0.0
            else:
                Alarm = 1.0
                if L_tank_new_val < self.HighLevel_Trip_cond.v:
                    Trip = 0.0
                else:
                    Trip = 1.0
            self.HighLevel_Alarm.v = Alarm  # High Level Alarm
            self.HighLevel_Trip.v = Trip    # High Level Trip

            # High Pressure Conditions
            if P_tank_new_val < self.HighPressure_Alarm_cond.v:
                Alarm = 0.0
                Trip = 0.0
            else:
                Alarm = 1.0
                if P_tank_new_val < self.HighPressure_Trip_cond.v:
                    Trip = 0.0
                else:
                    Trip = 1.0
            self.HighPressure_Alarm.v = Alarm  # High Pressure Alarm
            self.HighPressure_Trip.v = Trip    # High Pressure Trip

            # High Temp Diff between water entering drum and saturation temp
            T_tank_val = self.T_tank_new.v
            # CONVERTED-NEEDS UNITS CHECK: removed /1000 Pa->kPa and /1000 J/kg->kJ/kg; eeslib uses SI
            T_in_fw = fp.temperature("water", P=self.P_in.v, H=self.h_in.v)  # eeslib call for T from P,H for water
            Delta_T = abs(T_in_fw - T_tank_val)
            if Delta_T < self.HighDeltaT_Alarm_cond.v:
                Alarm = 0.0
                Trip = 0.0
            else:
                Alarm = 1.0
                if Delta_T < self.HighDeltaT_Trip_cond.v:
                    Trip = 0.0
                else:
                    Trip = 1.0
            self.HighDeltaT_Alarm.v = Alarm
            self.HighDeltaT_Trip.v = Trip

            # High HR in Steam Drum
            T_arr = np.empty(N_int)
            HR = 0.0
            for i in range(1, N_int + 1):
                T_arr[i - 1] = self._T_static[i - 1]
                if i > 1:
                    HR = HR + (T_arr[i - 1] - T_arr[i - 2]) / (self.model.timestep * 3600.0)
            HR = HR + (T_tank_val - T_arr[N_int - 1]) / (self.model.timestep * 3600.0)
            HR = HR / N_int * 60  # (Divide by number of HR solved for and multiply by 60 to get K/min)
            if abs(HR) < self.HighHR_Alarm_cond.v:
                Alarm = 0.0
                Trip = 0.0
            else:
                Alarm = 1.0
                if abs(HR) < self.HighHR_Trip_cond.v:
                    Trip = 0.0
                else:
                    Trip = 1.0
            self.HR_drum.v = HR        # Heating Rate
            self.HighHR_Alarm.v = Alarm  # Alarm for High Heating rate in the steam drum
            self.HighHR_Trip.v = Trip    # Trip for High Heating Rate in the Steam Drum
            # Reset Temperature Static Array Variables for next endoftimestep
            for i in range(1, N_int + 1):
                if i != N_int:
                    self._T_static[i - 1] = T_arr[i]  # shift forwards
                else:
                    self._T_static[i - 1] = T_tank_val

            # !!!! Evaporator Alarms and Trips !!!!
            # High temp diff between T_HTF_in and T_sat
            Delta_T = self.HTF_Temp_in.v - T_tank_val
            if Delta_T < self.HighTempDiffEvap_Alarm_cond.v:  # No alarm
                Alarm = 0.0
                Trip = 0.0
            else:
                Alarm = 1.0
                if Delta_T < self.HighTempDiffEvap_Trip_cond.v:
                    Trip = 0.0
                else:
                    Trip = 1.0
            self.HighTempDiffEvap_Alarm.v = Alarm  # Alarm for High Temp Diff between evap HTF in and T_sat
            self.HighTempDiffEvap_Trip.v = Trip    # Trip for High Temp Diff between evap HTF in and T_sat

            # High HTF HR at Evaporator Inlet
            T_htf_arr = np.empty(N_int)
            HR = 0.0
            for i in range(1, N_int + 1):
                T_htf_arr[i - 1] = self._T_static[N_int + i - 1]
                if i > 1:
                    HR = HR + (T_htf_arr[i - 1] - T_htf_arr[i - 2]) / (self.model.timestep * 3600.0)
            HR = HR + (self.HTF_Temp_in.v - T_htf_arr[N_int - 1]) / (self.model.timestep * 3600.0)
            HR = HR / N_int * 60  # (Divide by number of HR solved for and multiply by 60 to get K/min)
            if abs(HR) < self.HighHTF_HR_Alarm_cond.v:
                Alarm = 0.0
                Trip = 0.0
            else:
                Alarm = 1.0
                if abs(HR) < self.HighHTF_HR_Trip_cond.v:
                    Trip = 0.0
                else:
                    Trip = 1.0
            self.HR_HTF.v = HR               # HTF heating rate
            self.HighHTF_HR_Alarm.v = Alarm  # Alarm for High Heating rate at evaporator inlet
            self.HighHTF_HR_Trip.v = Trip    # Trip for High Heating Rate at evaporator inlet
            # Reset Temperature Static Array Variables for next endoftimestep
            for i in range(1, N_int + 1):
                if i != N_int:
                    self._T_static[N_int + i - 1] = T_htf_arr[i]  # shift forwards
                else:
                    self._T_static[N_int + i - 1] = self.HTF_Temp_in.v

            return

        # ------------------------------------------------------------------
        # Do All of the First Timestep Manipulations Here
        # There Are No Iterations at the Initial Time
        # ------------------------------------------------------------------
        if self.model.is_first_step:

            N = 60.0 / (self.model.timestep * 3600.0)
            N_int = math.ceil(N)
            self._N_int = N_int
            # Allocate static temperature history arrays (drum + HTF inlet), size 2*N_int
            self._T_static = np.zeros(2 * N_int)

            # CONVERTED-NEEDS UNITS CHECK: removed /1000 Pa->kPa and *1000 kJ/kg->J/kg; eeslib uses SI
            T_tank, h_out, rho_tank_g = (
                fp.t_sat("water", P=self.P_tank_init.v),   # eeslib call for T_sat from P,Q
                fp.enthalpy("water", P=self.P_tank_init.v, Q=1.0),       # eeslib call for h from P,Q
                fp.density("water", P=self.P_tank_init.v, Q=1.0),        # eeslib call for rho from P,Q
            )

            # Calculate current quality in the tank based on level
            # tank radius
            R_D = 1.0
            Area_liquid = (
                3.14 * self.D_tank.v ** 2.0 / 4.0
                * (0.5 - (math.asin(1.0 - 2.0 * self.L_tank_init.v / self.D_tank.v) / 3.14))
                - math.sqrt((self.D_tank.v / 2.0) ** 2.0 - (self.D_tank.v / 2.0 - self.L_tank_init.v) ** 2.0)
                * (self.D_tank.v / 2.0 - self.L_tank_init.v)
            )
            # total volume of liquid in the tank
            Vol_liquid = (
                Area_liquid * (self.Length_tank.v - self.D_tank.v / R_D)
                + 3.14 * (self.L_tank_init.v ** 2.0 * self.D_tank.v / 2.0 - self.L_tank_init.v ** 3.0 / 2.0) / R_D
            )
            # total volume of the tank
            Vol_tank = (
                3.14 * self.D_tank.v ** 2.0 / 4.0 * (self.Length_tank.v - self.D_tank.v)
                + 4.0 / 3.0 * 3.14 * (self.D_tank.v / 2.0) ** 3.0
            )

            # mass of liquid in the tank
            # density of saturated liquid water at this pressure
            # CONVERTED-NEEDS UNITS CHECK: removed /1000 Pa->kPa; eeslib uses SI
            rho_tank_f = fp.density("water", P=self.P_tank_init.v, Q=0.0)  # eeslib call for rho_f from P,Q
            m_tank_f = Vol_liquid * rho_tank_f  # total mass of liquid in the tank

            # mass of vapor in the tank
            Vol_vapor = Vol_tank - Vol_liquid
            m_tank_g = Vol_vapor * rho_tank_g

            # total mass in the tank
            m_tank_tot = m_tank_f + m_tank_g

            # quality of the tank based on initial conditions
            x_tank = m_tank_g / m_tank_tot  # starting quality in the tank

            # enthalpy of the tank based on initial conditions
            # CONVERTED-NEEDS UNITS CHECK: removed /1000 Pa->kPa and *1000 kJ/kg->J/kg; eeslib uses SI
            h_tank = fp.enthalpy("water", P=self.P_tank_init.v, Q=x_tank)  # eeslib call for h from P,Q

            # initial guess for steam flow leaving
            m_dot_superheat_init = 60.0
            Vol_dot_fw = m_dot_superheat_init / rho_tank_g

            # Volumetric flow rate for HTF
            rho_htf = fp.density(self.Fluid_ID.v, T=self.HTF_Temp_in.v, P=self.HTF_P_in.v)
            Vol_dot_HTF = self.HTF_mass_in.v / rho_htf

            # Set the Initial Values of the Outputs (#, Value)
            self.m_dot_steam.v = m_dot_superheat_init   # Mass of steam leaving the steam drum
            self.Vol_dot_fw.v = Vol_dot_fw              # Volumetric flow rate of steam leaving the steam drum
            self.P_steam_out.v = self.P_tank_init.v     # Pressure of steam leaving the steam drum
            self.T_steam_out.v = T_tank                 # Temperature of steam leaving the steam drum
            self.h_steam_out.v = h_out                  # Enthalpy of steam leaving the steam drum
            self.HTF_mass_out.v = self.HTF_mass_in.v    # HTF mass flow leaving the evaporator
            self.Vol_dot_HTF.v = Vol_dot_HTF            # HTF volumetric flow rate leaving the evaporator
            self.HTF_Temp_out.v = self.HTF_Temp_in.v    # HTF temperature leaving the evaporator
            self.HTF_P_out.v = self.HTF_P_in.v          # HTF pressure leaving the evaporator
            self.m_tank_new.v = m_tank_tot              # Total mass in the tank at the end of the timestep
            self.L_tank_new.v = self.L_tank_init.v      # Current Level in the tank at the end of the timestep
            self.P_tank_new.v = self.P_tank_init.v      # Current Pressure of tank at end of timestep
            self.T_tank_new.v = T_tank                  # Current Tank Temperature at end of timestep
            self.h_tank_new.v = h_tank                  # Current Tank enthalpy at end of timestep
            self.m_tank_prev.v = m_tank_tot             # Total mass in the tank at the beginning of the timestep
            self.L_tank_prev.v = self.L_tank_init.v     # Tank Level at beginning of timestep
            self.P_tank_prev.v = self.P_tank_init.v     # Tank Pressure at beginning of timestep
            self.h_tank_prev.v = h_tank                 # Enthalpy of Tank at beginning of timestep

            return

        # ------------------------------------------------------------------
        # Read the Inputs and Parameters (main iteration body)
        # ------------------------------------------------------------------

        # Evaporator calculations
        HTF_mass_in_val = self.HTF_mass_in.v
        if HTF_mass_in_val == 0.0:
            HTF_mass_in_val = 0.00001

        # CONVERTED-NEEDS UNITS CHECK: removed /1000 Pa->kPa; eeslib uses SI
        T_tank_prev, h_evap_out = (
            fp.temperature("water", P=self.P_tank_prev.v, Q=1.0),  # eeslib T from P,Q
            fp.enthalpy("water", P=self.P_tank_prev.v, Q=1.0),     # eeslib h from P,Q
        )
        # Surface area of the feedwater side of the heat exchanger
        A_s = (
            3.14
            * (self.Tube_OD.v - 2.0 * self.Tube_th.v)
            * self.Length_HX.v
            * self.No_tube_passes.v
            * self.No_tubes.v
        )
        UA_rated = self.Rated_heat_transfer.v * A_s
        # Off-design UA
        UA_OD = UA_rated * ((HTF_mass_in_val / self.Rated_htf_flow.v) ** self.Rated_exp.v)
        # specific heat of HTF Fluid entering
        cp_htf_max = fp.specheat(self.Fluid_ID.v, T=self.HTF_Temp_in.v, P=self.HTF_P_in.v)
        # lowest specific heat of HTF Fluid leaving
        cp_htf_min = fp.specheat(self.Fluid_ID.v, T=T_tank_prev, P=self.HTF_P_in.v)
        cp_htf_ave = (cp_htf_max + cp_htf_min) / 2.0
        NTU_OD = UA_OD / (HTF_mass_in_val * cp_htf_ave)
        # effectiveness of a one pass shell and tube heat exchanger when CR is equal to 0
        # (can never be equal to one otherwise will result in division by zero)
        eta_1pass = min(
            2.0 * (1.0 + (1.0 + math.exp(-NTU_OD)) / (1.0 - math.exp(-NTU_OD))) ** (-1.0),
            0.99999999,
        )
        CR = 0.0
        eta_od = (
            ((1.0 - eta_1pass * CR) / (1.0 - eta_1pass)) ** self.No_shell_passes.v - 1.0
        ) / (
            ((1.0 - eta_1pass * CR) / (1.0 - eta_1pass)) ** self.No_shell_passes.v - CR
        )

        # CONVERTED-NEEDS UNITS CHECK: removed /1000 Pa->kPa and *1000 kJ/kg->J/kg; eeslib uses SI
        h_evap_out = fp.enthalpy("water", P=self.P_tank_prev.v, Q=1.0)  # eeslib h from P,Q
        Q_dot_max = HTF_mass_in_val * cp_htf_ave * (self.HTF_Temp_in.v - T_tank_prev)
        Q_dot_actual = eta_od * Q_dot_max

        HTF_mass_out = HTF_mass_in_val
        # Energy balance on HTF through evaporator
        HTF_temp_out = (
            HTF_mass_in_val * cp_htf_ave * self.HTF_Temp_in.v - Q_dot_actual
        ) / (HTF_mass_in_val * cp_htf_ave)
        HTF_P_out = self.HTF_P_in.v  # No pressure drop across HX
        # CONVERTED-NEEDS UNITS CHECK: removed /1000 Pa->kPa and *1000 kJ/kg->J/kg; eeslib uses SI
        h_evap_in = fp.enthalpy("water", P=self.P_tank_prev.v, Q=0.0)  # eeslib h from P,Q
        # Energy Balance on Water through evaporator
        m_dot_evap_val = Q_dot_actual / (h_evap_out - h_evap_in)

        # FOR VERIFICATION ONLY - BLOWDOWN IS NOT IN MODEL
        # m_dot_blowdown = min(m_dot_in * P_tank_prev/10000000.0 * 0.07, m_dot_in * 0.07)
        m_dot_blowdown = 0.0

        # Tank Calculations
        ts = self.model.timestep * 3600.0  # converting timestep from hr to s
        dh = 1000.0
        dP = 1000.0
        # Use tank enthalpy at the beginning of the timestep
        h_tank_prev_val = self.h_tank_prev.v
        P_tank_prev_val = self.P_tank_prev.v
        m_tank_prev_val = self.m_tank_prev.v
        L_tank_prev_val = self.L_tank_prev.v
        Vol_tank = (
            3.14 * self.D_tank.v ** 2.0 / 4.0 * (self.Length_tank.v - self.D_tank.v)
            + 4.0 / 3.0 * 3.14 * (self.D_tank.v / 2.0) ** 3.0
        )

        # aa calculations
        drhodhcp_a = drhodhcp(P_tank=P_tank_prev_val, h_tank=h_tank_prev_val, dh=dh)
        drhodpch_a = drhodpch(P_tank=P_tank_prev_val, h_tank=h_tank_prev_val, dP=dP)
        dudhcp_a = dudhcp(P_tank=P_tank_prev_val, h_tank=h_tank_prev_val, dh=dh)
        dudpch_a = dudpch(P_tank=P_tank_prev_val, h_tank=h_tank_prev_val, dP=dP)
        # CONVERTED-NEEDS UNITS CHECK: removed /1000 Pa->kPa, /1000 J/kg->kJ/kg, and *1000 kJ/kg->J/kg; eeslib uses SI
        rho_tank = fp.density("water", P=P_tank_prev_val, H=h_tank_prev_val)  # eeslib rho from P,H
        u_tank = fp.internalenergy("water", P=P_tank_prev_val, H=h_tank_prev_val)  # eeslib u from P,H
        # CONVERTED-NEEDS UNITS CHECK: removed /1000 Pa->kPa and *1000 kJ/kg->J/kg; eeslib uses SI
        h_sat_f = fp.enthalpy("water", P=P_tank_prev_val, Q=0.0)  # eeslib h_f from P,Q
        h_sat_g = fp.enthalpy("water", P=P_tank_prev_val, Q=1.0)  # eeslib h_g from P,Q
        h_sat_fg = h_sat_g - h_sat_f
        h_out = h_sat_g

        # get around divide by 0 error
        denominator = drhodhcp_a * dudpch_a - drhodpch_a * dudhcp_a
        denominator = math.copysign(max(abs(denominator), 0.00001), denominator)
        drhodhcp_a = math.copysign(max(abs(drhodhcp_a), 0.00005), drhodhcp_a)

        tracker_1 = drhodhcp_a
        tracker_2 = drhodpch_a
        tracker_3 = dudhcp_a
        tracker_4 = dudpch_a
        tracker_5 = denominator

        dpdt_aa = (
            (
                (u_tank - h_sat_g) * self.m_dot_superheat.v
                + (u_tank - h_sat_f) * m_dot_blowdown
                + (self.h_in.v - u_tank) * self.m_dot_in.v
                + Q_dot_actual
            ) * drhodhcp_a
            + dudhcp_a * rho_tank * (m_dot_blowdown + self.m_dot_superheat.v - self.m_dot_in.v)
        ) / (Vol_tank * rho_tank * denominator)
        dhdt_aa = (
            (self.m_dot_in.v - self.m_dot_superheat.v - m_dot_blowdown) / Vol_tank
            - drhodpch_a * dpdt_aa
        ) / drhodhcp_a

        P_aa = P_tank_prev_val + dpdt_aa * ts / 2.0
        h_aa = h_tank_prev_val + dhdt_aa * ts / 2.0

        # bb calculations
        drhodhcp_a = drhodhcp(P_tank=P_aa, h_tank=h_aa, dh=dh)
        drhodpch_a = drhodpch(P_tank=P_aa, h_tank=h_aa, dP=dP)
        dudhcp_a = dudhcp(P_tank=P_aa, h_tank=h_aa, dh=dh)
        dudpch_a = dudpch(P_tank=P_aa, h_tank=h_aa, dP=dP)
        # CONVERTED-NEEDS UNITS CHECK: removed /1000 Pa->kPa, /1000 J/kg->kJ/kg, and *1000 kJ/kg->J/kg; eeslib uses SI
        rho_tank = fp.density("water", P=P_aa, H=h_aa)         # eeslib rho from P,H
        u_tank = fp.internalenergy("water", P=P_aa, H=h_aa)    # eeslib u from P,H
        # CONVERTED-NEEDS UNITS CHECK: removed /1000 Pa->kPa and *1000 kJ/kg->J/kg; eeslib uses SI
        h_sat_f = fp.enthalpy("water", P=P_aa, Q=0.0)  # eeslib h_f from P,Q
        h_sat_g = fp.enthalpy("water", P=P_aa, Q=1.0)  # eeslib h_g from P,Q
        h_sat_fg = h_sat_g - h_sat_f
        h_out = h_sat_g
        dpdt_bb = (
            (
                (u_tank - h_sat_g) * self.m_dot_superheat.v
                + (u_tank - h_sat_f) * m_dot_blowdown
                + (self.h_in.v - u_tank) * self.m_dot_in.v
                + Q_dot_actual
            ) * drhodhcp_a
            + dudhcp_a * rho_tank * (m_dot_blowdown + self.m_dot_superheat.v - self.m_dot_in.v)
        ) / (Vol_tank * rho_tank * denominator)
        dhdt_bb = (
            (self.m_dot_in.v - self.m_dot_superheat.v - m_dot_blowdown) / Vol_tank
            - drhodpch_a * dpdt_bb
        ) / drhodhcp_a

        P_bb = P_tank_prev_val + dpdt_bb * ts / 2.0
        h_bb = h_tank_prev_val + dhdt_bb * ts / 2.0

        denominator = drhodhcp_a * dudpch_a - drhodpch_a * dudhcp_a
        denominator = math.copysign(max(abs(denominator), 0.00001), denominator)
        drhodhcp_a = math.copysign(max(abs(drhodhcp_a), 0.00005), drhodhcp_a)

        tracker_6 = drhodhcp_a
        tracker_7 = drhodpch_a
        tracker_8 = dudhcp_a
        tracker_9 = dudpch_a
        tracker_10 = denominator

        # cc calculations
        drhodhcp_a = drhodhcp(P_tank=P_bb, h_tank=h_bb, dh=dh)
        drhodpch_a = drhodpch(P_tank=P_bb, h_tank=h_bb, dP=dP)
        dudhcp_a = dudhcp(P_tank=P_bb, h_tank=h_bb, dh=dh)
        dudpch_a = dudpch(P_tank=P_bb, h_tank=h_bb, dP=dP)
        # CONVERTED-NEEDS UNITS CHECK: removed /1000 Pa->kPa, /1000 J/kg->kJ/kg, and *1000 kJ/kg->J/kg; eeslib uses SI
        rho_tank = fp.density("water", P=P_bb, H=h_bb)         # eeslib rho from P,H
        u_tank = fp.internalenergy("water", P=P_bb, H=h_bb)    # eeslib u from P,H
        # CONVERTED-NEEDS UNITS CHECK: removed /1000 Pa->kPa and *1000 kJ/kg->J/kg; eeslib uses SI
        h_sat_f = fp.enthalpy("water", P=P_bb, Q=0.0)  # eeslib h_f from P,Q
        h_sat_g = fp.enthalpy("water", P=P_bb, Q=1.0)  # eeslib h_g from P,Q
        h_sat_fg = h_sat_g - h_sat_f
        h_out = h_sat_g
        dpdt_cc = (
            (
                (u_tank - h_sat_g) * self.m_dot_superheat.v
                + (u_tank - h_sat_f) * m_dot_blowdown
                + (self.h_in.v - u_tank) * self.m_dot_in.v
                + Q_dot_actual
            ) * drhodhcp_a
            + dudhcp_a * rho_tank * (m_dot_blowdown + self.m_dot_superheat.v - self.m_dot_in.v)
        ) / (Vol_tank * rho_tank * denominator)
        dhdt_cc = (
            (self.m_dot_in.v - self.m_dot_superheat.v - m_dot_blowdown) / Vol_tank
            - drhodpch_a * dpdt_cc
        ) / drhodhcp_a

        denominator = drhodhcp_a * dudpch_a - drhodpch_a * dudhcp_a
        denominator = math.copysign(max(abs(denominator), 0.00001), denominator)
        drhodhcp_a = math.copysign(max(abs(drhodhcp_a), 0.00005), drhodhcp_a)

        tracker_11 = drhodhcp_a
        tracker_12 = drhodpch_a
        tracker_13 = dudhcp_a
        tracker_14 = dudpch_a
        tracker_15 = denominator

        P_cc = P_tank_prev_val + dpdt_cc * ts
        h_cc = h_tank_prev_val + dhdt_cc * ts

        # dd calculations
        drhodhcp_a = drhodhcp(P_tank=P_cc, h_tank=h_cc, dh=dh)
        drhodpch_a = drhodpch(P_tank=P_cc, h_tank=h_cc, dP=dP)
        dudhcp_a = dudhcp(P_tank=P_cc, h_tank=h_cc, dh=dh)
        dudpch_a = dudpch(P_tank=P_cc, h_tank=h_cc, dP=dP)
        # CONVERTED-NEEDS UNITS CHECK: removed /1000 Pa->kPa, /1000 J/kg->kJ/kg, and *1000 kJ/kg->J/kg; eeslib uses SI
        rho_tank = fp.density("water", P=P_cc, H=h_cc)         # eeslib rho from P,H
        u_tank = fp.internalenergy("water", P=P_cc, H=h_cc)    # eeslib u from P,H
        # CONVERTED-NEEDS UNITS CHECK: removed /1000 Pa->kPa and *1000 kJ/kg->J/kg; eeslib uses SI
        h_sat_f = fp.enthalpy("water", P=P_cc, Q=0.0)  # eeslib h_f from P,Q
        h_sat_g = fp.enthalpy("water", P=P_cc, Q=1.0)  # eeslib h_g from P,Q
        h_sat_fg = h_sat_g - h_sat_f
        h_out = h_sat_g

        denominator = drhodhcp_a * dudpch_a - drhodpch_a * dudhcp_a
        denominator = math.copysign(max(abs(denominator), 0.00001), denominator)
        drhodhcp_a = math.copysign(max(abs(drhodhcp_a), 0.00005), drhodhcp_a)

        tracker_16 = drhodhcp_a
        tracker_17 = drhodpch_a
        tracker_18 = dudhcp_a
        tracker_19 = dudpch_a
        tracker_20 = denominator

        dpdt_dd = (
            (
                (u_tank - h_sat_g) * self.m_dot_superheat.v
                + (u_tank - h_sat_f) * m_dot_blowdown
                + (self.h_in.v - u_tank) * self.m_dot_in.v
                + Q_dot_actual
            ) * drhodhcp_a
            + dudhcp_a * rho_tank * (m_dot_blowdown + self.m_dot_superheat.v - self.m_dot_in.v)
        ) / (Vol_tank * rho_tank * denominator)
        dhdt_dd = (
            (self.m_dot_in.v - self.m_dot_superheat.v - m_dot_blowdown) / Vol_tank
            - drhodpch_a * dpdt_dd
        ) / drhodhcp_a

        # End of timestep Pressure and Enthalpy (4th-order Runge-Kutta)
        P_tank_new_val = P_tank_prev_val + (dpdt_aa + 2.0 * dpdt_bb + 2.0 * dpdt_cc + dpdt_dd) * ts / 6.0
        h_tank_new_val = h_tank_prev_val + (dhdt_aa + 2.0 * dhdt_bb + 2.0 * dhdt_cc + dhdt_dd) * ts / 6.0
        rho_tank_new = (m_tank_prev_val + (self.m_dot_in.v - self.m_dot_superheat.v) * ts) / Vol_tank
        # CONVERTED-NEEDS UNITS CHECK: removed /1000 Pa->kPa and /1000 J/kg->kJ/kg; eeslib uses SI
        x_tank_new = fp.quality("water", P=P_tank_new_val, H=h_tank_new_val)   # eeslib quality from P,H
        rho_tank_new = fp.density("water", P=P_tank_new_val, H=h_tank_new_val) # eeslib rho from P,H
        T_tank_new_val = fp.temperature("water", P=P_tank_new_val, H=h_tank_new_val)  # eeslib T from P,H

        # Solve for new level in the tank
        # total mass of water in the tank
        m_tank_new_val = m_tank_prev_val + (self.m_dot_in.v - self.m_dot_superheat.v - m_dot_blowdown) * ts
        m_tank_g_new = m_tank_new_val * x_tank_new     # total mass of vapor in the tank
        m_tank_f_new = m_tank_new_val - m_tank_g_new   # total mass of liquid in the tank
        # CONVERTED-NEEDS UNITS CHECK: removed /1000 Pa->kPa; eeslib uses SI
        rho_f_new = fp.density("water", P=P_tank_new_val, Q=0.0)  # eeslib rho_f from P,Q
        Vol_liquid = m_tank_f_new / rho_f_new
        level_tol = 0.000001  # [m]
        L_tank_new_val = tank_level(Vol_liquid, self.D_tank.v, self.Length_tank.v, L_tank_prev_val, level_tol)

        # enthalpy of steam leaving the steam drum
        # CONVERTED-NEEDS UNITS CHECK: removed /1000 Pa->kPa and *1000 kJ/kg->J/kg; eeslib uses SI
        h_out = fp.enthalpy("water", P=P_tank_new_val, Q=1.0)    # eeslib h_g from P,Q
        rho_g_new = fp.density("water", P=P_tank_new_val, Q=1.0) # eeslib rho_g from P,Q

        # calculating volumetric flow rates
        Vol_dot_fw = self.m_dot_superheat.v / rho_g_new
        rho_htf = fp.density(self.Fluid_ID.v, T=HTF_temp_out, P=HTF_P_out)
        Vol_dot_htf = HTF_mass_out / rho_htf

        # ------------------------------------------------------------------
        # Set the Outputs from this Model (#, Value)
        # ------------------------------------------------------------------
        self.m_dot_steam.v = self.m_dot_superheat.v  # mass flow rate of steam leaving the steam drum
        self.Vol_dot_fw.v = Vol_dot_fw               # Volumetric flow rate of steam leaving the steam drum
        # Pressure of superheat leaving the steam drum, constant during timestep iterations,
        # only changes at beginning of each timestep
        self.P_steam_out.v = P_tank_prev_val
        # Temperature of superheat leaving the steam drum
        self.T_steam_out.v = T_tank_prev
        # enthalpy of superheat leaving the steam drum, constant during timestep iterations,
        # only changes at the beginning of each timestep
        self.h_steam_out.v = h_out
        self.HTF_mass_out.v = HTF_mass_out           # mass flow rate of htf leaving the evaporator
        self.Vol_dot_HTF.v = Vol_dot_htf             # Volumetric flow rate of htf leaving the evaporator
        self.HTF_Temp_out.v = HTF_temp_out           # Temperature of htf leaving the evaporator
        self.HTF_P_out.v = HTF_P_out                 # Pressure of HTF leaving the evaporator
        self.m_tank_new.v = m_tank_new_val           # new mass in tank at the end of the timestep
        self.L_tank_new.v = L_tank_new_val           # New Tank Level at end of timestep
        self.P_tank_new.v = P_tank_new_val           # New Tank Pressure at end of timestep
        self.T_tank_new.v = T_tank_new_val           # New Tank Temperature at end of timestep
        self.h_tank_new.v = h_tank_new_val           # New Tank enthalpy at end of timestep
        self.eta_OD.v = eta_od                       # off-design effectiveness of the evaporator
        self.m_dot_evap.v = m_dot_evap_val           # amount of water that was evaporated in the steam drum during this timestep
