"""Turbines and bypass network component model (Type 6028)."""

import math

from eeslib import fluid_properties as fp
from esclab.simulate import Component
from esclab.components.flownetwork.esol6015_helpers import (
    drhodhcp, drhodpch, dudhcp, dudpch,
    f_cp_water,
    PB_CV_data,
    VP_new,
    StodolaStage,
    h_lpt_stage,
    eta_SCC_hpt, eta_SCC_HPT,
    eta_SCC_lpt, eta_SCC_LPT,
    convection_dynamicpipe,
    valve_massflow,
    specheat,
    density,
)


class TurbinesAndBypassNetwork(Component):
    """
    Object: ESOL6028-Turbines&BypassNetwork
    Simulation Studio Model: ESOL6028-Turbines&BypassNetwork

    Models the combined turbines and bypass network for a steam power cycle,
    including high-pressure and low-pressure turbines, steam separator, reheater,
    main piping dynamic pressure and enthalpy calculations, and alarm/trip logic.
    """

    # -------------------------------------------------------------------------
    # Parameters
    # -------------------------------------------------------------------------
    P_HPmain_ini = Component.Parameter()       # Initial Pressure in the main piping line to the HP turbine [Pa]
    T_HPmain_ini = Component.Parameter()       # Initial Temperature in the main piping line to the HP turbine [K]
    Length_HPmain = Component.Parameter()      # Length of the main piping line to the HP turbine [m]
    D_HPmain = Component.Parameter()           # Diameter of the main piping line to the HP turbine [m]
    P_LPmain_ini = Component.Parameter()       # Initial Pressure of the main piping line to the LP turbine [Pa]
    T_LPmain_ini = Component.Parameter()       # Initial Temperature of the main piping line to the LP turbine [K]
    Length_LPmain = Component.Parameter()      # Length of the main piping line to the LP turbine [m]
    D_LPmain = Component.Parameter()           # Diameter of the main piping line to the LP turbine [m]
    P_Aux_ini = Component.Parameter()          # Initial Pressure of the bypass auxillary line [Pa]
    T_Aux_ini = Component.Parameter()          # Initial Temperature of the bypass auxillary line [K]
    Length_AuxLine = Component.Parameter()     # Length of the bypass auxillary line [m]
    D_AuxLine = Component.Parameter()          # Diameter of the bypass auxillary line [m]
    mc_HPmain_pipe = Component.Parameter()     # Mass and specific heat of metal HP piping line [J/K]
    mc_LPmain_pipe = Component.Parameter()     # Mass and specific heat of metal LP piping line [J/K]
    mc_Aux_pipe = Component.Parameter()        # Mass and specific heat of metal Aux piping line [J/K]
    HP_Bypass_d = Component.Parameter()
    HP_Bypass_vs = Component.Parameter()
    HP_Bypass_vt = Component.Parameter()
    HP_Aux_d = Component.Parameter()
    HP_Aux_vs = Component.Parameter()
    HP_Aux_vt = Component.Parameter()
    HP_warmup_d = Component.Parameter()
    HP_warmup_vs = Component.Parameter()
    HP_warmup_vt = Component.Parameter()
    HP_drain_d = Component.Parameter()
    HP_drain_vs = Component.Parameter()
    HP_drain_vt = Component.Parameter()
    Aux_DA_d = Component.Parameter()
    Aux_DA_vs = Component.Parameter()
    Aux_DA_vt = Component.Parameter()
    LP_Bypass_d = Component.Parameter()
    LP_Bypass_vs = Component.Parameter()
    LP_Bypass_vt = Component.Parameter()
    LP_Aux_d = Component.Parameter()
    LP_Aux_vs = Component.Parameter()
    LP_Aux_vt = Component.Parameter()
    LP_warmup_d = Component.Parameter()
    LP_warmup_vs = Component.Parameter()
    LP_warmup_vt = Component.Parameter()
    LP_drain_d = Component.Parameter()
    LP_drain_vs = Component.Parameter()
    LP_drain_vt = Component.Parameter()
    m_dot_turbine_seals = Component.Parameter()
    P_ts_req = Component.Parameter()
    HPT_no_GS = Component.Parameter()          # Number of governing stages in the HPT (must be 1 or 2)
    HPT_parallel_sects = Component.Parameter()  # Number of parallel sections in the HPT
    HPT_CV_NUMBER = Component.Parameter()       # Number of control valves before the HPT
    HPT_CV_D = Component.Parameter()            # Diameter of the control valves before the HPT [m]
    HPT_CV_VS = Component.Parameter()           # Valve speed [deg/s]
    HPT_CV_VPD = Component.Parameter()          # Design Control Valve Position for HPT
    GS_diameter = Component.Parameter()         # HPT governing stage pitch diameter [m]
    m_dot_HPT_d = Component.Parameter()         # Design mass flow rate into HPT [kg/s]
    P_HPT_in_d = Component.Parameter()          # Design pressure into HPT [Pa]
    T_HPT_in_d = Component.Parameter()          # Design temperature into HPT [K]
    P_HPT1_d = Component.Parameter()            # Design Pressure of steam extraction leaving HPT stage 1 [Pa]
    P_HPT_exh_d = Component.Parameter()         # Design exhaust pressure leaving HPT [Pa]
    m_dot_HPT1_d = Component.Parameter()        # Design mass flow rate of steam extraction 1 leaving HPT [kg/s]
    m_dot_HPT2_d = Component.Parameter()        # Design mass flow rate of steam extraction 2 leaving HPT [kg/s]
    HX_UA_d = Component.Parameter()             # Design UA for reheater [W/K]
    m_dot_fw_HX_d = Component.Parameter()       # design mass flow rate of feedwater through reheater [kg/s]
    m_dot_HTF_HX_d = Component.Parameter()      # design mass flow rate of HTF through reheater [kg/s]
    HX_exp = Component.Parameter()              # heat transfer off-design exponent
    HX_no_shell = Component.Parameter()         # Number of shells in reheater
    HX_length = Component.Parameter()           # Length of reheater [m]
    HX_tube_OD = Component.Parameter()          # Reheater tubes outer diameter [m]
    HX_tube_th = Component.Parameter()          # Reheater tube thickness [m]
    HX_No_tubes = Component.Parameter()         # Number of tubes in reheater
    Fluid_ID = Component.Parameter()            # Fluid ID for HTF
    LPT_parallel_sects = Component.Parameter()  # Number of parallel sections in LPT
    LPT_EXP_A0 = Component.Parameter()          # LPT expansion coefficient A
    LPT_EXP_A1 = Component.Parameter()          # LPT expansion coefficient B
    LPT_EXP_A2 = Component.Parameter()          # LPT expansion coefficient C
    m_dot_LPT_d = Component.Parameter()         # Design mass flow rate through LPT [kg/s]
    P_LPT_in_d = Component.Parameter()          # Design pressure entering LPT [Pa]
    T_LPT_in_d = Component.Parameter()          # Design temperature entering LPT [K]
    P_LPT1_d = Component.Parameter()            # Design pressure of steam extraction 1 from LPT [Pa]
    P_LPT2_d = Component.Parameter()            # Design pressure of steam extraction 2 from LPT [Pa]
    P_LPT3_d = Component.Parameter()            # Design pressure of steam extraction 3 from LPT [Pa]
    P_LPT4_d = Component.Parameter()            # Design pressure of steam extraction 4 from LPT [Pa]
    m_dot_LPT1_d = Component.Parameter()        # Design mass flow rate of steam extraction 1 from LPT [kg/s]
    m_dot_LPT2_d = Component.Parameter()        # Design mass flow rate of steam extraction 2 from LPT [kg/s]
    m_dot_LPT3_d = Component.Parameter()        # Design mass flow rate of steam extraction 3 from LPT [kg/s]
    m_dot_LPT4_d = Component.Parameter()        # Design mass flow rate of steam extraction 4 from LPT [kg/s]
    P_cond_d = Component.Parameter()            # Design exhaust pressure of LPT (condenser) [Pa]
    HPT_SH_Alarm_FL = Component.Parameter()     # Superheat Alarm for HPT at Full Load [Delta C]
    HPT_SH_Trip_FL = Component.Parameter()      # Superheat Trip for HPT at Full Load [Delta C]
    HPT_SH_Alarm_PL = Component.Parameter()     # Superheat Alarm for HPT at Partial Load [Delta C]
    # NOTE: Fortran code incorrectly assigned HTP_SH_ALARM_PL twice, leaving HPT_SH_Trip_PL unassigned. This is fixed in the conversion
    HPT_SH_Trip_PL = Component.Parameter()      # Superheat Trip for HPT at Partial Load [Delta C]
    Partial_Load = Component.Parameter()        # Partial Load Defining Factor [W]
    HPT_HighTemp_Alarm = Component.Parameter()  # High Temperature entering HPT Alarm [K]
    HPT_HighTemp_TimedTrip = Component.Parameter()  # High Temperature entering HPT Timed Trip [K]
    HPT_HighTemp_Trip = Component.Parameter()   # High Temperature entering HPT Trip [K]
    HPT_ExhPres_Alarm = Component.Parameter()   # High Exhaust Pressure Leaving HPT Alarm [Pa]
    HPT_ExhPres_Trip = Component.Parameter()    # High Exhaust Pressure Leaving HPT Trip [Pa]
    LPT_SH_Alarm = Component.Parameter()        # Superheat Alarm for LPT [Delta C]
    LPT_SH_Trip = Component.Parameter()         # Superheat Trip for LPT [Delta C]
    LPT_HighTemp_Alarm = Component.Parameter()  # High Temperature entering LPT Alarm [K]
    LPT_HighTemp_TimedTrip = Component.Parameter()  # High Temperature entering LPT Timed Trip [K]
    LPT_HighTemp_Trip = Component.Parameter()   # High Temperature entering LPT Trip [K]
    T_HPpipe_ini = Component.Parameter()        # Initial temperature of HP pipe metal [K]
    T_LPpipe_ini = Component.Parameter()        # Initial temperature of LP pipe metal [K]
    T_Auxpipe_ini = Component.Parameter()       # Initial temperature of Aux pipe metal [K]
    HighHTF_Temp_Alarm = Component.Parameter()  # High HTF Temperature entering reheater Alarm [K]
    HighHTF_Temp_Trip = Component.Parameter()   # High HTF Temperature entering reheater Trip [K]
    LowHTF_TempOut_Alarm = Component.Parameter()  # Low HTF Temperature leaving reheater Alarm [K]
    LowHTF_TempOut_Trip = Component.Parameter()   # Low HTF Temperature leaving reheater Trip [K]
    HighHTF_Flow_Alarm = Component.Parameter()  # High HTF Flow entering reheater Alarm [kg/s]
    HighHTF_Flow_Trip = Component.Parameter()   # High HTF Flow entering reheater Trip [kg/s]
    HighHTF_Pres_Alarm = Component.Parameter()  # High HTF Pressure entering reheaters Alarm [Pa]
    HighHTF_Pres_Trip = Component.Parameter()   # High HTF Pressure entering reheaters Trip [Pa]
    HighDT_HX_Alarm = Component.Parameter()     # High differential temperature between HTF inlet and steam inlet Alarm [Delta K]
    HighDT_HX_Trip = Component.Parameter()      # High differential temperature between HTF inlet and steam inlet Trip [Delta K]
    HighHTF_HR_Alarm = Component.Parameter()    # High HTF heating rate alarm [K/min]
    HighHTF_HR_Trip = Component.Parameter()     # High HTF heating rate trip [K/min]
    HighSteam_HR_Alarm = Component.Parameter()  # High Steam heating rate alarm [K/min]
    HighSteam_HR_Trip = Component.Parameter()   # High Steam heating rate trip [K/min]

    # -------------------------------------------------------------------------
    # Inputs
    # -------------------------------------------------------------------------
    Turbine_ON = Component.Input()       # Turbine is on signal (if on == 1)
    ifAuto_HPBypass = Component.Input()  # Signal for HPBypass Valve if acting in Auto mode or manual mode
    ifAuto_LPBypass = Component.Input()  # Signal for LPBypass Valve if acting in Auto mode
    ifAuto_HPAux = Component.Input()     # Signal for HPAux Valve if acting in auto or manual
    ifAuto_LPAux = Component.Input()     # Signal for LPAux Valve if acting in auto or manual
    ifAuto_DAaux = Component.Input()     # Signal for AUXDA Valve if acting in auto or manual
    ifAuto_HPTCV = Component.Input()     # Signal for HPT control valve if acting in auto or manual
    P_DA = Component.Input()             # Current Pressure of the Deaerator for this timestep [Pa]
    HPT_CV_VPi = Component.Input()       # Valve Position for high pressure turbine control valves
    HP_Bypass_VPi = Component.Input()    # Valve Position for high pressure bypass line
    HP_Aux_VPi = Component.Input()       # Valve Position for valve between hp-main piping and the aux-piping
    HP_warmup_VPi = Component.Input()    # Valve Position for high pressure warm-up line
    HP_drain_VPi = Component.Input()     # Valve Position of drain valve connected to hp-piping
    Aux_DA_VPi = Component.Input()       # Valve Position of aux valve leading to DA
    LP_Bypass_VPi = Component.Input()    # Valve Position of low pressure bypass line
    LP_Aux_VPi = Component.Input()       # Valve Position for valve between lp main piping and the aux piping
    LP_warmup_VPi = Component.Input()    # Valve Position for low pressure warmup line
    LP_drain_VPi = Component.Input()     # Valve Position for low pressure drain line
    P_SGT_in = Component.Input()         # Pressure of steam into the hpt system [Pa]
    h_SGT_in = Component.Input()         # Enthalpy of steam into the hpt system [J/kg]
    m_dot_HTF = Component.Input()        # mass flow rate of htf into the reheater [kg/s]
    HTF_P_in = Component.Input()         # pressure of htf into the reheater [Pa]
    HTF_T_in = Component.Input()         # temperature of htf into the reheater [K]
    P_cond = Component.Input()           # condenser pressure [Pa]
    m_dot_HPT1 = Component.Input()       # steam extraction 1 leaving high pressure turbine [kg/s]
    m_dot_HPT2 = Component.Input()       # steam extraction 2 leaving high pressure turbine [kg/s]
    m_dot_LPT1 = Component.Input()       # steam extraction 1 leaving low pressure turbine [kg/s]
    m_dot_LPT2 = Component.Input()       # steam extraction 2 leaving low pressure turbine [kg/s]
    m_dot_LPT3 = Component.Input()       # steam extraction 3 leaving low pressure turbine [kg/s]
    m_dot_LPT4 = Component.Input()       # steam extraction 4 leaving low pressure turbine [kg/s]
    PID_HPbypass = Component.Input()     # PID control signal to HP bypass valve
    PID_LPbypass = Component.Input()     # PID control signal to LP Bypass valve
    PID_HPaux = Component.Input()        # PID control signal to HP Aux
    PID_LPaux = Component.Input()        # PID control signal to LP Aux
    PID_DAaux = Component.Input()        # PID control signal to DA-Aux valve
    PID_HPTCV = Component.Input()        # PID control signal to HPT_CV valve

    # -------------------------------------------------------------------------
    # Outputs
    # -------------------------------------------------------------------------
    Turbine_ON_out = Component.Output()           # HPT & LPT Turbine On
    HPT_CV_VPo = Component.Output()               # HPT Control Valve Position output
    HP_Bypass_VPo = Component.Output()            # High Pressure Bypass Valve Position output
    HP_Aux_VPo = Component.Output()
    HP_warmup_VPo = Component.Output()
    HP_drain_VPo = Component.Output()
    Aux_DA_VPo = Component.Output()
    LP_Bypass_VPo = Component.Output()
    LP_Aux_VPo = Component.Output()
    LP_warmup_VPo = Component.Output()
    LP_drain_VPo = Component.Output()
    m_dot_SGT = Component.Output()                # output to SD type (ESOL6003) [kg/s]
    m_dot_cond = Component.Output()               # output to condenser type (ESOL6007) [kg/s]
    h_cond = Component.Output()                   # output to condenser type (ESOL6007) [J/kg]
    m_dot_DA = Component.Output()                 # output to DA type (ESOL6011) [kg/s]
    h_DA = Component.Output()                     # output to DA type (ESOL6011) [J/kg]
    W_dot_tot = Component.Output()                # Total Turbine Power [W]
    m_dot_HTF_out = Component.Output()            # Reheater - HTF Mass Flow Out [kg/s]
    Vol_dot_HTF = Component.Output()              # Reheater - HTF Volumetric Flow out [m3/s]
    HTF_P_out = Component.Output()                # Reheater - HTF Pressure out [Pa]
    T_HTF_out = Component.Output()                # Reheater - HTF Temperature out [K]
    m_dot_HPT_in = Component.Output()             # Mass Flow rate entering HPT [kg/s]
    P_GS_out = Component.Output()                 # HPT - GS Outlet Pressure [Pa]
    h_GS_out = Component.Output()                 # HPT - GS Outlet Enthalpy [J/kg]
    m_dot_HPTS1 = Component.Output()              # HPT-Stage 1 Inlet Mass flow [kg/s]
    P_HPT1 = Component.Output()                   # HPT-Stage 1 Outlet Pressure [Pa]
    H_HPT1 = Component.Output()                   # HPT-Stage 1 Outlet Enthalpy [J/kg]
    m_dot_HPTS2 = Component.Output()              # HPT-Stage 2 Inlet Mass Flow [kg/s]
    P_HPT2 = Component.Output()                   # HPT-Stage 2 Outlet Pressure [Pa]
    H_HPT2 = Component.Output()                   # HPT-Stage 2 Outlet Enthalpy [J/kg]
    m_dot_HPT_exh = Component.Output()            # HPT-Exhaust Mass Flow [kg/s]
    P_HPT_exh = Component.Output()                # HPT-Exhaust Pressure [Pa]
    H_HPT_exh = Component.Output()                # HPT-Exhaust Enthalpy [J/kg]
    T_HPT_exh = Component.Output()                # HPT-Exhaust Temperature [K]
    m_dot_SS_drain = Component.Output()           # SS-Drain flow to Deaerator [kg/s]
    Vol_dot_SS_drain = Component.Output()         # SS-Drain Volumetric Flow to Deaerator [m3/s]
    P_SS_drain = Component.Output()               # SS-Drain Pressure to Deaerator [Pa]
    h_SS_drain = Component.Output()               # SS-Drain Enthalpy to Deaerator [J/kg]
    T_SS_drain = Component.Output()               # SS-Drain Temperature to Deaerator [K]
    m_dot_SS_steam = Component.Output()           # SS-Steam Mass Flow to Reheater [kg/s]
    Vol_dot_SS_steam = Component.Output()         # SS-Steam Volumetric Flow to Reheater [m3/s]
    P_SS_steam = Component.Output()               # SS-Steam Pressure to Reheater [Pa]
    h_SS_steam = Component.Output()               # SS-Steam Enthalpy to Reheater [J/kg]
    T_SS_steam = Component.Output()               # SS-Steam Temperature to Reheater [K]
    m_dot_RH_steam = Component.Output()           # Reheater - Steam Mass Flow Rate Leaving Reheater [kg/s]
    Vol_dot_RH_steam = Component.Output()         # Reheater - Steam Volumetric Flow Rate Leaving Reheater [m3/s]
    T_HX_out = Component.Output()                 # Reheater - Steam Temperature Leaving Reheater [K]
    P_RH_steam = Component.Output()               # Reheater - Steam Pressure Leaving Reheater [Pa]
    Q_dot_HX = Component.Output()                 # Reheater-Total Heat Transfer [W]
    Eta_OD = Component.Output()                   # Reheater-Effectiveness [-]
    m_dot_LPT_in = Component.Output()             # LPT-Stage 1 Mass Flow [kg/s]
    P_LPT1 = Component.Output()                   # LPT-Stage 1 Outlet Pressure [Pa]
    T_LPT1 = Component.Output()                   # LPT-Stage 1 Outlet Temperature [K]
    H_LPT1 = Component.Output()                   # LPT-Stage 1 Outlet Enthalpy [J/kg]
    m_dot_LPTS2 = Component.Output()              # LPT-Stage 2 Mass Flow [kg/s]
    P_LPT2 = Component.Output()                   # LPT-Stage 2 Outlet Pressure [Pa]
    T_LPT2 = Component.Output()                   # LPT-Stage 2 Outlet Temperature [K]
    H_LPT2 = Component.Output()                   # LPT-Stage 2 Outlet Enthalpy [J/kg]
    m_dot_LPTS3 = Component.Output()              # LPT-Stage 3 Mass Flow [kg/s]
    P_LPT3 = Component.Output()                   # LPT-Stage 3 Outlet Pressure [Pa]
    T_LPT3 = Component.Output()                   # LPT-Stage 3 Outlet Temperature [K]
    H_LPT3 = Component.Output()                   # LPT-Stage 3 Outlet Enthalpy [J/kg]
    m_dot_LPTS4 = Component.Output()              # LPT-Stage 4 Mass Flow [kg/s]
    T_LPT4 = Component.Output()                   # LPT-Stage 4 Outlet Temperature [K]
    P_LPT4 = Component.Output()                   # LPT-Stage 4 Outlet Pressure [Pa]
    H_LPT4 = Component.Output()                   # LPT-Stage 4 Outlet Enthalpy [J/kg]
    m_dot_LPT_exh = Component.Output()            # LPT Exhaust Mass Flow [kg/s]
    Vol_dot_LPT_exh = Component.Output()          # LPT Exhaust Volumetric Flow [m3/s]
    T_LPT_exh = Component.Output()                # LPT Exhaust Temperature [K]
    H_LPT_exh = Component.Output()                # LPT Exhaust Enthalpy [J/kg]
    P_HPmain = Component.Output()                 # HP main piping pressure [Pa]
    T_HPmain = Component.Output()                 # HP main piping temperature [K]
    x_HPmain = Component.Output()                 # HP main piping quality [-]
    T_HP_pipe = Component.Output()                # HP pipe metal temperature [K]
    P_LPmain = Component.Output()                 # LP main piping pressure [Pa]
    T_LPmain = Component.Output()                 # LP main piping temperature [K]
    x_LPmain = Component.Output()                 # LP main piping quality [-]
    T_LP_pipe = Component.Output()                # LP pipe metal temperature [K]
    P_AUX = Component.Output()                    # AUX piping pressure [Pa]
    T_AUX = Component.Output()                    # AUX piping temperature [K]
    x_AUX = Component.Output()                    # AUX piping quality [-]
    T_AUX_pipe = Component.Output()               # AUX pipe metal temperature [K]
    m_dot_HP_bypass = Component.Output()          # HP bypass mass flow [kg/s]
    m_dot_HP_AUX = Component.Output()             # HP aux mass flow [kg/s]
    m_dot_HP_drain = Component.Output()           # HP drain mass flow [kg/s]
    m_dot_HP_warmup = Component.Output()          # HP warmup mass flow [kg/s]
    m_dot_LP_AUX = Component.Output()             # LP aux mass flow [kg/s]
    m_dot_LP_bypass = Component.Output()          # LP bypass mass flow [kg/s]
    m_dot_LP_drain = Component.Output()           # LP drain mass flow [kg/s]
    m_dot_LP_warmup = Component.Output()          # LP warmup mass flow [kg/s]
    m_dot_AUX_DA = Component.Output()             # AUX-DA mass flow [kg/s]
    ff_HPmain = Component.Output()                # Friction factor HP main pipe [-]
    ff_LPmain = Component.Output()                # Friction factor LP main pipe [-]
    ff_AUX = Component.Output()                   # Friction factor AUX pipe [-]
    TRIP_TS = Component.Output()                  # Trip signal for turbine seals
    HPT_SH_Alarm_state = Component.Output()       # Alarm State for HPT Superheat Alarm
    HPT_SH_Trip_state = Component.Output()        # Trip State for HPT Superheat Trip
    HPT_HighTemp_Alarm_state = Component.Output()  # Alarm State for High Temperature Steam entering the HPT
    HPT_HighTemp_TimedTrip_state = Component.Output()  # Timed Trip State for High Temperature Steam entering the HPT
    HPT_HighTemp_Trip_state = Component.Output()  # Trip State for High Temperature Steam entering the HPT
    HPT_ExhPres_Alarm_state = Component.Output()  # Alarm State for High Pressure Leaving HPT
    HPT_ExhPres_Trip_state = Component.Output()   # Trip State for High Pressure Leaving the HPT
    LPT_SH_Alarm_state = Component.Output()       # Alarm State for Low Superheat Levels Entering LPT
    LPT_SH_Trip_state = Component.Output()        # Trip State for Low Superheat Levels Entering LPT
    LPT_HighTemp_Alarm_state = Component.Output()  # Alarm State for High Temperature entering LPT
    LPT_HighTemp_TimedTrip_state = Component.Output()  # Timed Trip State for High Temperature Entering LPT
    LPT_HighTemp_Trip_state = Component.Output()  # Trip State for High Temperature Entering LPT
    # NOTE: outputs 108 and 109 also serve as ifAuto_HPBypass / ifAuto_LPBypass at the very first
    # iteration in the Fortran getIsStartTime block – replicated here as heating-rate outputs only.
    HR_HTF = Component.Output()                   # HTF Heating Rate [Delta K/min]
    HR_Steam = Component.Output()                 # Steam Heating Rate [Delta K/min]
    HighHTF_Temp_Alarm_state = Component.Output()  # Alarm State for High HTF Temperature entering reheater
    HighHTF_Temp_Trip_state = Component.Output()   # Trip State for High HTF Temperature entering reheater
    LowHTF_TempOut_Alarm_state = Component.Output()  # Alarm State for Low HTF Temperature leaving reheater
    LowHTF_TempOut_Trip_state = Component.Output()   # Trip State for Low HTF Temperature leaving reheaters
    HighHTF_Flow_Alarm_state = Component.Output()  # Alarm State for High HTF Flow entering reheaters
    HighHTF_Flow_Trip_state = Component.Output()   # Trip State for High HTF Flow entering reheaters
    HighHTF_Pres_Alarm_state = Component.Output()  # Alarm State for High HTF Pressure entering reheaters
    HighHTF_Pres_Trip_state = Component.Output()   # Trip state for High HTF Pressure entering reheaters
    HighDT_HX_Alarm_state = Component.Output()    # Alarm state for High differential temperature between HTF inlet and Steam inlet
    HighDT_HX_Trip_state = Component.Output()     # Trip state for High differential temperature between HTF inlet and Steam inlet
    HighHTF_HR_Alarm_state = Component.Output()   # Alarm state for High heating rate of HTF entering reheater
    HighHTF_HR_Trip_state = Component.Output()    # Trip state for High heating rate of HTF entering reheater
    HighSteam_HR_Alarm_state = Component.Output()  # Alarm state for High heating rate of Steam leaving reheater
    HighSteam_HR_Trip_state = Component.Output()   # Trip state for High heating rate of Steam leaving reheater

    # -------------------------------------------------------------------------

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        # Dynamic array state (replaces TRNSYS getDynamicArrayValueLastTimestep)
        self._m_HPmain_prev: float = 0.0   # mass of water in HP main piping [kg]
        self._P_HPmain_prev: float = 0.0   # pressure in HP main piping [Pa]
        self._h_HPmain_prev: float = 0.0   # enthalpy in HP main piping [J/kg]
        self._T_HP_pipe: float = 0.0       # temperature of HP pipe metal [K]
        self._m_LPmain_prev: float = 0.0
        self._P_LPmain_prev: float = 0.0
        self._h_LPmain_prev: float = 0.0
        self._T_LP_pipe: float = 0.0
        self._m_Aux_prev: float = 0.0
        self._P_Aux_prev: float = 0.0
        self._h_Aux_prev: float = 0.0
        self._T_Aux_pipe: float = 0.0
        # Static array state (replaces TRNSYS getStaticArrayValue/SetStaticArrayValue)
        # Design conditions stored at first timestep
        self._h_HPT_in_d: float = 0.0
        self._v_HPT_in_d: float = 0.0
        self._s_HPT_in_d: float = 0.0
        self._P_GS_out_d: float = 0.0
        self._h_GS_out_d: float = 0.0
        self._v_GS_out_d: float = 0.0
        self._s_GS_out_d: float = 0.0
        self._h_HPT1_d: float = 0.0
        self._v_HPT1_d: float = 0.0
        self._s_HPT1_d: float = 0.0
        self._h_HPT_exh_d: float = 0.0
        self._v_HPT_exh_d: float = 0.0
        self._s_HPT_exh_d: float = 0.0
        self._h_LPT_in_d: float = 0.0
        self._v_LPT_in_d: float = 0.0
        self._s_LPT_in_d: float = 0.0
        self._h_LPT1_d: float = 0.0
        self._v_LPT1_d: float = 0.0
        self._s_LPT1_d: float = 0.0
        self._h_LPT2_d: float = 0.0
        self._v_LPT2_d: float = 0.0
        self._s_LPT2_d: float = 0.0
        self._h_LPT3_d: float = 0.0
        self._v_LPT3_d: float = 0.0
        self._s_LPT3_d: float = 0.0
        self._h_LPT4_d: float = 0.0
        self._v_LPT4_d: float = 0.0
        self._s_LPT4_d: float = 0.0
        self._h_LPT_exh_d: float = 0.0
        self._v_LPT_exh_d: float = 0.0
        self._s_LPT_exh_d: float = 0.0
        # Temperature history arrays for heating rate alarm calculations
        # Size is 2*N_int; initialized in the first-step block when N_int is known.
        self._N_int: int = 1
        self._T_HTF_history: list = [0.0]  # HTF inlet temperature history (N_int values)
        self._T_Steam_history: list = [0.0]  # steam outlet temperature history (N_int values)

    def initialize(self):
        # convert timestep from hours to seconds
        ts = self.model.timestep * 3600.0
        # find how many seconds are in a timestep
        N = 60.0 / ts
        # Amount of information we need to solve for each heating rate
        self._N_int = math.ceil(N)

        # ---- Read inputs used at start time ----
        HPT_CV_VPi = self.HPT_CV_VPi.v       # Valve Position for high pressure turbine control valves
        HP_Bypass_VPi = self.HP_Bypass_VPi.v  # Valve Position for high pressure bypass line
        HP_Aux_VPi = self.HP_Aux_VPi.v        # Valve Position for valve between hp-main piping and the aux-piping
        HP_warmup_VPi = self.HP_warmup_VPi.v  # Valve Position for high pressure warm-up line
        HP_drain_VPi = self.HP_drain_VPi.v    # Valve Position of drain valve connected to hp-piping
        Aux_DA_VPi = self.Aux_DA_VPi.v        # Valve Position of aux valve leading to DA
        LP_Bypass_VPi = self.LP_Bypass_VPi.v  # Valve Position of low pressure bypass line
        LP_Aux_VPi = self.LP_Aux_VPi.v        # Valve Position for valve between lp main piping and the aux piping
        LP_warmup_VPi = self.LP_warmup_VPi.v  # Valve Position for low pressure warmup line
        LP_drain_VPi = self.LP_drain_VPi.v    # Valve Position for low pressure drain line

        # ---- Calculate the design efficiencies for HPT and LPT using Spencer Cotton Cannon Efficiency Calculations ----
        # CONVERTED-NEEDS UNITS CHECK: eeslib uses SI throughout; all P in Pa and h/s in J/kg
        h_HPT_in_d = fp.enthalpy("water", T=self.T_HPT_in_d.v, P=self.P_HPT_in_d.v)  # CONVERTED-NEEDS UNITS CHECK: eeslib returns J/kg; no *1000 needed
        v_HPT_in_d = fp.spec_vol("water", T=self.T_HPT_in_d.v, P=self.P_HPT_in_d.v)  # CONVERTED-NEEDS UNITS CHECK: eeslib uses Pa; no /1000 needed
        s_HPT_in_d = fp.entropy("water", T=self.T_HPT_in_d.v, P=self.P_HPT_in_d.v)   # CONVERTED-NEEDS UNITS CHECK: eeslib uses Pa; no /1000 needed
        rho_steam = fp.density("water", T=self.T_HPT_in_d.v, P=self.P_HPT_in_d.v)    # CONVERTED-NEEDS UNITS CHECK: eeslib uses Pa; no /1000 needed
        eta_HPT_d = eta_SCC_hpt(self.m_dot_HPT_d.v, self.m_dot_HPT_d.v, self.HPT_parallel_sects.v, self.GS_diameter.v, self.P_HPT_in_d.v, v_HPT_in_d, v_HPT_in_d, self.P_HPT_exh_d.v, self.HPT_CV_NUMBER.v, self.HPT_no_GS.v)

        # ---- Find Control Valve Pressure Losses ----
        # Find the total flow going through one control valve
        m_dot_HPT_CV = self.m_dot_HPT_d.v / self.HPT_CV_NUMBER.v
        CV_controlvalves = PB_CV_data(1, self.HPT_CV_D.v, self.HPT_CV_VPD.v)
        # Volumetric flow entering the valve [m^3/s]
        HPT_Vol_in = m_dot_HPT_CV * v_HPT_in_d
        # Volumetric flow entering valve [GPM]
        HPT_Vol_in_gpm = HPT_Vol_in * 15850.3
        spec_grav = rho_steam / 1000.0
        # pressure drop in psi across valve
        DELTA_P_psi = spec_grav / (CV_controlvalves / HPT_Vol_in_gpm) ** 2.0
        # pressure drop in Pa across valve
        DELTA_P_CV = DELTA_P_psi * 6894.76
        P_GS_out_d = self.P_HPT_in_d.v - DELTA_P_CV
        h_GS_out_d = h_HPT_in_d
        s_GS_out_d = fp.entropy("water", P=P_GS_out_d, h=h_GS_out_d)   # CONVERTED-NEEDS UNITS CHECK: eeslib uses SI; inputs confirmed in Pa/J/kg
        v_GS_out_d = fp.spec_vol("water", P=P_GS_out_d, h=h_GS_out_d)  # CONVERTED-NEEDS UNITS CHECK: eeslib uses SI; inputs confirmed in Pa/J/kg

        # ---- Find Exhaust (Stage 2 Outlet) design pressure, enthalpy, specific volume and entropy ----
        h_HPT_exh_s = fp.enthalpy("water", P=self.P_HPT_exh_d.v, s=s_GS_out_d)  # CONVERTED-NEEDS UNITS CHECK: eeslib uses SI; inputs confirmed in Pa/J/kg
        # CONVERTED-NEEDS UNITS CHECK: eeslib returns J/kg; no *1000 needed
        h_HPT_exh_d = h_HPT_in_d - eta_HPT_d * (h_HPT_in_d - h_HPT_exh_s)
        s_HPT_exh_d = fp.entropy("water", P=self.P_HPT_exh_d.v, h=h_HPT_exh_d)    # CONVERTED-NEEDS UNITS CHECK: eeslib uses SI; inputs confirmed in Pa/J/kg
        v_HPT_exh_d = fp.spec_vol("water", P=self.P_HPT_exh_d.v, h=h_HPT_exh_d)   # CONVERTED-NEEDS UNITS CHECK: eeslib uses SI; inputs confirmed in Pa/J/kg
        x_HPT_exh_d = fp.quality("water", P=self.P_HPT_exh_d.v, h=h_HPT_exh_d)    # CONVERTED-NEEDS UNITS CHECK: eeslib uses SI; inputs confirmed in Pa/J/kg

        # ---- Find HPT Stage 1 Outlet enthalpy, specific volume and entropy ----
        m_elep_d = (h_GS_out_d - h_HPT_exh_d) / (s_GS_out_d - s_HPT_exh_d)
        b_elep_d = h_GS_out_d - m_elep_d * s_GS_out_d
        s_HPT1_d = (s_GS_out_d + s_HPT_exh_d) / 2.0
        s_max = s_HPT_exh_d
        s_min = s_GS_out_d
        LR = 0.25
        whileiterations3 = 0.0
        tol2 = 10.0
        error2 = tol2 + 1.0
        s_HPT1_d_prev = s_HPT1_d
        error2_prev = 0.0
        h_HPT1_d = 0.0
        while abs(error2) > tol2:
            whileiterations3 = whileiterations3 + 1.0
            h_elep_d = m_elep_d * s_HPT1_d + b_elep_d
            h_HPT1_d = fp.enthalpy("water", P=self.P_HPT1_d.v, s=s_HPT1_d)  # CONVERTED-NEEDS UNITS CHECK: eeslib uses SI; inputs confirmed in Pa/J/kg
            error2 = h_elep_d - h_HPT1_d
            if abs(error2) > tol2:
                if whileiterations3 > 1.0:
                    if s_HPT1_d != s_HPT1_d_prev:
                        m_guess = (error2 - error2_prev) / (s_HPT1_d - s_HPT1_d_prev)
                        s_HPT1_d_prev = s_HPT1_d
                        error2_prev = error2
                        if m_guess != 0.0:
                            y_int = error2 - m_guess * s_HPT1_d
                            s_new = max(min(-y_int / m_guess, s_max), s_min)
                            s_HPT1_d = s_HPT1_d + (s_new - s_HPT1_d) * LR
                        else:
                            if error2 > 0.0:
                                s_HPT1_d = min(s_HPT1_d + 1.0, s_max)
                            else:
                                s_HPT1_d = max(s_HPT1_d - 1.0, s_min)
                    else:
                        s_HPT1_d_prev = s_HPT1_d
                        error2_prev = error2
                        if error2 > 0.0:
                            s_HPT1_d = min(s_HPT1_d + 1.0, s_max)
                        else:
                            s_HPT1_d = max(s_HPT1_d - 1.0, s_min)
                else:
                    s_HPT1_d_prev = s_HPT1_d
                    error2_prev = error2
                    if error2 > 0.0:
                        s_HPT1_d = min(s_HPT1_d + 1.0, s_max)
                    else:
                        s_HPT1_d = max(s_HPT1_d - 1.0, s_min)

        v_HPT1_d = fp.spec_vol("water", P=self.P_HPT1_d.v, h=h_HPT1_d)   # CONVERTED-NEEDS UNITS CHECK: eeslib uses SI; inputs confirmed in Pa/J/kg
        x_HPT_exh_d = fp.quality("water", P=self.P_HPT1_d.v, h=h_HPT1_d)  # CONVERTED-NEEDS UNITS CHECK: eeslib uses SI; inputs confirmed in Pa/J/kg

        m_dot_HPTS1 = self.m_dot_HPT_d.v
        m_dot_HPTS2 = m_dot_HPTS1 - self.m_dot_HPT1_d.v

        w_dot_HPT = self.m_dot_HPT_d.v * (h_GS_out_d - h_HPT1_d) + m_dot_HPTS1 * (h_HPT1_d - h_HPT_exh_d)
        m_dot_HPT_exh_d = m_dot_HPTS2 - self.m_dot_HPT2_d.v

        # ---- Steam Separator Calculations ----
        if x_HPT_exh_d < 1.0:
            m_dot_LPT_d = m_dot_HPT_exh_d * x_HPT_exh_d
        else:
            m_dot_LPT_d = m_dot_HPT_exh_d

        # ---- LPT Inlet Conditions ----
        # LPT inlet pressure set equal to HPT exhaust pressure
        P_LPT_in_d = self.P_HPT_exh_d.v
        h_LPT_in_d = fp.enthalpy("water", P=P_LPT_in_d, T=self.T_LPT_in_d.v)   # CONVERTED-NEEDS UNITS CHECK: eeslib uses SI; inputs confirmed in Pa/J/kg
        s_LPT_in_d = fp.entropy("water", P=P_LPT_in_d, T=self.T_LPT_in_d.v)    # CONVERTED-NEEDS UNITS CHECK: eeslib returns J/(kg·K); no *1000 needed
        v_LPT_in_d = fp.spec_vol("water", P=P_LPT_in_d, T=self.T_LPT_in_d.v)   # CONVERTED-NEEDS UNITS CHECK: eeslib uses SI; inputs confirmed in Pa/J/kg

        # ---- LPT Efficiency Calculations ----
        eta_LPT_d = eta_SCC_lpt(m_dot_LPT_d, m_dot_LPT_d, self.LPT_parallel_sects.v, P_LPT_in_d, self.T_LPT_in_d.v, v_LPT_in_d, v_LPT_in_d, self.P_cond_d.v)
        h_LPT_exh_s = fp.enthalpy("water", P=self.P_cond_d.v, s=s_LPT_in_d)  # CONVERTED-NEEDS UNITS CHECK: eeslib uses SI; inputs confirmed in Pa/J/kg
        h_LPT_exh_d = h_LPT_in_d - eta_LPT_d * (h_LPT_in_d - h_LPT_exh_s)
        s_LPT_exh_d = fp.entropy("water", P=self.P_cond_d.v, h=h_LPT_exh_d)   # CONVERTED-NEEDS UNITS CHECK: eeslib returns J/(kg·K); no *1000 needed
        v_LPT_exh_d = fp.spec_vol("water", P=self.P_cond_d.v, h=h_LPT_exh_d)  # CONVERTED-NEEDS UNITS CHECK: eeslib uses SI; inputs confirmed in Pa/J/kg

        h_guess = (h_LPT_in_d - h_LPT_exh_d) / 2.0
        tol = 10.0

        # finding enthalpy out of stage 1
        h_LPT1_d = h_lpt_stage(h_guess, h_LPT_exh_d, h_LPT_in_d, self.P_LPT1_d.v, 0.0, self.LPT_EXP_A0.v, self.LPT_EXP_A1.v, self.LPT_EXP_A2.v, tol)
        s_LPT1_d = fp.entropy("water", P=self.P_LPT1_d.v, h=h_LPT1_d)   # CONVERTED-NEEDS UNITS CHECK: eeslib returns J/(kg·K); no *1000 needed
        v_LPT1_d = fp.spec_vol("water", P=self.P_LPT1_d.v, h=h_LPT1_d)  # CONVERTED-NEEDS UNITS CHECK: eeslib uses SI; inputs confirmed in Pa/J/kg

        # finding enthalpy out of stage 2
        h_LPT2_d = h_lpt_stage(h_guess, h_LPT_exh_d, h_LPT_in_d, self.P_LPT2_d.v, 0.0, self.LPT_EXP_A0.v, self.LPT_EXP_A1.v, self.LPT_EXP_A2.v, tol)
        s_LPT2_d = fp.entropy("water", P=self.P_LPT2_d.v, h=h_LPT2_d)   # CONVERTED-NEEDS UNITS CHECK: eeslib returns J/(kg·K); no *1000 needed
        v_LPT2_d = fp.spec_vol("water", P=self.P_LPT2_d.v, h=h_LPT2_d)  # CONVERTED-NEEDS UNITS CHECK: eeslib uses SI; inputs confirmed in Pa/J/kg

        # finding enthalpy out of stage 3
        h_LPT3_d = h_lpt_stage(h_guess, h_LPT_exh_d, h_LPT_in_d, self.P_LPT3_d.v, 0.0, self.LPT_EXP_A0.v, self.LPT_EXP_A1.v, self.LPT_EXP_A2.v, tol)
        s_LPT3_d = fp.entropy("water", P=self.P_LPT3_d.v, h=h_LPT3_d)   # CONVERTED-NEEDS UNITS CHECK: eeslib returns J/(kg·K); no *1000 needed
        v_LPT3_d = fp.spec_vol("water", P=self.P_LPT3_d.v, h=h_LPT3_d)  # CONVERTED-NEEDS UNITS CHECK: eeslib uses SI; inputs confirmed in Pa/J/kg

        # finding enthalpy out of stage 4
        h_LPT4_d = h_lpt_stage(h_guess, h_LPT_exh_d, h_LPT_in_d, self.P_LPT4_d.v, 0.0, self.LPT_EXP_A0.v, self.LPT_EXP_A1.v, self.LPT_EXP_A2.v, tol)
        s_LPT4_d = fp.entropy("water", P=self.P_LPT4_d.v, h=h_LPT4_d)   # CONVERTED-NEEDS UNITS CHECK: eeslib returns J/(kg·K); no *1000 needed
        v_LPT4_d = fp.spec_vol("water", P=self.P_LPT4_d.v, h=h_LPT4_d)  # CONVERTED-NEEDS UNITS CHECK: eeslib uses SI; inputs confirmed in Pa/J/kg

        # Stage Design Inlet Mass Flows
        m_dot_LPTS1 = m_dot_LPT_d
        m_dot_LPTS2 = m_dot_LPT_d - self.m_dot_LPT1_d.v
        m_dot_LPTS3 = m_dot_LPTS2 - self.m_dot_LPT2_d.v
        m_dot_LPTS4 = m_dot_LPTS3 - self.m_dot_LPT3_d.v
        m_dot_LPT_exh = m_dot_LPTS4 - self.m_dot_LPT4_d.v

        # Finding rated power
        w_dot_LPT = m_dot_LPTS1 * (h_LPT_in_d - h_LPT1_d)
        w_dot_LPT = w_dot_LPT + m_dot_LPTS2 * (h_LPT1_d - h_LPT2_d)
        w_dot_LPT = w_dot_LPT + m_dot_LPTS3 * (h_LPT2_d - h_LPT3_d)
        w_dot_LPT = w_dot_LPT + m_dot_LPTS4 * (h_LPT3_d - h_LPT4_d)
        w_dot_LPT = w_dot_LPT + m_dot_LPT_exh * (h_LPT4_d - h_LPT_exh_d)
        W_dot_total = w_dot_LPT + w_dot_HPT

        # ---- High pressure main initial enthalpy ----
        # CONVERTED-NEEDS UNITS CHECK: eeslib uses Pa; no /1000 needed
        T_sat_HP = fp.temperature("water", P=self.P_HPmain_ini.v, Q=1.0)  # CONVERTED-NEEDS UNITS CHECK: eeslib uses Pa; no /1000 needed
        h_sat_HP = fp.enthalpy("water", P=self.P_HPmain_ini.v, Q=1.0)     # CONVERTED-NEEDS UNITS CHECK: eeslib returns J/kg; no *1000 needed
        rho_HP = fp.density("water", P=self.P_HPmain_ini.v, Q=1.0)        # CONVERTED-NEEDS UNITS CHECK: eeslib uses SI; inputs confirmed in Pa/J/kg
        if self.T_HPmain_ini.v > T_sat_HP:
            h_hpmain_ini = fp.enthalpy("water", P=self.P_HPmain_ini.v, T=self.T_HPmain_ini.v)  # CONVERTED-NEEDS UNITS CHECK: eeslib returns J/kg; no *1000 needed
            rho_HP = fp.density("water", P=self.P_HPmain_ini.v, T=self.T_HPmain_ini.v)         # CONVERTED-NEEDS UNITS CHECK: eeslib uses SI; inputs confirmed in Pa/J/kg
            m_hpmain_ini = rho_HP * math.pi / 4.0 * self.D_HPmain.v ** 2.0 * self.Length_HPmain.v
        else:
            h_hpmain_ini = h_sat_HP
            m_hpmain_ini = rho_HP * math.pi / 4.0 * self.D_HPmain.v ** 2.0 * self.Length_HPmain.v

        # ---- Low pressure main initial enthalpy and mass ----
        T_sat_LP = fp.temperature("water", P=self.P_LPmain_ini.v, Q=1.0)  # CONVERTED-NEEDS UNITS CHECK: eeslib uses Pa; no /1000 needed
        h_sat_LP = fp.enthalpy("water", P=self.P_LPmain_ini.v, Q=1.0)     # CONVERTED-NEEDS UNITS CHECK: eeslib returns J/kg; no *1000 needed
        rho_LP = fp.density("water", P=self.P_LPmain_ini.v, Q=1.0)        # CONVERTED-NEEDS UNITS CHECK: eeslib uses SI; inputs confirmed in Pa/J/kg
        if self.T_LPmain_ini.v > T_sat_LP:
            h_lpmain_ini = fp.enthalpy("water", P=self.P_LPmain_ini.v, T=self.T_LPmain_ini.v)  # CONVERTED-NEEDS UNITS CHECK: eeslib returns J/kg; no *1000 needed
            rho_LP = fp.density("water", P=self.P_LPmain_ini.v, T=self.T_LPmain_ini.v)         # CONVERTED-NEEDS UNITS CHECK: eeslib uses SI; inputs confirmed in Pa/J/kg
            m_lpmain_ini = rho_LP * math.pi / 4.0 * self.D_LPmain.v ** 2.0 * self.Length_LPmain.v
        else:
            h_lpmain_ini = h_sat_LP
            m_lpmain_ini = rho_LP * math.pi / 4.0 * self.D_LPmain.v ** 2.0 * self.Length_LPmain.v

        # ---- Aux line initial enthalpy and mass ----
        T_sat_Aux = fp.temperature("water", P=self.P_Aux_ini.v, Q=1.0)  # CONVERTED-NEEDS UNITS CHECK: eeslib uses Pa; no /1000 needed
        h_sat_Aux = fp.enthalpy("water", P=self.P_Aux_ini.v, Q=1.0)     # CONVERTED-NEEDS UNITS CHECK: eeslib returns J/kg; no *1000 needed
        rho_Aux = fp.density("water", P=self.P_Aux_ini.v, Q=1.0)        # CONVERTED-NEEDS UNITS CHECK: eeslib uses SI; inputs confirmed in Pa/J/kg
        if self.T_Aux_ini.v > T_sat_Aux:
            h_aux_ini = fp.enthalpy("water", P=self.P_Aux_ini.v, T=self.T_Aux_ini.v)  # CONVERTED-NEEDS UNITS CHECK: eeslib returns J/kg; no *1000 needed
            rho_Aux = fp.density("water", P=self.P_Aux_ini.v, T=self.T_Aux_ini.v)     # CONVERTED-NEEDS UNITS CHECK: eeslib uses SI; inputs confirmed in Pa/J/kg
            m_aux_ini = rho_Aux * math.pi / 4.0 * self.D_AuxLine.v ** 2.0 * self.Length_AuxLine.v
        else:
            h_aux_ini = h_sat_Aux
            # NOTE: Fortran uses D_lpmain and Length_lPmain here (appears to be a bug in the Fortran source)
            m_aux_ini = rho_Aux * math.pi / 4.0 * self.D_LPmain.v ** 2.0 * self.Length_LPmain.v

        # ---- Store design values in static instance attributes ----
        self._h_HPT_in_d = h_HPT_in_d    # Design Specific Enthalpy entering the HPT [J/kg]
        self._v_HPT_in_d = v_HPT_in_d    # Design Specific Volume entering the HPT [m^3/kg]
        self._s_HPT_in_d = s_HPT_in_d    # Design Specific Entropy entering the HPT [J/kg-K]
        self._P_GS_out_d = P_GS_out_d    # Design outlet pressure through the HPT control valves [Pa]
        self._h_GS_out_d = h_GS_out_d    # Design Specific enthalpy out of the control valves [J/kg]
        self._v_GS_out_d = v_GS_out_d    # Design Specific volume out of the control valves [m^3/kg]
        self._s_GS_out_d = s_GS_out_d    # Design Specific entropy out of the control valves [J/kg-K]
        self._h_HPT1_d = h_HPT1_d        # Design specific enthalpy out of high pressure turbine stage 1 [J/kg]
        self._v_HPT1_d = v_HPT1_d        # Design specific volume out of high pressure turbine stage 1 [m^3/kg]
        self._s_HPT1_d = s_HPT1_d        # Design specific entropy out of high pressure turbine stage 1 [J/kg-K]
        self._h_HPT_exh_d = h_HPT_exh_d  # Design specific enthalpy out of the high pressure turbine [J/kg]
        self._v_HPT_exh_d = v_HPT_exh_d  # Design specific volume out of high pressure turbine [m^3/kg]
        self._s_HPT_exh_d = s_HPT_exh_d  # Design specific entropy out of high pressure turbine [J/kg-K]
        self._h_LPT_in_d = h_LPT_in_d    # Design specific enthalpy into the low pressure turbine [J/kg]
        self._v_LPT_in_d = v_LPT_in_d    # Design specific volume into low pressure turbine [m^3/kg]
        self._s_LPT_in_d = s_LPT_in_d    # Design specific entropy into the low pressure turbine [J/kg-K]
        self._h_LPT1_d = h_LPT1_d        # Design specific enthalpy out of low pressure turbine stage 1 [J/kg]
        self._v_LPT1_d = v_LPT1_d        # Design specific volume out of low pressure turbine stage 1 [m^3/kg]
        self._s_LPT1_d = s_LPT1_d        # Design specific entropy out of low pressure turbine stage 1 [J/kg-K]
        self._h_LPT2_d = h_LPT2_d        # Design specific enthalpy out of low pressure turbine stage 2 [J/kg]
        self._v_LPT2_d = v_LPT2_d        # Design specific volume out of low pressure turbine stage 2 [m^3/kg]
        self._s_LPT2_d = s_LPT2_d        # Design specific entropy out of low pressure turbine stage 2 [J/kg-K]
        self._h_LPT3_d = h_LPT3_d        # Design specific enthalpy out of low pressure turbine stage 3 [J/kg]
        self._v_LPT3_d = v_LPT3_d        # Design specific volume out of low pressure turbine stage 3 [m^3/kg]
        self._s_LPT3_d = s_LPT3_d        # Design specific entropy out of low pressure turbine stage 3 [J/kg-K]
        self._h_LPT4_d = h_LPT4_d        # Design specific enthalpy out of low pressure turbine stage 4 [J/kg]
        self._v_LPT4_d = v_LPT4_d        # Design specific volume out of low pressure turbine stage 4 [m^3/kg]
        self._s_LPT4_d = s_LPT4_d        # Design specific entropy out of low pressure turbine stage 4 [J/kg-K]
        self._h_LPT_exh_d = h_LPT_exh_d  # Design specific enthalpy out of low pressure turbine [J/kg]
        self._v_LPT_exh_d = v_LPT_exh_d  # Design specific volume out of low pressure turbine [m^3/kg]
        self._s_LPT_exh_d = s_LPT_exh_d  # Design specific entropy out of low pressure turbine [J/kg-K]

        # ---- Store dynamic (piping state) initial values ----
        self._m_HPmain_prev = m_hpmain_ini   # mass of water in the main piping line before the HP turbine
        self._P_HPmain_prev = self.P_HPmain_ini.v  # Initial pressure in the main piping line before the HP turbine
        self._h_HPmain_prev = h_hpmain_ini   # Initial enthalpy in the main piping line before the HP turbine
        self._T_HP_pipe = self.T_HPpipe_ini.v  # Temperature of the HP main metal piping at the start of the simulation
        self._m_LPmain_prev = m_lpmain_ini   # mass of water in the main piping line before the LP turbine
        self._P_LPmain_prev = self.P_LPmain_ini.v  # Initial pressure in the main piping line before the LP turbine
        self._h_LPmain_prev = h_lpmain_ini   # Initial enthalpy of the LP main metal piping at the start of the simulation
        self._T_LP_pipe = self.T_LPpipe_ini.v  # Temperature of the LP main metal piping at the start of the simulation
        self._m_Aux_prev = m_aux_ini         # mass of water in the aux piping line leading to the turbine seals and DA
        self._P_Aux_prev = self.P_Aux_ini.v  # Initial pressure in the aux piping line leading to the turbine seals and DA
        self._h_Aux_prev = h_aux_ini         # Initial temperature in the aux piping line leading to the turbine seals and DA
        self._T_Aux_pipe = self.T_Auxpipe_ini.v  # Temperature of the aux piping line at the beginning of the simulation

        # ---- Initialize temperature history arrays for heating rate calculations ----
        self._T_HTF_history = [self.HTF_T_in.v] * self._N_int
        self._T_Steam_history = [self.T_Aux_ini.v] * self._N_int

        # ---- Set initial output values ----
        self.Turbine_ON_out.v = self.Turbine_ON.v           # HPT & LPT Turbine On
        self.HPT_CV_VPo.v = HPT_CV_VPi                      # HPT Control Valve Position output
        self.HP_Bypass_VPo.v = HP_Bypass_VPi                # High Pressure Bypass Valve Position output
        self.HP_Aux_VPo.v = HP_Aux_VPi
        self.HP_warmup_VPo.v = HP_warmup_VPi
        self.HP_drain_VPo.v = HP_drain_VPi
        self.Aux_DA_VPo.v = Aux_DA_VPi
        self.LP_Bypass_VPo.v = LP_Bypass_VPi
        self.LP_Aux_VPo.v = LP_Aux_VPi
        self.LP_warmup_VPo.v = LP_warmup_VPi
        self.LP_drain_VPo.v = LP_drain_VPi
        self.m_dot_SGT.v = 60.0                              # output to input to SD type (ESOL6003)
        self.m_dot_cond.v = 0.0                              # output to condenser type (ESOL6007)
        self.h_cond.v = 0.0                                  # output to condenser type (ESOL6007)
        self.m_dot_DA.v = 0.0                                # output to DA type (ESOL6011)
        self.h_DA.v = 0.0                                    # output to DA type (ESOL6011)
        self.W_dot_tot.v = W_dot_total                       # Actual Turbine Outlet Power
        self.m_dot_HTF_out.v = self.m_dot_HTF.v             # Reheater - HTF Mass Flow Out
        self.Vol_dot_HTF.v = 0.0                             # Reheater - HTF Volumetric Flow out
        self.HTF_P_out.v = self.HTF_P_in.v                  # Reheater - HTF Pressure out
        self.T_HTF_out.v = self.HTF_T_in.v                  # Reheater - HTF Temperature out
        self.m_dot_HPT_in.v = self.m_dot_HPT_d.v            # Mass Flow rate entering HPT
        self.P_GS_out.v = P_GS_out_d                        # HPT - GS Outlet Pressure
        self.h_GS_out.v = h_GS_out_d                        # HPT - GS Outlet Enthalpy
        self.m_dot_HPTS1.v = m_dot_HPTS1                    # HPT-Stage 1 Inlet Mass flow
        self.P_HPT1.v = self.P_HPT1_d.v                     # HPT-Stage 1 Outlet Pressure
        self.H_HPT1.v = h_HPT1_d                            # HPT-Stage 1 Outlet Enthalpy
        self.m_dot_HPTS2.v = m_dot_HPTS2                    # HPT-Stage 2 Inlet Mass Flow
        self.P_HPT2.v = 0.0                                  # HPT-Stage 2 Outlet Pressure (P_HPT2_d uninitialized in Fortran)
        self.H_HPT2.v = 0.0                                  # HPT-Stage 2 Outlet Enthalpy (H_HPT2_d uninitialized in Fortran)
        self.m_dot_HPT_exh.v = m_dot_HPT_exh_d              # HPT-Exhaust Mass Flow
        self.P_HPT_exh.v = self.P_HPT_exh_d.v               # HPT-Exhaust Pressure
        self.H_HPT_exh.v = h_HPT_exh_d                      # HPT-Exhaust Enthalpy
        self.T_HPT_exh.v = 0.0                               # HPT-Exhaust Temperature
        self.m_dot_LPT_in.v = m_dot_LPTS1                   # LPT-Stage 1 Mass Flow
        self.P_LPT1.v = self.P_LPT1_d.v                     # LPT-Stage 1 Outlet Pressure
        self.T_LPT1.v = 0.0                                  # LPT-Stage 1 Outlet Temperature
        self.H_LPT1.v = h_LPT1_d                            # LPT-Stage 1 Outlet Enthalpy
        self.m_dot_LPTS2.v = m_dot_LPTS2                    # LPT-Stage 2 Mass Flow
        self.P_LPT2.v = self.P_LPT2_d.v                     # LPT-Stage 2 Outlet Pressure
        self.T_LPT2.v = 0.0                                  # LPT-Stage 2 Outlet Temperature
        self.H_LPT2.v = h_LPT2_d                            # LPT-Stage 2 Outlet Enthalpy
        self.m_dot_LPTS3.v = m_dot_LPTS3                    # LPT-Stage 3 Mass Flow
        self.P_LPT3.v = self.P_LPT3_d.v                     # LPT-Stage 3 Outlet Pressure
        self.T_LPT3.v = 0.0                                  # LPT-Stage 3 Outlet Temperature
        self.H_LPT3.v = h_LPT3_d                            # LPT-Stage 3 Outlet Enthalpy
        self.m_dot_LPTS4.v = m_dot_LPTS4                    # LPT-Stage 4 Mass Flow
        self.T_LPT4.v = 0.0                                  # LPT-Stage 4 Outlet Temperature
        self.P_LPT4.v = self.P_LPT4_d.v                     # LPT-Stage 4 Outlet Pressure
        self.H_LPT4.v = h_LPT4_d                            # LPT-Stage 4 Outlet Enthalpy
        self.m_dot_LPT_exh.v = m_dot_LPT_exh                # LPT Exhaust Mass Flow
        self.Vol_dot_LPT_exh.v = 0.0                         # LPT Exhaust Volumetric Flow
        self.T_LPT_exh.v = 0.0                               # LPT Exhaust Temperature
        self.H_LPT_exh.v = h_LPT_exh_d                      # LPT Exhaust Enthalpy
        self.P_HPmain.v = self.P_HPmain_ini.v
        self.P_LPmain.v = self.P_LPmain_ini.v
        self.P_AUX.v = self.P_Aux_ini.v
        self.ff_HPmain.v = 0.1                               # for HPmain
        self.ff_LPmain.v = 0.1                               # for LPmain
        self.ff_AUX.v = 0.1                                  # for AUX

    def calculate(self):
        ts = self.model.timestep * 3600.0  # convert timestep from hours to seconds
        N = 60.0 / ts                       # how many timesteps fit in a minute
        N_int = math.ceil(N)               # number of intervals for heating-rate history

        # -------------------------------------------------------------------------
        # END-OF-TIMESTEP BLOCK (replaces Fortran getIsEndOfTimestep())
        # -------------------------------------------------------------------------
        if self.model.is_converged:
            # ----------------------------------------------------------------
            # END-OF-TIMESTEP ALARM AND TRIP CHECKS
            # (replaces Fortran getIsEndOfTimestep() block, lines 113–470)
            # ----------------------------------------------------------------

            # Initialise history arrays on first step
            if self.model.is_first_step:
                self._T_HTF_history = [self.HTF_T_in.v] * N_int
                self._T_Steam_history = [self.T_SS_steam.v] * N_int

            # ---- Turbine alarms (Fortran lines 160–302) ----
            if self.Turbine_ON.v != 1.0:
                # Turbine OFF — clear all turbine alarm/trip outputs
                self.HPT_SH_Alarm_state.v = 0.0
                self.HPT_SH_Trip_state.v = 0.0
                self.HPT_HighTemp_Alarm_state.v = 0.0
                self.HPT_HighTemp_TimedTrip_state.v = 0.0
                self.HPT_HighTemp_Trip_state.v = 0.0
                self.HPT_ExhPres_Alarm_state.v = 0.0
                self.HPT_ExhPres_Trip_state.v = 0.0
                self.LPT_SH_Alarm_state.v = 0.0
                self.LPT_SH_Trip_state.v = 0.0
                self.LPT_HighTemp_Alarm_state.v = 0.0
                self.LPT_HighTemp_TimedTrip_state.v = 0.0
                self.LPT_HighTemp_Trip_state.v = 0.0
            else:
                # Turbine ON — check alarms and trips
                P_HPmain_prev_eon = self.P_HPmain.v    # output 71
                T_hpmain_prev_eon = self.T_HPmain.v    # output 72
                P_LPmain_prev_eon = self.P_LPmain.v    # output 75
                T_LPmain_eon = self.T_LPmain.v         # output 76
                W_dot_total_eon = self.W_dot_tot.v

                # HPT Superheat alarm/trip (Fortran lines 170–217)
                T_sat_eon = fp.temperature("water", P=P_HPmain_prev_eon, Q=0.0)  # CONVERTED-NEEDS UNITS CHECK: eeslib uses SI; inputs confirmed in Pa/J/kg
                Superheat_eon = T_hpmain_prev_eon - T_sat_eon
                if W_dot_total_eon > self.Partial_Load.v:
                    # Full load conditions
                    if Superheat_eon > self.HPT_SH_Alarm_FL.v:
                        self.HPT_SH_Alarm_state.v = 0.0
                        self.HPT_SH_Trip_state.v = 0.0
                    else:
                        self.HPT_SH_Alarm_state.v = 1.0
                        self.HPT_SH_Trip_state.v = 0.0 if Superheat_eon > self.HPT_SH_Trip_FL.v else 1.0
                else:
                    # Partial load conditions
                    if Superheat_eon > self.HPT_SH_Alarm_PL.v:
                        self.HPT_SH_Alarm_state.v = 0.0
                        self.HPT_SH_Trip_state.v = 0.0
                    else:
                        self.HPT_SH_Alarm_state.v = 1.0
                        self.HPT_SH_Trip_state.v = 0.0 if Superheat_eon > self.HPT_SH_Trip_PL.v else 1.0

                # HPT High temperature alarm/timed-trip/trip (Fortran lines 219–241)
                if T_hpmain_prev_eon < self.HPT_HighTemp_Alarm.v:
                    self.HPT_HighTemp_Alarm_state.v = 0.0
                    self.HPT_HighTemp_TimedTrip_state.v = 0.0
                    self.HPT_HighTemp_Trip_state.v = 0.0
                else:
                    self.HPT_HighTemp_Alarm_state.v = 1.0
                    if T_hpmain_prev_eon < self.HPT_HighTemp_TimedTrip.v:
                        self.HPT_HighTemp_TimedTrip_state.v = 0.0
                        self.HPT_HighTemp_Trip_state.v = 0.0
                    else:
                        self.HPT_HighTemp_TimedTrip_state.v = 1.0
                        self.HPT_HighTemp_Trip_state.v = (
                            0.0 if T_hpmain_prev_eon < self.HPT_HighTemp_Trip.v else 1.0
                        )

                # HPT Exhaust pressure alarm/trip (Fortran lines 243–255)
                if P_LPmain_prev_eon < self.HPT_ExhPres_Alarm.v:
                    self.HPT_ExhPres_Alarm_state.v = 0.0
                    self.HPT_ExhPres_Trip_state.v = 0.0
                else:
                    self.HPT_ExhPres_Alarm_state.v = 1.0
                    self.HPT_ExhPres_Trip_state.v = (
                        0.0 if P_LPmain_prev_eon < self.HPT_ExhPres_Trip.v else 1.0
                    )

                # LPT Superheat alarm/trip (Fortran lines 257–276)
                T_sat_lpt_eon = fp.temperature("water", P=P_LPmain_prev_eon, Q=0.0)  # CONVERTED-NEEDS UNITS CHECK: eeslib uses SI; inputs confirmed in Pa/J/kg
                Superheat_lpt = T_LPmain_eon - T_sat_lpt_eon
                if Superheat_lpt > self.LPT_SH_Alarm.v:
                    self.LPT_SH_Alarm_state.v = 0.0
                    self.LPT_SH_Trip_state.v = 0.0
                else:
                    self.LPT_SH_Alarm_state.v = 1.0
                    self.LPT_SH_Trip_state.v = 0.0 if Superheat_lpt > self.LPT_SH_Trip.v else 1.0

                # LPT High temperature alarm/timed-trip/trip (Fortran lines 278–300)
                if T_LPmain_eon < self.LPT_HighTemp_Alarm.v:
                    self.LPT_HighTemp_Alarm_state.v = 0.0
                    self.LPT_HighTemp_TimedTrip_state.v = 0.0
                    self.LPT_HighTemp_Trip_state.v = 0.0
                else:
                    self.LPT_HighTemp_Alarm_state.v = 1.0
                    if T_LPmain_eon < self.LPT_HighTemp_TimedTrip.v:
                        self.LPT_HighTemp_TimedTrip_state.v = 0.0
                        self.LPT_HighTemp_Trip_state.v = 0.0
                    else:
                        self.LPT_HighTemp_TimedTrip_state.v = 1.0
                        self.LPT_HighTemp_Trip_state.v = (
                            0.0 if T_LPmain_eon < self.LPT_HighTemp_Trip.v else 1.0
                        )

            # ---- Reheater alarms (always, Fortran lines 308–470) ----
            HTF_T_in_eon = self.HTF_T_in.v
            T_SS_eon = self.T_SS_steam.v

            # High HTF Temperature entering reheater
            if HTF_T_in_eon < self.HighHTF_Temp_Alarm.v:
                self.HighHTF_Temp_Alarm_state.v = 0.0
                self.HighHTF_Temp_Trip_state.v = 0.0
            else:
                self.HighHTF_Temp_Alarm_state.v = 1.0
                self.HighHTF_Temp_Trip_state.v = (
                    0.0 if HTF_T_in_eon < self.HighHTF_Temp_Trip.v else 1.0
                )

            # Low HTF Temperature leaving reheater
            T_HTF_out_eon = self.T_HTF_out.v
            if T_HTF_out_eon > self.LowHTF_TempOut_Alarm.v:
                self.LowHTF_TempOut_Alarm_state.v = 0.0
                self.LowHTF_TempOut_Trip_state.v = 0.0
            else:
                self.LowHTF_TempOut_Alarm_state.v = 1.0
                self.LowHTF_TempOut_Trip_state.v = (
                    0.0 if T_HTF_out_eon > self.LowHTF_TempOut_Trip.v else 1.0
                )

            # High HTF Flow entering reheater
            m_dot_HTF_eon = self.m_dot_HTF.v
            if m_dot_HTF_eon < self.HighHTF_Flow_Alarm.v:
                self.HighHTF_Flow_Alarm_state.v = 0.0
                self.HighHTF_Flow_Trip_state.v = 0.0
            else:
                self.HighHTF_Flow_Alarm_state.v = 1.0
                self.HighHTF_Flow_Trip_state.v = (
                    0.0 if m_dot_HTF_eon < self.HighHTF_Flow_Trip.v else 1.0
                )

            # High HTF Pressure entering reheater
            HTF_P_in_eon = self.HTF_P_in.v
            if HTF_P_in_eon < self.HighHTF_Pres_Alarm.v:
                self.HighHTF_Pres_Alarm_state.v = 0.0
                self.HighHTF_Pres_Trip_state.v = 0.0
            else:
                self.HighHTF_Pres_Alarm_state.v = 1.0
                self.HighHTF_Pres_Trip_state.v = (
                    0.0 if HTF_P_in_eon < self.HighHTF_Pres_Trip.v else 1.0
                )

            # High differential temperature (HTF inlet vs steam inlet)
            DT_HX_eon = HTF_T_in_eon - T_SS_eon
            if DT_HX_eon < self.HighDT_HX_Alarm.v:
                self.HighDT_HX_Alarm_state.v = 0.0
                self.HighDT_HX_Trip_state.v = 0.0
            else:
                self.HighDT_HX_Alarm_state.v = 1.0
                self.HighDT_HX_Trip_state.v = (
                    0.0 if DT_HX_eon < self.HighDT_HX_Trip.v else 1.0
                )

            # ---- Heating rate calcs (Fortran lines 395–470) ----
            if not self.model.is_first_step:
                # HTF heating rate [K/min]
                T_hist_HTF = self._T_HTF_history   # N_int stored values
                HR_HTF = 0.0
                for _idx in range(1, N_int):
                    HR_HTF += (T_hist_HTF[_idx] - T_hist_HTF[_idx - 1]) / ts
                HR_HTF += (HTF_T_in_eon - T_hist_HTF[N_int - 1]) / ts
                HR_HTF = HR_HTF / N_int * 60.0     # [K/min]
                self.HR_HTF.v = HR_HTF
                if abs(HR_HTF) < self.HighHTF_HR_Alarm.v:
                    self.HighHTF_HR_Alarm_state.v = 0.0
                    self.HighHTF_HR_Trip_state.v = 0.0
                else:
                    self.HighHTF_HR_Alarm_state.v = 1.0
                    self.HighHTF_HR_Trip_state.v = (
                        0.0 if abs(HR_HTF) < self.HighHTF_HR_Trip.v else 1.0
                    )
                # Shift HTF history: drop oldest, append current
                self._T_HTF_history = self._T_HTF_history[1:] + [HTF_T_in_eon]

                # Steam heating rate [K/min]
                T_hist_Steam = self._T_Steam_history
                HR_Steam = 0.0
                for _idx in range(1, N_int):
                    HR_Steam += (T_hist_Steam[_idx] - T_hist_Steam[_idx - 1]) / ts
                HR_Steam += (T_SS_eon - T_hist_Steam[N_int - 1]) / ts
                HR_Steam = HR_Steam / N_int * 60.0
                self.HR_Steam.v = HR_Steam
                if abs(HR_Steam) < self.HighSteam_HR_Alarm.v:
                    self.HighSteam_HR_Alarm_state.v = 0.0
                    self.HighSteam_HR_Trip_state.v = 0.0
                else:
                    self.HighSteam_HR_Alarm_state.v = 1.0
                    self.HighSteam_HR_Trip_state.v = (
                        0.0 if abs(HR_Steam) < self.HighSteam_HR_Trip.v else 1.0
                    )
                # Shift Steam history
                self._T_Steam_history = self._T_Steam_history[1:] + [T_SS_eon]

            return

        # -------------------------------------------------------------------------
        # READ INPUTS
        # -------------------------------------------------------------------------
        Turbine_ON = self.Turbine_ON.v
        ifAuto_HPBypass = self.ifAuto_HPBypass.v
        ifAuto_LPBypass = self.ifAuto_LPBypass.v
        ifAuto_HPAux = self.ifAuto_HPAux.v
        ifAuto_LPAux = self.ifAuto_LPAux.v
        ifAuto_DAaux = self.ifAuto_DAaux.v
        ifAuto_HPTCV = self.ifAuto_HPTCV.v
        P_DA = self.P_DA.v
        HPT_CV_VPi = self.HPT_CV_VPi.v
        HP_Bypass_VPi = self.HP_Bypass_VPi.v
        HP_Aux_VPi = self.HP_Aux_VPi.v
        HP_warmup_VPi = self.HP_warmup_VPi.v
        HP_drain_VPi = self.HP_drain_VPi.v
        Aux_DA_VPi = self.Aux_DA_VPi.v
        LP_Bypass_VPi = self.LP_Bypass_VPi.v
        LP_Aux_VPi = self.LP_Aux_VPi.v
        LP_warmup_VPi = self.LP_warmup_VPi.v
        LP_drain_VPi = self.LP_drain_VPi.v
        P_SGT_in = self.P_SGT_in.v
        h_SGT_in = self.h_SGT_in.v
        m_dot_HTF = self.m_dot_HTF.v
        HTF_P_in = self.HTF_P_in.v
        HTF_T_in = self.HTF_T_in.v
        P_cond = self.P_cond.v
        m_dot_HPT1 = self.m_dot_HPT1.v
        m_dot_HPT2 = self.m_dot_HPT2.v
        m_dot_LPT1 = self.m_dot_LPT1.v
        m_dot_LPT2 = self.m_dot_LPT2.v
        m_dot_LPT3 = self.m_dot_LPT3.v
        m_dot_LPT4 = self.m_dot_LPT4.v
        PID_HPbypass = self.PID_HPbypass.v
        PID_LPbypass = self.PID_LPbypass.v
        PID_HPaux = self.PID_HPaux.v
        PID_LPaux = self.PID_LPaux.v
        PID_DAaux = self.PID_DAaux.v
        PID_HPTCV = self.PID_HPTCV.v

        # Commonly used constants
        P_ts_min = 500000.0   # [Pa]
        P_atm = 101325.0      # [Pa]
        dh = 25000.0          # enthalpy perturbation for partial derivatives [J/kg]
        dP = 25000.0          # pressure perturbation for partial derivatives [Pa]

        # -------------------------------------------------------------------------
        # FIRST ITERATION: UPDATE VALVE POSITIONS
        # (replaces Fortran getTimestepIteration() == 0 block, lines 1049–1111)
        # -------------------------------------------------------------------------
        # Retrieve current valve-position outputs (carried from previous iteration/timestep)
        HPT_CV_VPo = self.HPT_CV_VPo.v
        HP_Bypass_VPo = self.HP_Bypass_VPo.v
        HP_Aux_VPo = self.HP_Aux_VPo.v
        HP_warmup_VPo = self.HP_warmup_VPo.v
        HP_drain_VPo = self.HP_drain_VPo.v
        Aux_DA_VPo = self.Aux_DA_VPo.v
        LP_Bypass_VPo = self.LP_Bypass_VPo.v
        LP_Aux_VPo = self.LP_Aux_VPo.v
        LP_warmup_VPo = self.LP_warmup_VPo.v
        LP_drain_VPo = self.LP_drain_VPo.v

        if self.model.iteration == 0:
            # Ramp all valve positions toward targets
            # (replaces Fortran VP_new block, lines 1057–1111)
            if ifAuto_HPTCV != 1.0:
                HPT_CV_VPo = VP_new(HPT_CV_VPo, HPT_CV_VPi, self.HPT_CV_VS.v, ts)
            else:
                HPT_CV_VPo = VP_new(HPT_CV_VPo, PID_HPTCV, self.HPT_CV_VS.v, ts)

            if ifAuto_HPBypass != 1.0:
                HP_Bypass_VPo = VP_new(HP_Bypass_VPo, HP_Bypass_VPi, self.HP_Bypass_vs.v, ts)
            else:
                HP_Bypass_VPo = VP_new(HP_Bypass_VPo, PID_HPbypass, self.HP_Bypass_vs.v, ts)

            if ifAuto_LPBypass != 1.0:
                LP_Bypass_VPo = VP_new(LP_Bypass_VPo, LP_Bypass_VPi, self.LP_Bypass_vs.v, ts)
            else:
                LP_Bypass_VPo = VP_new(LP_Bypass_VPo, PID_LPbypass, self.LP_Bypass_vs.v, ts)

            if ifAuto_HPAux != 1.0:
                HP_Aux_VPo = VP_new(HP_Aux_VPo, HP_Aux_VPi, self.HP_Aux_vs.v, ts)
            else:
                HP_Aux_VPo = VP_new(HP_Aux_VPo, PID_HPaux, self.HP_Aux_vs.v, ts)

            if ifAuto_LPAux != 1.0:
                LP_Aux_VPo = VP_new(LP_Aux_VPo, LP_Aux_VPi, self.LP_Aux_vs.v, ts)
            else:
                LP_Aux_VPo = VP_new(LP_Aux_VPo, PID_LPaux, self.LP_Aux_vs.v, ts)

            if ifAuto_DAaux != 1.0:
                Aux_DA_VPo = VP_new(Aux_DA_VPo, Aux_DA_VPi, self.Aux_DA_vs.v, ts)
            else:
                Aux_DA_VPo = VP_new(Aux_DA_VPo, PID_DAaux, self.Aux_DA_vs.v, ts)

            # Warmup and drain valves always run in manual mode (Fortran lines 1097–1107)
            HP_warmup_VPo = VP_new(HP_warmup_VPo, HP_warmup_VPi, self.HP_warmup_vs.v, ts)
            HP_drain_VPo  = VP_new(HP_drain_VPo,  HP_drain_VPi,  self.HP_drain_vs.v,  ts)
            LP_warmup_VPo = VP_new(LP_warmup_VPo, LP_warmup_VPi, self.LP_warmup_vs.v, ts)
            LP_drain_VPo  = VP_new(LP_drain_VPo,  LP_drain_VPi,  self.LP_drain_vs.v,  ts)

        # Write valve position outputs (these are updated only at iteration 0 but always assigned)
        self.HPT_CV_VPo.v = HPT_CV_VPo
        self.HP_Bypass_VPo.v = HP_Bypass_VPo
        self.HP_Aux_VPo.v = HP_Aux_VPo
        self.HP_warmup_VPo.v = HP_warmup_VPo
        self.HP_drain_VPo.v = HP_drain_VPo
        self.Aux_DA_VPo.v = Aux_DA_VPo
        self.LP_Bypass_VPo.v = LP_Bypass_VPo
        self.LP_Aux_VPo.v = LP_Aux_VPo
        self.LP_warmup_VPo.v = LP_warmup_VPo
        self.LP_drain_VPo.v = LP_drain_VPo

        # -------------------------------------------------------------------------
        # INITIALISE SUB-TIMESTEP ACCUMULATORS
        # -------------------------------------------------------------------------
        m_dot_HPT_tot = 0.0
        W_dot_HPT_tot = 0.0
        m_dot_LPT_tot = 0.0
        W_dot_LPT_tot = 0.0
        m_dot_LPT_exh_tot = 0.0
        h_LPT_exh_tot = 0.0
        m_dot_HPT_exh_tot = 0.0
        m_dot_HP_bypass_tot = 0.0
        m_dot_HP_warmup_tot = 0.0
        m_dot_HP_drain_tot = 0.0
        m_dot_HP_AUX_tot = 0.0
        m_dot_AUX_ts_tot = 0.0
        m_dot_AUX_DA_tot = 0.0
        m_dot_LP_bypass_tot = 0.0
        m_dot_LP_warmup_tot = 0.0
        m_dot_LP_AUX_tot = 0.0
        m_dot_LP_drain_tot = 0.0
        m_dot_SS_drain_tot = 0.0
        m_dot_SGT_tot = 0.0
        h_LP_bypass_tot = 0.0
        h_LP_warmup_tot = 0.0
        h_HP_warmup_tot = 0.0
        h_AUX_DA_tot = 0.0
        h_SS_drain_tot = 0.0

        # Sub-timestep size (target 0.1 s, rounded so ts is evenly divisible)
        t_crit_base = 0.1
        ts_sub_i = math.ceil(ts / t_crit_base)
        t_crit = ts / float(ts_sub_i)

        # Working copies of dynamic piping state — will be updated each sub-step
        m_HPmain_prev = self._m_HPmain_prev
        P_HPmain_prev = self._P_HPmain_prev
        h_HPmain_prev = self._h_HPmain_prev
        T_HP_pipe = self._T_HP_pipe
        m_LPmain_prev = self._m_LPmain_prev
        P_LPmain_prev = self._P_LPmain_prev
        h_LPmain_prev = self._h_LPmain_prev
        T_LP_pipe = self._T_LP_pipe
        m_Aux_prev = self._m_Aux_prev
        P_Aux_prev = self._P_Aux_prev
        h_Aux_prev = self._h_Aux_prev
        T_Aux_pipe = self._T_Aux_pipe

        # Working copies of turbine state carried across sub-steps
        m_dot_HPT_in = self.m_dot_HPT_in.v
        m_dot_LPT_in = self.m_dot_LPT_in.v
        h_HPT1 = self.H_HPT1.v
        h_LPT1 = self.H_LPT1.v
        h_LPT2 = self.H_LPT2.v
        h_LPT3 = self.H_LPT3.v
        h_LPT4 = self.H_LPT4.v

        # Retrieved design-point values from static state
        h_HPT_in_d = self._h_HPT_in_d
        v_HPT_in_d = self._v_HPT_in_d
        s_HPT_in_d = self._s_HPT_in_d
        P_GS_out_d = self._P_GS_out_d
        h_GS_out_d = self._h_GS_out_d
        h_HPT1_d = self._h_HPT1_d
        h_LPT_in_d = self._h_LPT_in_d
        v_LPT_in_d = self._v_LPT_in_d
        s_LPT_in_d = self._s_LPT_in_d
        h_LPT1_d = self._h_LPT1_d
        h_LPT2_d = self._h_LPT2_d
        h_LPT3_d = self._h_LPT3_d
        h_LPT4_d = self._h_LPT4_d

        # Piping output friction factors carried across iterations
        ff_HP = self.ff_HPmain.v
        ff_LP = self.ff_LPmain.v
        ff_AUX = self.ff_AUX.v

        # Trip-to-turbine-seals state (initialise; will be set inside loop)
        TRIP_TS = 0.0

        # =========================================================================
        # SUB-TIMESTEP LOOP
        # =========================================================================
        for i in range(1, ts_sub_i + 1):

            # ------------------------------------------------------------------
            # i == 1: dynamic-array recall already done above (Python instance attrs
            # replace getDynamicArrayValueLastTimestep).  Nothing extra needed here.
            # ------------------------------------------------------------------

            # ------------------------------------------------------------------
            # TURBINE OFF BRANCH
            # ------------------------------------------------------------------
            if Turbine_ON != 1.0:
                m_dot_HPT_in = 0.0
                m_dot_HPT_exh = 0.0
                m_dot_LPT_in = 0.0
                m_dot_LPT_exh = 0.0
                h_LPT_exh = 0.0
                W_dot_LPT = 0.0
                W_dot_HPT = 0.0
                W_dot_total = 0.0
                H_HPT_in = 0.0
                H_LPT_in = 0.0
                h_HPT1 = 0.0
                H_HPT2 = 0.0
                P_HPT1 = 0.0
                P_HPT2 = 0.0
                P_GS_out = 0.0
                H_GS_out = 0.0
                m_dot_HPTS2 = 0.0
                m_dot_HPT_exh = 0.0
                P_HPT_exh = 0.0
                m_dot_LPTS1 = m_dot_LPT_in
                m_dot_LPTS2 = m_dot_LPT_in
                m_dot_LPTS3 = m_dot_LPT_in
                m_dot_LPTS4 = m_dot_LPT_in
                P_LPT1 = 0.0
                P_LPT2 = 0.0
                P_LPT3 = 0.0
                P_LPT4 = 0.0
                P_LPT_exh = 0.0
                T_LPT_in = 0.0
                # Turbine outputs not explicitly reset in Fortran OFF branch —
                # outputs persist from last turbine-ON iteration.

            else:
                # --------------------------------------------------------------
                # TURBINE ON: HIGH PRESSURE TURBINE MASS FLOW ITERATION
                # (replaces Fortran outer do-while around label 10, lines 1171–1361)
                # --------------------------------------------------------------
                pressure_tol_HPT = 1000.0
                pressure_error = pressure_tol_HPT + 1.0
                whileiterations1 = 0.0
                LR = 0.01
                m_dot_max = self.m_dot_HPT_d.v * 1.25
                m_dot_min = 0.0000001
                if m_dot_HPT_in == 0.0:
                    m_dot_HPT_in = 1.0
                SD_boundary_check = 0.0
                times_negative = 0.0
                times_positive = 0.0
                turbine_upper_bound = m_dot_max
                small_step = 0.000001
                enthalpy_ss = 10000.0
                enthalpy_LR = 0.1
                P_minimum_HPT = P_LPmain_prev - 1000.0
                m_dot_HPT_prev = m_dot_HPT_in
                pressure_error_prev = 0.0

                # Local turbine-state variables (initialised before entering loops)
                P_HPT1 = 0.0
                P_HPT2 = 0.0
                H_HPT2 = 0.0
                T_HPT2 = 0.0
                P_GS_out = 0.0
                H_GS_out = 0.0
                m_dot_HPTS2 = 0.0
                m_dot_HPT_exh = 0.0

                # OUTER pressure-convergence while (replaces Fortran outer do-while / label 10)
                while abs(pressure_error) >= pressure_tol_HPT or pressure_error < 0.0:

                    # -- Outer secant/bisection step (Fortran lines 1183–1256) --
                    if whileiterations1 > 1.0:
                        if pressure_error > 0.0:
                            times_negative = 0.0
                            times_positive += 1.0
                            if times_positive >= 50.0:
                                if times_positive == 50.0:
                                    if turbine_upper_bound < m_dot_max:
                                        m_dot_HPT_prev = m_dot_HPT_in
                                        pressure_error_prev = pressure_error
                                        m_dot_HPT_in = turbine_upper_bound
                                    else:
                                        times_positive = 0.0
                                        m_dot_HPT_in += small_step
                                else:  # times_positive > 50 — increase upper bound
                                    turbine_upper_bound = min(turbine_upper_bound * 1.001, m_dot_max)
                                    times_positive = 0.0
                                    m_dot_HPT_prev = m_dot_HPT_in
                                    pressure_error_prev = pressure_error  # Fortran uses turbine_error (typo), use pressure_error
                                    m_dot_HPT_in += small_step
                            elif whileiterations1 == 2.0:
                                m_dot_HPT_prev = m_dot_HPT_in
                                pressure_error_prev = pressure_error
                                m_dot_HPT_in = min(m_dot_HPT_in * 1.01, turbine_upper_bound)
                            else:
                                if pressure_error_prev < 0.0:
                                    m_dot_new = min(m_dot_HPT_in * 1.001, turbine_upper_bound)
                                    m_dot_HPT_in += (m_dot_new - m_dot_HPT_in) * LR
                                elif m_dot_HPT_in != m_dot_HPT_prev:
                                    flow_m = (pressure_error - pressure_error_prev) / (m_dot_HPT_in - m_dot_HPT_prev)
                                    flow_int = pressure_error - flow_m * m_dot_HPT_in
                                    if flow_m != 0.0:
                                        m_dot_new = max(min(-flow_int / flow_m, turbine_upper_bound), m_dot_min)
                                        delta_m = m_dot_new - m_dot_HPT_in
                                        if delta_m > 0.0:
                                            m_dot_new = m_dot_HPT_in + delta_m * LR
                                        else:
                                            m_dot_new = min(m_dot_HPT_in * 1.01, turbine_upper_bound)
                                    else:
                                        m_dot_new = min(m_dot_HPT_in * 1.01, turbine_upper_bound)
                                    m_dot_HPT_in += (m_dot_new - m_dot_HPT_in) * LR
                                else:
                                    m_dot_HPT_in = min(m_dot_HPT_in * 1.01, turbine_upper_bound)
                            if (turbine_upper_bound - m_dot_HPT_in) <= 0.0001:
                                m_dot_HPT_in = m_dot_HPT_prev + (turbine_upper_bound - m_dot_HPT_prev) * LR
                            m_dot_HPT_prev = m_dot_HPT_in
                            pressure_error_prev = pressure_error
                        else:  # pressure_error <= 0 — take step back
                            times_negative += 1.0
                            times_positive = 0.0
                            pressure_error_prev = pressure_error
                            if times_negative == 1.0:
                                turbine_upper_bound = m_dot_HPT_in
                                m_dot_HPT_prev = m_dot_HPT_in
                                m_dot_HPT_in *= 0.995
                            else:
                                m_dot_HPT_prev = m_dot_HPT_in
                                m_dot_HPT_in = max(m_dot_HPT_in - m_dot_HPT_in * 0.01, m_dot_min)

                    # -- INNER enthalpy-convergence while (Fortran lines 1258–1361) --
                    enthalpy_tol = 100.0
                    enthalpy_error = enthalpy_tol + 1.0
                    whileiterations2 = 0.0
                    goto10_flag = False

                    while abs(enthalpy_error) > enthalpy_tol:
                        whileiterations2 += 1.0

                        # Solve for HPT efficiency via Spencer-Cotton-Cannon
                        P_HPT_in = P_HPmain_prev
                        H_HPT_in = h_HPmain_prev
                        v_HPT_in = fp.spec_vol("water", P=P_HPT_in, h=H_HPT_in)  # CONVERTED-NEEDS UNITS CHECK: eeslib uses SI; inputs confirmed in Pa/J/kg
                        eta_HPT = eta_SCC_HPT(
                            m_dot_HPT_in, self.m_dot_HPT_d.v, self.HPT_parallel_sects.v,
                            self.GS_diameter.v, P_HPT_in, v_HPT_in, v_HPT_in_d,
                            self.P_HPT_exh_d.v, self.HPT_CV_NUMBER.v, self.HPT_no_GS.v
                        )

                        # Control-valve pressure drop
                        v_ref = 0.001  # specific volume of water at reference state [m^3/kg]
                        spec_grav = v_ref / v_HPT_in
                        CV_controlvalves = PB_CV_data(1, self.HPT_CV_D.v, HPT_CV_VPo)
                        m_dot_HPT_CV = m_dot_HPT_in / self.HPT_CV_NUMBER.v
                        HPT_Vol_in = m_dot_HPT_CV * v_HPT_in
                        HPT_Vol_in_gpm = HPT_Vol_in * 15850.3
                        DELTA_P_psi = spec_grav / (CV_controlvalves / HPT_Vol_in_gpm) ** 2.0
                        DELTA_P_CV = DELTA_P_psi * 6894.76
                        P_GS_out = P_HPT_in - DELTA_P_CV
                        H_GS_out = H_HPT_in  # no enthalpy change through control valves
                        s_GS_out = fp.entropy("water", P=P_GS_out, h=H_GS_out)  # CONVERTED-NEEDS UNITS CHECK: eeslib uses SI; inputs confirmed in Pa/J/kg

                        # HPT Stage 1 pressure via Stodola's Ellipse
                        P_HPT1 = StodolaStage(
                            P_GS_out_d, self.P_HPT1_d.v, P_GS_out,
                            h_GS_out_d, H_GS_out, self.m_dot_HPT_d.v, m_dot_HPT_in
                        )
                        if P_HPT1 <= P_minimum_HPT:
                            pressure_error = -10000.0
                            goto10_flag = True
                            break  # ← GO TO 10

                        # HPT extraction mass flows for stage 2
                        m_dot_HPTS2_d = self.m_dot_HPT_d.v - self.m_dot_HPT1_d.v
                        m_dot_HPTS2 = m_dot_HPT_in - m_dot_HPT1

                        # HPT Stage 2 pressure via Stodola's Ellipse
                        P_HPT2 = StodolaStage(
                            self.P_HPT1_d.v, self.P_HPT_exh_d.v, P_HPT1,
                            h_HPT1_d, h_HPT1, m_dot_HPTS2_d, m_dot_HPTS2
                        )
                        if P_HPT2 <= P_minimum_HPT:
                            pressure_error = -10000.0
                            goto10_flag = True
                            break  # ← GO TO 10

                        # Enthalpy and entropy leaving HPT stage 2
                        H_HPT2_s = fp.enthalpy("water", P=P_HPT2, s=s_GS_out)  # CONVERTED-NEEDS UNITS CHECK: eeslib uses SI; inputs confirmed in Pa/J/kg
                        H_HPT2 = H_GS_out - eta_HPT * (H_GS_out - H_HPT2_s)
                        s_HPT2 = fp.entropy("water", P=P_HPT2, h=H_HPT2)  # CONVERTED-NEEDS UNITS CHECK: eeslib uses SI; inputs confirmed in Pa/J/kg

                        # Expansion line slope and intercept
                        if s_HPT2 != s_GS_out:
                            m_elep = (H_HPT2 - H_GS_out) / (s_HPT2 - s_GS_out)
                        else:
                            m_elep = 1000000.0
                        b_elep = H_GS_out - s_GS_out * m_elep

                        # Current entropy at HPT stage 1 exit from previous h_HPT1
                        s_HPT1 = fp.entropy("water", P=P_HPT1, h=h_HPT1)  # CONVERTED-NEEDS UNITS CHECK: eeslib uses SI; inputs confirmed in Pa/J/kg
                        s_max = s_HPT2
                        s_min = s_GS_out
                        h_HPT1_prev = h_HPT1

                        # Nested entropy-iteration to find h_HPT1 on expansion line
                        LR2 = 0.25
                        whileiterations3 = 0.0
                        tol2 = 10.0
                        error2 = tol2 + 1.0
                        s_HPT1_prev_iter = s_HPT1
                        error2_prev = 0.0
                        while abs(error2) > tol2:
                            whileiterations3 += 1.0
                            h_elep = m_elep * s_HPT1 + b_elep
                            h_HPT1 = fp.enthalpy("water", P=P_HPT1, s=s_HPT1)  # CONVERTED-NEEDS UNITS CHECK: eeslib uses SI; inputs confirmed in Pa/J/kg
                            error2 = h_elep - h_HPT1
                            if abs(error2) > tol2:
                                if whileiterations3 > 1.0:
                                    if s_HPT1 != s_HPT1_prev_iter:
                                        m_guess = (error2 - error2_prev) / (s_HPT1 - s_HPT1_prev_iter)
                                        s_HPT1_prev_iter = s_HPT1
                                        error2_prev = error2
                                        if m_guess != 0.0:
                                            y_int = error2 - m_guess * s_HPT1
                                            s_new = max(min(-y_int / m_guess, s_max), s_min)
                                            s_HPT1 += (s_new - s_HPT1) * LR2
                                        else:
                                            if error2 > 0.0:
                                                s_HPT1 = min(s_HPT1 + 1.0, s_max)
                                            else:
                                                s_HPT1 = max(s_HPT1 - 1.0, s_min)
                                    else:
                                        s_HPT1_prev_iter = s_HPT1
                                        error2_prev = error2
                                        if error2 > 0.0:
                                            s_HPT1 = min(s_HPT1 + 1.0, s_max)
                                        else:
                                            s_HPT1 = max(s_HPT1 - 1.0, s_min)
                                else:
                                    s_HPT1_prev_iter = s_HPT1
                                    error2_prev = error2
                                    if error2 > 0.0:
                                        s_HPT1 = min(s_HPT1 + 1.0, s_max)
                                    else:
                                        s_HPT1 = max(s_HPT1 - 1.0, s_min)

                        # Update h_HPT1 with learning rate
                        enthalpy_error = h_HPT1 - h_HPT1_prev
                        if abs(enthalpy_error) <= enthalpy_tol:
                            h_HPT1 = min(h_HPT1 + enthalpy_error * enthalpy_LR, h_HPT1 + enthalpy_ss)
                        else:
                            h_HPT1 = max(h_HPT1 + enthalpy_error * enthalpy_LR, h_HPT1 - enthalpy_ss)
                        # end of inner enthalpy while

                    # Replaces GO TO 10: if a pressure-minimum violation was flagged,
                    # skip the pressure-error update and restart the outer while.
                    if goto10_flag:
                        whileiterations1 += 1.0
                        continue  # ← GO TO 10: restart outer pressure while

                    # After inner while converges, evaluate outer pressure error
                    pressure_error = P_HPT2 - P_LPmain_prev
                    whileiterations1 += 1.0
                    # end of HPT outer while

                P_HPT_exh = P_HPT2
                m_dot_HPT_exh = m_dot_HPTS2 - m_dot_HPT2  # subtract HP FWH2 extraction

                # Accumulate HPT mass-flow total (W_dot_HPT is accumulated in [E] after LPT while)
                m_dot_HPT_exh_tot += m_dot_HPT_exh * t_crit

                # --------------------------------------------------------------
                # LOW PRESSURE TURBINE MASS FLOW ITERATION
                # (replaces Fortran outer do-while around label 20, lines 1393–1635)
                # --------------------------------------------------------------
                P_LPT_in = P_LPmain_prev
                h_LPT_in = h_LPmain_prev
                pressure_tol_LPT = 2000.0
                pressure_error = pressure_tol_LPT + 1.0
                P_minimum_LPT = max(P_cond - 1000.0, 1000.0)
                turbine_upper_bound = m_dot_max
                if m_dot_LPT_in == 0.0:
                    m_dot_LPT_in = 1.0
                whileiterations1 = 0.0
                times_negative = 0.0
                times_positive = 0.0
                m_dot_LPT_prev = m_dot_LPT_in
                pressure_error_prev = 0.0

                # LPT stage design mass flows
                m_dot_LPTS1_d = self.m_dot_LPT_d.v
                m_dot_LPTS2_d = m_dot_LPTS1_d - self.m_dot_LPT1_d.v
                m_dot_LPTS3_d = m_dot_LPTS2_d - self.m_dot_LPT2_d.v
                m_dot_LPTS4_d = m_dot_LPTS3_d - self.m_dot_LPT3_d.v
                m_dot_LPT_exh_d = m_dot_LPTS4_d - self.m_dot_LPT4_d.v

                # LPT actual stage mass flows (based on current m_dot_LPT_in / extraction inputs)
                m_dot_LPTS1 = m_dot_LPT_in
                m_dot_LPTS2 = m_dot_LPT_in - m_dot_LPT1
                m_dot_LPTS3 = m_dot_LPTS2 - m_dot_LPT2
                m_dot_LPTS4 = m_dot_LPTS3 - m_dot_LPT3
                m_dot_LPT_exh = m_dot_LPTS4 - m_dot_LPT4

                # LPT stage exit pressures and enthalpies (initialise)
                P_LPT1 = 0.0
                P_LPT2 = 0.0
                P_LPT3 = 0.0
                P_LPT4 = 0.0
                P_LPT_exh = 0.0
                h_LPT_exh = 0.0
                T_LPT_in = 0.0
                v_LPT_in = 0.0
                s_LPT_in = 0.0

                # OUTER pressure-convergence while (replaces Fortran outer do-while / label 20)
                while abs(pressure_error) >= pressure_tol_LPT or pressure_error < 0.0:
                    # -- OUTER STEP: update m_dot_LPT_in via secant/bisection
                    #    (replaces Fortran lines 1404–1475; mirrors HPT outer step)
                    whileiterations1 += 1.0
                    if whileiterations1 > 1.0:
                        if pressure_error > 0.0:
                            times_negative = 0.0
                            times_positive = times_positive + 1.0
                            if times_positive >= 50.0:
                                if times_positive == 50.0:
                                    if turbine_upper_bound < m_dot_max:
                                        m_dot_LPT_prev = m_dot_LPT_in
                                        pressure_error_prev = pressure_error
                                        m_dot_LPT_in = turbine_upper_bound
                                    else:
                                        times_positive = 0.0
                                        m_dot_LPT_in = m_dot_LPT_in + small_step
                                else:
                                    turbine_upper_bound = min(turbine_upper_bound * 1.001, m_dot_max)
                                    times_positive = 0.0
                                    m_dot_LPT_prev = m_dot_LPT_in
                                    pressure_error_prev = pressure_error  # NOTE: Fortran uses turbine_error — typo; use pressure_error
                                    m_dot_LPT_in = m_dot_LPT_in + small_step
                            elif whileiterations1 == 2.0:
                                m_dot_LPT_prev = m_dot_LPT_in
                                pressure_error_prev = pressure_error
                                m_dot_LPT_in = min(m_dot_LPT_in * 1.01, turbine_upper_bound)
                            else:
                                if pressure_error_prev < 0.0:
                                    m_dot_new = min(m_dot_LPT_in * 1.001, turbine_upper_bound)
                                    m_dot_LPT_in = m_dot_LPT_in + (m_dot_new - m_dot_LPT_in) * LR
                                elif m_dot_LPT_in != m_dot_LPT_prev:
                                    flow_m = (pressure_error - pressure_error_prev) / (m_dot_LPT_in - m_dot_LPT_prev)
                                    flow_int = pressure_error - flow_m * m_dot_LPT_in
                                    if flow_m != 0.0:
                                        m_dot_new = max(min(-flow_int / flow_m, turbine_upper_bound), m_dot_min)
                                        delta_m = m_dot_new - m_dot_LPT_in
                                        if delta_m > 0.0:
                                            m_dot_new = m_dot_LPT_in + delta_m * LR
                                        else:
                                            m_dot_new = min(m_dot_LPT_in * 1.01, turbine_upper_bound)
                                    else:
                                        m_dot_new = min(m_dot_LPT_in * 1.01, turbine_upper_bound)
                                    m_dot_LPT_in = m_dot_LPT_in + (m_dot_new - m_dot_LPT_in) * LR
                                else:
                                    m_dot_LPT_in = min(m_dot_LPT_in * 1.01, turbine_upper_bound)
                            if (turbine_upper_bound - m_dot_LPT_in) <= 0.0001:
                                m_dot_LPT_in = m_dot_LPT_prev + (turbine_upper_bound - m_dot_LPT_prev) * LR
                            m_dot_LPT_prev = m_dot_LPT_in
                            pressure_error_prev = pressure_error
                        else:  # pressure_error <= 0: step back
                            times_negative = times_negative + 1.0
                            times_positive = 0.0
                            pressure_error_prev = pressure_error
                            if times_negative == 1.0:
                                turbine_upper_bound = m_dot_LPT_in
                                m_dot_LPT_prev = m_dot_LPT_in
                                m_dot_LPT_in = m_dot_LPT_in * 0.995
                            else:
                                m_dot_LPT_prev = m_dot_LPT_in
                                m_dot_LPT_in = max(m_dot_LPT_in - m_dot_LPT_in * 0.01, m_dot_min)

                    # Update actual stage mass flows with new m_dot_LPT_in
                    m_dot_LPTS1 = m_dot_LPT_in
                    m_dot_LPTS2 = m_dot_LPT_in - m_dot_LPT1
                    m_dot_LPTS3 = m_dot_LPTS2 - m_dot_LPT2
                    m_dot_LPTS4 = m_dot_LPTS3 - m_dot_LPT3
                    m_dot_LPT_exh = m_dot_LPTS4 - m_dot_LPT4

                    # FIT_PH for LPT inlet (Fortran lines 1481–1483)
                    v_LPT_in = fp.spec_vol("water", P=P_LPT_in, h=h_LPT_in)       # CONVERTED-NEEDS UNITS CHECK: eeslib uses SI; inputs confirmed in Pa/J/kg
                    s_LPT_in = fp.entropy("water", P=P_LPT_in, h=h_LPT_in)        # CONVERTED-NEEDS UNITS CHECK: eeslib uses SI; inputs confirmed in Pa/J/kg
                    T_LPT_in = fp.temperature("water", P=P_LPT_in, h=h_LPT_in)    # CONVERTED-NEEDS UNITS CHECK: eeslib uses SI; inputs confirmed in Pa/J/kg

                    # -- INNER ENTHALPY-CONVERGENCE WHILE (Fortran inner do-while, lines 1489–1636) --
                    # GO TO 20 inside body: set goto20_flag=True; break → outer while continues
                    enthalpy_tol = 100.0
                    enthalpy_error = enthalpy_tol + 1.0
                    enthalpy_error_prev = enthalpy_tol + 1.0  # init so 2x check passes first time
                    whileiterations2 = 0.0
                    goto20_flag = False

                    while abs(enthalpy_error) > enthalpy_tol:
                        whileiterations2 += 1.0

                        # LPT efficiency
                        eta_LPT = eta_SCC_LPT(
                            m_dot_LPT_in, self.m_dot_LPT_d.v,
                            self.LPT_parallel_sects.v,
                            P_LPT_in, T_LPT_in, v_LPT_in, self._v_LPT_in_d,
                            self.P_cond_d.v
                        )

                        # Stage exit pressures via Stodola
                        P_LPT1 = StodolaStage(
                            self._P_LPT_in_d, self._P_LPT1_d,
                            P_LPT_in, self._h_LPT_in_d, h_LPT_in,
                            m_dot_LPTS1_d, m_dot_LPTS1
                        )
                        if P_LPT1 <= P_minimum_LPT:
                            pressure_error = -10000.0
                            goto20_flag = True
                            break

                        P_LPT2 = StodolaStage(
                            self._P_LPT1_d, self._P_LPT2_d,
                            P_LPT1, self._h_LPT1_d, h_LPT1,
                            m_dot_LPTS2_d, m_dot_LPTS2
                        )
                        if P_LPT2 <= P_minimum_LPT:
                            pressure_error = -10000.0
                            goto20_flag = True
                            break

                        P_LPT3 = StodolaStage(
                            self._P_LPT2_d, self._P_LPT3_d,
                            P_LPT2, self._h_LPT2_d, h_LPT2,
                            m_dot_LPTS3_d, m_dot_LPTS3
                        )
                        if P_LPT3 <= P_minimum_LPT:
                            pressure_error = -10000.0
                            goto20_flag = True
                            break

                        P_LPT4 = StodolaStage(
                            self._P_LPT3_d, self._P_LPT4_d,
                            P_LPT3, self._h_LPT3_d, h_LPT3,
                            m_dot_LPTS4_d, m_dot_LPTS4
                        )
                        if P_LPT4 <= P_minimum_LPT:
                            pressure_error = -10000.0
                            goto20_flag = True
                            break

                        P_LPT_exh = StodolaStage(
                            self._P_LPT4_d, self.P_cond_d.v,
                            P_LPT4, self._h_LPT4_d, h_LPT4,
                            m_dot_LPT_exh_d, m_dot_LPT_exh
                        )
                        if P_LPT_exh <= P_minimum_LPT:
                            pressure_error = -10000.0
                            goto20_flag = True
                            break

                        # Pressure error set INSIDE inner while for LPT (difference from HPT)
                        pressure_error = P_LPT_exh - P_cond

                        # Isentropic exhaust enthalpy + actual exhaust enthalpy
                        s_LPT_exh_s = s_LPT_in
                        h_LPT_exh_s = fp.enthalpy("water", P=P_LPT_exh, s=s_LPT_exh_s)  # CONVERTED-NEEDS UNITS CHECK: eeslib uses SI; inputs confirmed in Pa/J/kg
                        h_LPT_exh = h_LPT_in - eta_LPT * (h_LPT_in - h_LPT_exh_s)
                        s_LPT_exh = fp.entropy("water", P=P_LPT_exh, h=h_LPT_exh)       # CONVERTED-NEEDS UNITS CHECK: eeslib uses SI; inputs confirmed in Pa/J/kg
                        s_prime = self.LPT_EXP_A0.v + self.LPT_EXP_A1.v * h_LPT_exh + self.LPT_EXP_A2.v * h_LPT_exh ** 2.0
                        DELTA_S = s_LPT_exh - s_prime

                        # Stage enthalpy updates via h_lpt_stage
                        h_LPT1_prev = h_LPT1
                        h_LPT1_new = h_lpt_stage(
                            h_LPT1, h_LPT_exh, h_LPT_in, P_LPT1,
                            DELTA_S, self.LPT_EXP_A0.v, self.LPT_EXP_A1.v, self.LPT_EXP_A2.v, 0.0005
                        )
                        enthalpy_error_1 = h_LPT1_new - h_LPT1_prev
                        if abs(enthalpy_error_1) <= enthalpy_tol:
                            h_LPT1 = min(h_LPT1 + enthalpy_error_1 * enthalpy_LR, h_LPT1 + enthalpy_ss)
                        else:
                            h_LPT1 = max(h_LPT1 + enthalpy_error_1 * enthalpy_LR, h_LPT1 - enthalpy_ss)

                        h_LPT2_prev = h_LPT2
                        h_LPT2_new = h_lpt_stage(
                            h_LPT2, h_LPT_exh, h_LPT_in, P_LPT2,
                            DELTA_S, self.LPT_EXP_A0.v, self.LPT_EXP_A1.v, self.LPT_EXP_A2.v, 0.0005
                        )
                        enthalpy_error_2 = h_LPT2_new - h_LPT2_prev
                        if enthalpy_error_2 > 0.0:
                            h_LPT2 = min(h_LPT2 + enthalpy_error_2 * enthalpy_LR, h_LPT2 + enthalpy_ss)
                        else:
                            h_LPT2 = max(h_LPT2 + enthalpy_error_2 * enthalpy_LR, h_LPT2 - enthalpy_ss)

                        h_LPT3_prev = h_LPT3
                        h_LPT3_new = h_lpt_stage(
                            h_LPT3, h_LPT_exh, h_LPT_in, P_LPT3,
                            DELTA_S, self.LPT_EXP_A0.v, self.LPT_EXP_A1.v, self.LPT_EXP_A2.v, 0.0005
                        )
                        enthalpy_error_3 = h_LPT3_new - h_LPT3_prev
                        if enthalpy_error_3 > 0.0:
                            h_LPT3 = min(h_LPT3 + enthalpy_error_3 * enthalpy_LR, h_LPT3 + enthalpy_ss)
                        else:
                            h_LPT3 = max(h_LPT3 + enthalpy_error_3 * enthalpy_LR, h_LPT3 - enthalpy_ss)

                        h_LPT4_prev = h_LPT4
                        h_LPT4_new = h_lpt_stage(
                            h_LPT4, h_LPT_exh, h_LPT_in, P_LPT4,
                            DELTA_S, self.LPT_EXP_A0.v, self.LPT_EXP_A1.v, self.LPT_EXP_A2.v, 0.0005
                        )
                        enthalpy_error_4 = h_LPT4_new - h_LPT4_prev
                        if enthalpy_error_4 > 0.0:
                            h_LPT4 = min(h_LPT4 + enthalpy_error_4 * enthalpy_LR, h_LPT4 + enthalpy_ss)
                        else:
                            h_LPT4 = max(h_LPT4 + enthalpy_error_4 * enthalpy_LR, h_LPT4 - enthalpy_ss)

                        # Combined enthalpy error; dampen if diverging
                        enthalpy_error_prev = enthalpy_error
                        enthalpy_error = (abs(enthalpy_error_4) + abs(enthalpy_error_3)
                                          + abs(enthalpy_error_2) + abs(enthalpy_error_1))
                        if enthalpy_error > 2.0 * enthalpy_error_prev:
                            enthalpy_error = enthalpy_tol / 2.0

                        # On convergence, revert stage enthalpies to previous (pre-update) values
                        if abs(enthalpy_error) < enthalpy_tol:
                            h_LPT4 = h_LPT4_prev
                            h_LPT3 = h_LPT3_prev
                            h_LPT2 = h_LPT2_prev
                            h_LPT1 = h_LPT1_prev
                        # end inner LPT while

                    # Replaces GO TO 20: restart outer while if any P_LPTx undershot
                    if goto20_flag:
                        continue  # ← restart outer pressure while

                    # pressure_error already set inside inner while (P_LPT_exh - P_cond)
                    # end of LPT outer while

                # ---- Turbine work and stage temperatures (Fortran lines 1646–1680) ----
                T_HPT1 = fp.temperature("water", P=P_HPT1, h=h_HPT1)    # CONVERTED-NEEDS UNITS CHECK: eeslib uses SI; inputs confirmed in Pa/J/kg
                T_HPT2 = fp.temperature("water", P=P_HPT2, h=H_HPT2)    # CONVERTED-NEEDS UNITS CHECK: eeslib uses SI; inputs confirmed in Pa/J/kg
                T_LPT1 = fp.temperature("water", P=P_LPT1, h=h_LPT1)    # CONVERTED-NEEDS UNITS CHECK: eeslib uses SI; inputs confirmed in Pa/J/kg
                T_LPT2 = fp.temperature("water", P=P_LPT2, h=h_LPT2)    # CONVERTED-NEEDS UNITS CHECK: eeslib uses SI; inputs confirmed in Pa/J/kg
                T_LPT3 = fp.temperature("water", P=P_LPT3, h=h_LPT3)    # CONVERTED-NEEDS UNITS CHECK: eeslib uses SI; inputs confirmed in Pa/J/kg
                T_LPT4 = fp.temperature("water", P=P_LPT4, h=h_LPT4)    # CONVERTED-NEEDS UNITS CHECK: eeslib uses SI; inputs confirmed in Pa/J/kg
                T_LPT_exh = fp.temperature("water", P=P_LPT_exh, h=h_LPT_exh)  # CONVERTED-NEEDS UNITS CHECK: eeslib uses SI; inputs confirmed in Pa/J/kg

                W_dot_HPT = (m_dot_HPT_in * (H_HPT_in - h_HPT1)
                             + m_dot_HPTS2 * (h_HPT1 - H_HPT2))
                W_dot_LPT = (m_dot_LPTS1 * (h_LPT_in - h_LPT1)
                             + m_dot_LPTS2 * (h_LPT1 - h_LPT2)
                             + m_dot_LPTS3 * (h_LPT2 - h_LPT3))
                W_dot_LPT += (m_dot_LPTS4 * (h_LPT3 - h_LPT4)
                              + m_dot_LPTS4 * (h_LPT4 - h_LPT_exh))  # NOTE: Fortran bug uses m_dot_LPTS4 twice
                W_dot_total = W_dot_HPT + W_dot_LPT

                # Accumulate turbine totals (Fortran sub-timestep accumulators)
                W_dot_HPT_tot += W_dot_HPT * t_crit
                m_dot_LPT_exh_tot += m_dot_LPT_exh * t_crit
                h_LPT_exh_tot += m_dot_LPT_exh * h_LPT_exh * t_crit
                W_dot_LPT_tot += W_dot_LPT * t_crit

                # Write per-iteration turbine outputs (Fortran lines 1682–1725)
                self.W_dot_tot.v = W_dot_total
                self.P_GS_out.v = P_GS_out
                self.h_GS_out.v = H_GS_out
                self.m_dot_HPT_in.v = m_dot_HPT_in
                self.m_dot_HPTS1.v = m_dot_HPT_in
                self.P_HPT1.v = P_HPT1
                self.H_HPT1.v = h_HPT1
                self.m_dot_HPTS2.v = m_dot_HPTS2
                self.P_HPT2.v = P_HPT2
                self.H_HPT2.v = H_HPT2
                self.m_dot_HPT_exh.v = m_dot_HPT_exh
                self.P_HPT_exh.v = P_HPT2   # NOTE: Fortran output 32 = P_HPT2 (not P_HPT_exh_d)
                self.H_HPT_exh.v = H_HPT2   # NOTE: Fortran output 33 = H_HPT2
                self.T_HPT_exh.v = T_HPT2
                self.m_dot_LPT_in.v = m_dot_LPT_in
                self.P_LPT1.v = P_LPT1
                self.T_LPT1.v = T_LPT1
                self.H_LPT1.v = h_LPT1
                self.m_dot_LPTS2.v = m_dot_LPTS2
                self.P_LPT2.v = P_LPT2
                self.T_LPT2.v = T_LPT2
                self.H_LPT2.v = h_LPT2
                self.m_dot_LPTS3.v = m_dot_LPTS3
                self.P_LPT3.v = P_LPT3
                self.T_LPT3.v = T_LPT3
                self.H_LPT3.v = h_LPT3
                self.m_dot_LPTS4.v = m_dot_LPTS4
                self.T_LPT4.v = T_LPT4
                self.P_LPT4.v = P_LPT4
                self.H_LPT4.v = h_LPT4
                self.m_dot_LPT_exh.v = m_dot_LPT_exh
                self.Vol_dot_LPT_exh.v = 0.0
                self.T_LPT_exh.v = T_LPT_exh
                self.H_LPT_exh.v = h_LPT_exh

            # end Turbine ON/OFF branch

            # ------------------------------------------------------------------
            # BYPASS & DRAIN VALVE FLOWS + ENTHALPIES
            # (replaces Fortran lines 1731–1853)
            # ------------------------------------------------------------------
            # HP bypass (Fortran lines 1706–1712)
            if self.HP_bypass_VPo.v != 0.0:
                HP_bypass_CV = PB_CV_data(self.HP_bypass_vt.v, self.HP_bypass_d.v, self.HP_bypass_VPo.v)
                m_dot_HP_bypass = valve_massflow(HP_bypass_CV, P_SGT_in, h_SGT_in, P_LPmain_prev)
                h_HP_bypass = h_SGT_in
            else:
                m_dot_HP_bypass = 0.0
                h_HP_bypass = h_SGT_in

            # HP AUX (Fortran lines 1714–1720)
            if self.HP_AUX_VPo.v != 0.0:
                HP_Aux_CV = PB_CV_data(self.HP_AUX_vt.v, self.HP_AUX_d.v, self.HP_AUX_VPo.v)
                m_dot_HP_AUX = valve_massflow(HP_Aux_CV, P_SGT_in, h_SGT_in, P_Aux_prev)
                h_HP_AUX = h_SGT_in
            else:
                m_dot_HP_AUX = 0.0
                h_HP_AUX = h_SGT_in

            # HP warmup (Fortran lines 1722–1731)
            if self.HP_Warmup_VPo.v != 0.0:
                HP_warmup_CV = PB_CV_data(self.HP_Warmup_vt.v, self.HP_Warmup_d.v, self.HP_Warmup_VPo.v)
                m_dot_HP_warmup = valve_massflow(HP_warmup_CV, P_HPmain_prev, h_HPmain_prev, P_cond)
                h_HP_warmup = h_HPmain_prev
            else:
                m_dot_HP_warmup = 0.0
                h_HP_warmup = h_HPmain_prev

            # HP drain — two-phase check (Fortran lines 1733–1754)
            if self.HP_drain_VPo.v != 0.0:
                HP_drain_CV = PB_CV_data(self.HP_drain_vt.v, self.HP_drain_d.v, self.HP_drain_VPo.v)
                h_sat_g_HPdrain = fp.enthalpy("water", P=P_HPmain_prev, Q=1.0)   # CONVERTED-NEEDS UNITS CHECK: eeslib uses SI; inputs confirmed in Pa/J/kg
                if h_HPmain_prev >= h_sat_g_HPdrain:
                    m_dot_HP_drain = 0.0
                    h_HP_drain = h_HPmain_prev
                else:
                    x_HPdrain = fp.quality("water", P=P_HPmain_prev, h=h_HPmain_prev)  # CONVERTED-NEEDS UNITS CHECK: eeslib uses SI; inputs confirmed in Pa/J/kg
                    h_sat_f_HPdrain = fp.enthalpy("water", P=P_HPmain_prev, Q=0.0)     # CONVERTED-NEEDS UNITS CHECK: eeslib uses SI; inputs confirmed in Pa/J/kg
                    m_dot_HP_drain_max = m_HPmain_prev * (1.0 - x_HPdrain)
                    m_dot_HP_drain = min(
                        valve_massflow(HP_drain_CV, P_HPmain_prev, h_sat_f_HPdrain, P_atm),
                        m_dot_HP_drain_max
                    )
                    h_HP_drain = h_sat_f_HPdrain
            else:
                m_dot_HP_drain = 0.0
                h_HP_drain = 0.0

            # AUX to DA (Fortran lines 1756–1764)
            if self.AUX_DA_VPo.v != 0.0:
                Aux_DA_CV = PB_CV_data(self.AUX_DA_vt.v, self.AUX_DA_d.v, self.AUX_DA_VPo.v)
                m_dot_AUX_DA = valve_massflow(Aux_DA_CV, P_Aux_prev, h_Aux_prev, P_DA)
                h_AUX_DA = h_Aux_prev
            else:
                m_dot_AUX_DA = 0.0
                h_AUX_DA = 0.0

            # LP bypass (Fortran lines 1768–1775)
            if self.LP_Bypass_VPo.v != 0.0:
                LP_bypass_CV = PB_CV_data(self.LP_bypass_vt.v, self.LP_bypass_d.v, self.LP_Bypass_VPo.v)
                m_dot_LP_bypass = valve_massflow(LP_bypass_CV, P_LPmain_prev, h_LPmain_prev, P_cond)
                h_LP_bypass = h_LPmain_prev
            else:
                m_dot_LP_bypass = 0.0
                h_LP_bypass = h_LPmain_prev

            # LP AUX (Fortran lines 1777–1784)
            if self.LP_AUX_VPo.v != 0.0:
                LP_Aux_CV = PB_CV_data(self.LP_AUX_vt.v, self.LP_AUX_d.v, self.LP_AUX_VPo.v)
                m_dot_LP_AUX = valve_massflow(LP_Aux_CV, P_LPmain_prev, h_LPmain_prev, P_Aux_prev)
                h_LP_AUX = h_LPmain_prev
            else:
                m_dot_LP_AUX = 0.0
                h_LP_AUX = h_LPmain_prev

            # LP warmup (Fortran lines 1786–1794)
            if self.LP_Warmup_VPo.v != 0.0:
                LP_warmup_CV = PB_CV_data(self.LP_Warmup_vt.v, self.LP_Warmup_d.v, self.LP_Warmup_VPo.v)
                m_dot_LP_warmup = valve_massflow(LP_warmup_CV, P_LPmain_prev, h_LPmain_prev, P_cond)
                h_LP_warmup = h_LPmain_prev
            else:
                m_dot_LP_warmup = 0.0
                h_LP_warmup = h_LPmain_prev

            # LP drain — two-phase check (Fortran lines 1796–1817)
            if self.LP_Drain_VPo.v != 0.0:
                LP_drain_CV = PB_CV_data(self.LP_drain_vt.v, self.LP_drain_d.v, self.LP_Drain_VPo.v)
                h_sat_g_LPdrain = fp.enthalpy("water", P=P_LPmain_prev, Q=1.0)   # CONVERTED-NEEDS UNITS CHECK: eeslib uses SI; inputs confirmed in Pa/J/kg
                if h_LPmain_prev >= h_sat_g_LPdrain:
                    m_dot_LP_drain = 0.0
                    h_LP_drain = h_LPmain_prev
                else:
                    x_LPdrain = fp.quality("water", P=P_LPmain_prev, h=h_LPmain_prev)  # CONVERTED-NEEDS UNITS CHECK: eeslib uses SI; inputs confirmed in Pa/J/kg
                    h_sat_f_LPdrain = fp.enthalpy("water", P=P_LPmain_prev, Q=0.0)     # CONVERTED-NEEDS UNITS CHECK: eeslib uses SI; inputs confirmed in Pa/J/kg
                    m_dot_LP_drain_max = m_LPmain_prev * (1.0 - x_LPdrain)
                    m_dot_LP_drain = min(
                        valve_massflow(LP_drain_CV, P_LPmain_prev, h_sat_f_LPdrain, P_atm),
                        m_dot_LP_drain_max
                    )
                    h_LP_drain = h_sat_f_LPdrain
            else:
                m_dot_LP_drain = 0.0
                h_LP_drain = 0.0

            # Turbine seal flow (Fortran lines 1819–1829)
            if self.Turbine_On.v != 1.0:
                if P_Aux_prev > self.P_ts_req.v:
                    m_dot_ts = self.m_dot_turbine_seals.v
                    TRIP_TS = 0.0  # flow going to turbine seals
                else:
                    m_dot_ts = max(
                        self.m_dot_turbine_seals.v * (P_Aux_prev - P_ts_min) / (self.P_ts_req.v - P_ts_min),
                        0.0
                    )
                    TRIP_TS = 1.0  # no flow to turbine seals
            else:
                m_dot_ts = 0.0
                TRIP_TS = 0.0

            # Accumulate valve-flow totals
            m_dot_HP_bypass_tot += m_dot_HP_bypass * t_crit
            m_dot_LP_bypass_tot += m_dot_LP_bypass * t_crit
            m_dot_HP_AUX_tot += m_dot_HP_AUX * t_crit
            m_dot_LP_AUX_tot += m_dot_LP_AUX * t_crit
            m_dot_HP_drain_tot += m_dot_HP_drain * t_crit
            m_dot_LP_drain_tot += m_dot_LP_drain * t_crit
            m_dot_HP_warmup_tot += m_dot_HP_warmup * t_crit
            m_dot_LP_warmup_tot += m_dot_LP_warmup * t_crit
            m_dot_AUX_DA_tot += m_dot_AUX_DA * t_crit
            m_dot_AUX_ts_tot += m_dot_ts * t_crit
            h_HP_warmup_tot += m_dot_HP_warmup * t_crit * h_HPmain_prev
            h_AUX_DA_tot += m_dot_AUX_DA * t_crit * h_Aux_prev
            h_LP_bypass_tot += m_dot_LP_bypass * t_crit * h_LPmain_prev
            h_LP_warmup_tot += m_dot_LP_warmup * t_crit * h_LPmain_prev

            # Write per-iteration bypass/drain outputs
            self.m_dot_HP_bypass.v = m_dot_HP_bypass
            self.m_dot_HP_AUX.v = m_dot_HP_AUX
            self.m_dot_HP_drain.v = m_dot_HP_drain
            self.m_dot_HP_warmup.v = m_dot_HP_warmup
            self.m_dot_LP_AUX.v = m_dot_LP_AUX
            self.m_dot_LP_bypass.v = m_dot_LP_bypass
            self.m_dot_LP_drain.v = m_dot_LP_drain
            self.m_dot_LP_warmup.v = m_dot_LP_warmup
            self.m_dot_AUX_DA.v = m_dot_AUX_DA

            # ------------------------------------------------------------------
            # STEAM SEPARATOR
            # (replaces Fortran lines 1855–1935)
            # ------------------------------------------------------------------
            m_dot_SS_in = m_dot_HPT_exh + m_dot_HP_bypass
            if m_dot_SS_in > 0.0:
                # Weighted inlet enthalpy (Fortran line 1851)
                h_SS_in = (m_dot_HPT_exh * H_HPT2 + m_dot_HP_bypass * h_HP_bypass) / m_dot_SS_in
                # Flash the mixture at LP-main pressure
                x_SS_in   = fp.quality("water", P=P_LPmain_prev, h=h_SS_in)    # CONVERTED-NEEDS UNITS CHECK: eeslib uses SI; inputs confirmed in Pa/J/kg
                T_SS_in   = fp.temperature("water", P=P_LPmain_prev, h=h_SS_in) # CONVERTED-NEEDS UNITS CHECK: eeslib uses SI; inputs confirmed in Pa/J/kg
                rho_SS_in = fp.density("water", P=P_LPmain_prev, h=h_SS_in)     # CONVERTED-NEEDS UNITS CHECK: eeslib uses SI; inputs confirmed in Pa/J/kg
                # Saturated vapour properties
                h_sat_g   = fp.enthalpy("water", P=P_LPmain_prev, Q=1.0)        # CONVERTED-NEEDS UNITS CHECK: eeslib uses SI; inputs confirmed in Pa/J/kg
                rho_sat_g = fp.density("water", P=P_LPmain_prev, Q=1.0)         # CONVERTED-NEEDS UNITS CHECK: eeslib uses SI; inputs confirmed in Pa/J/kg
                T_sat     = fp.temperature("water", P=P_LPmain_prev, Q=1.0)     # CONVERTED-NEEDS UNITS CHECK: eeslib uses SI; inputs confirmed in Pa/J/kg

                if h_SS_in >= h_sat_g:  # all steam
                    m_dot_SS_drain    = 0.0
                    Vol_dot_SS_drain  = 0.0
                    h_SS_drain        = 0.0
                    T_SS_drain        = 0.0
                    m_dot_SS_steam    = m_dot_SS_in
                    Vol_dot_SS_steam  = m_dot_SS_in / rho_SS_in
                    h_SS_steam        = h_SS_in
                    T_SS_steam        = T_SS_in
                else:  # two-phase: separate drain (liquid) from steam
                    # Round quality to 3 dp to aid convergence (Fortran: nint(x*1000)/1000)
                    x_SS_in = round(x_SS_in, 3)
                    h_sat_f   = fp.enthalpy("water", P=P_LPmain_prev, Q=0.0)    # CONVERTED-NEEDS UNITS CHECK: eeslib uses SI; inputs confirmed in Pa/J/kg
                    rho_sat_f = fp.density("water", P=P_LPmain_prev, Q=0.0)     # CONVERTED-NEEDS UNITS CHECK: eeslib uses SI; inputs confirmed in Pa/J/kg
                    # Drain (liquid) outputs
                    m_dot_SS_drain   = m_dot_SS_in * (1.0 - x_SS_in)
                    Vol_dot_SS_drain = m_dot_SS_drain / rho_sat_f
                    h_SS_drain       = h_sat_f
                    T_SS_drain       = T_sat
                    # Steam outputs
                    m_dot_SS_steam   = m_dot_SS_in - m_dot_SS_drain
                    Vol_dot_SS_steam = m_dot_SS_steam / rho_sat_g
                    h_SS_steam       = h_sat_g
                    T_SS_steam       = T_sat
            else:
                m_dot_SS_drain = 0.0
                Vol_dot_SS_drain = 0.0
                h_SS_drain = 0.0
                T_SS_drain = 0.0
                m_dot_SS_steam = 0.0
                Vol_dot_SS_steam = 0.0
                h_SS_steam = 0.0
                T_SS_steam = 0.0

            m_dot_SS_drain_tot += m_dot_SS_drain * t_crit
            h_SS_drain_tot += m_dot_SS_drain * t_crit * h_SS_drain

            # Write per-iteration steam separator outputs
            self.m_dot_SS_drain.v = m_dot_SS_drain
            self.Vol_dot_SS_drain.v = Vol_dot_SS_drain
            self.P_SS_drain.v = P_LPmain_prev
            self.h_SS_drain.v = h_SS_drain
            self.T_SS_drain.v = T_SS_drain
            self.m_dot_SS_steam.v = m_dot_SS_steam
            self.Vol_dot_SS_steam.v = Vol_dot_SS_steam
            self.P_SS_steam.v = P_LPmain_prev
            self.h_SS_steam.v = h_SS_steam
            self.T_SS_steam.v = T_SS_steam

            # ------------------------------------------------------------------
            # REHEATER (NTU-effectiveness)
            # (replaces Fortran lines 1937–2020)
            # ------------------------------------------------------------------
            h_HX_out = h_SS_steam
            T_HX_out = T_SS_steam
            Q_dot_HX = 0.0
            Eta_OD = 0.0
            T_HTF_out = HTF_T_in
            Vol_dot_HTF = 0.0

            if m_dot_SS_steam > 0.0:
                if m_dot_HTF > 0.0:
                    # HTF specific heat
                    cp_HTF = specheat(fnumd=self.Fluid_ID.v, T=HTF_T_in, P=HTF_P_in)
                    # Steam average temperature for cp estimate
                    T_avg_steam = (HTF_T_in + T_SS_steam) / 2.0
                    cp_steam = f_cp_water(P=P_LPmain_prev, T=T_avg_steam)
                    # Heat exchanger geometry: surface area of steam/water side
                    A_s = (math.pi * (self.HX_tube_OD.v - 2.0 * self.HX_tube_th.v)
                           * self.HX_length.v * self.HX_no_shell.v * self.HX_No_tubes.v)
                    UA_rated = self.HX_UA_d.v * A_s
                    UA_OD = UA_rated * (m_dot_HTF / self.m_dot_HTF_HX_d.v) ** self.HX_exp.v
                    # Capacity rates and NTU
                    Cap_steam = m_dot_SS_steam * cp_steam
                    Cap_HTF   = m_dot_HTF * cp_HTF
                    Cap_min   = min(Cap_steam, Cap_HTF)
                    CR        = max(Cap_min / max(Cap_steam, Cap_HTF), 0.001)
                    NTU_OD    = UA_OD / Cap_min
                    # One-pass shell-and-tube effectiveness
                    _sqrt_term = math.sqrt(1.0 + CR ** 2.0)
                    _exp_term  = math.exp(-NTU_OD * _sqrt_term)
                    eta_1pass  = (2.0 * (1.0 + CR + _sqrt_term
                                         * (1.0 + _exp_term) / (1.0 - _exp_term)) ** (-1.0))
                    # Multi-shell effectiveness
                    _ratio = (1.0 - eta_1pass * CR) / (1.0 - eta_1pass)
                    Eta_OD = (_ratio ** self.HX_no_shell.v - 1.0) / (_ratio ** self.HX_no_shell.v - CR)
                    # Heat transfer — clipped by temperature direction
                    if T_SS_steam < HTF_T_in:  # normal: HTF heats steam
                        h_HX_out_s = fp.enthalpy("water", P=P_LPmain_prev, T=HTF_T_in)  # CONVERTED-NEEDS UNITS CHECK: eeslib uses SI; inputs confirmed in Pa/J/kg
                        Q_dot_HX = min(
                            Eta_OD * m_dot_SS_steam * (h_HX_out_s - h_SS_steam),
                            Eta_OD * m_dot_HTF * cp_HTF * (HTF_T_in - T_SS_steam)
                        )
                        h_HX_out  = (m_dot_SS_steam * h_SS_steam + Q_dot_HX) / m_dot_SS_steam
                        T_HTF_out = (m_dot_HTF * cp_HTF * HTF_T_in - Q_dot_HX) / (m_dot_HTF * cp_HTF)
                    else:  # heat going wrong direction
                        h_HX_out_s = fp.enthalpy("water", P=P_LPmain_prev, T=HTF_T_in)  # CONVERTED-NEEDS UNITS CHECK: eeslib uses SI; inputs confirmed in Pa/J/kg
                        Q_dot_HX = min(
                            Eta_OD * m_dot_SS_steam * (h_SS_steam - h_HX_out_s),
                            Eta_OD * m_dot_HTF * cp_HTF * (T_SS_steam - HTF_T_in)
                        )
                        h_HX_out  = (m_dot_SS_steam * h_SS_steam - Q_dot_HX) / m_dot_SS_steam
                        T_HTF_out = (m_dot_HTF * cp_HTF * HTF_T_in + Q_dot_HX) / (m_dot_HTF * cp_HTF)
                    # HTF volumetric flow
                    rho_htf    = density(fnumd=self.Fluid_ID.v, T=T_HTF_out, P=HTF_P_in)
                    Vol_dot_HTF = m_dot_HTF / rho_htf
                    # Update steam exit temperature and volumetric flow
                    T_HX_out  = fp.temperature("water", P=P_LPmain_prev, h=h_HX_out)  # CONVERTED-NEEDS UNITS CHECK: eeslib uses SI; inputs confirmed in Pa/J/kg
                    rho_steam = fp.density("water", P=P_LPmain_prev, h=h_HX_out)      # CONVERTED-NEEDS UNITS CHECK: eeslib uses SI; inputs confirmed in Pa/J/kg
                    if rho_steam == 0.0:
                        rho_steam = 0.6  # [kg/m³] fallback
                    Vol_dot_SS_steam = m_dot_SS_steam / rho_steam
                else:  # no HTF flow
                    Vol_dot_HTF = 0.0
                    T_HTF_out   = HTF_T_in
                    T_HX_out    = T_SS_steam
                    Q_dot_HX    = 0.0
                    Eta_OD      = 0.0
            else:  # no steam flow
                Q_dot_HX    = 0.0
                Eta_OD      = 0.0
                rho_htf     = density(fnumd=self.Fluid_ID.v, T=HTF_T_in, P=HTF_P_in)
                Vol_dot_HTF = m_dot_HTF / rho_htf if rho_htf != 0.0 else 0.0
                T_HTF_out   = HTF_T_in
                h_HX_out    = h_SS_steam
                T_HX_out    = T_SS_steam

            # Write per-iteration reheater outputs
            self.m_dot_HTF_out.v = m_dot_HTF
            self.Vol_dot_HTF.v = Vol_dot_HTF
            self.HTF_P_out.v = HTF_P_in
            self.T_HTF_out.v = T_HTF_out
            self.m_dot_RH_steam.v = m_dot_SS_steam
            self.Vol_dot_RH_steam.v = Vol_dot_SS_steam
            self.T_HX_out.v = T_HX_out
            self.P_RH_steam.v = P_LPmain_prev
            self.Q_dot_HX.v = Q_dot_HX
            self.Eta_OD.v = Eta_OD

            # ------------------------------------------------------------------
            # HP MAIN PIPING — RK4 INTEGRATION
            # (replaces Fortran lines 2022–2200; GO TO 50 → goto50_HP_flag)
            # ------------------------------------------------------------------
            # Initialise piping solver
            pressure_tol_pipe = 10.0
            pressure_error = pressure_tol_pipe + 1.0
            LR = 0.2
            m_dot_min_pipe = 0.0
            m_dot_max_pipe = self.m_dot_HPT_d.v * 1.25
            whileiterations1 = 0.0
            vol_pipe_HP = math.pi / 4.0 * self.D_HPmain.v ** 2.0 * self.Length_HPmain.v

            m_dot_out_HP = m_dot_HP_warmup + m_dot_HP_drain + m_dot_HPT_in
            m_dot_in_HP = m_dot_out_HP   # initial guess: mass entering ≈ mass leaving
            h_in_HP = h_SGT_in
            m_dot_prev_HP = m_dot_in_HP
            pressure_error_prev = 0.0

            T_hpmain_prev = fp.temperature("water", P=P_HPmain_prev, h=h_HPmain_prev)  # CONVERTED-NEEDS UNITS CHECK: eeslib uses SI; inputs confirmed in Pa/J/kg

            # Initialise RK4 derivatives/heat-transfer rates (used in label-50 even if loop skipped)
            Q_dot_aa = Q_dot_bb = Q_dot_cc = Q_dot_dd = 0.0
            dPdt_aa = dPdt_bb = dPdt_cc = dPdt_dd = 0.0
            dhdt_aa = dhdt_bb = dhdt_cc = dhdt_dd = 0.0

            # RK4 state for HP pipe; updated in converged solution
            P_HPmain_new = P_HPmain_prev
            h_HPmain_new = h_HPmain_prev

            goto50_HP_flag = False  # flag to skip RK4 iteration (replaces GO TO 50)

            while abs(pressure_error) > pressure_tol_pipe:
                whileiterations1 += 1.0

                # Guard: no backflow allowed
                if m_dot_in_HP == 0.0 and pressure_error > 0.0:
                    goto50_HP_flag = True
                    break   # ← replaces GO TO 50

                # Guard: cap at max flow
                if m_dot_in_HP == m_dot_max_pipe and pressure_error < 0.0:
                    goto50_HP_flag = True
                    break   # ← replaces GO TO 50

                # Secant step for m_dot_in_HP (Fortran lines 2020-2087)
                if whileiterations1 > 2.0:
                    if m_dot_in_HP != m_dot_prev_HP:
                        _m_sl = (pressure_error - pressure_error_prev) / (m_dot_in_HP - m_dot_prev_HP)
                        _y_i  = pressure_error - m_dot_in_HP * _m_sl
                        if _m_sl != 0.0:
                            m_dot_prev_HP = m_dot_in_HP
                            pressure_error_prev = pressure_error
                            _m_new = -_y_i / _m_sl
                            m_dot_in_HP = m_dot_in_HP + (_m_new - m_dot_in_HP) * LR
                            m_dot_in_HP = max(min(m_dot_in_HP, m_dot_max_pipe), m_dot_min_pipe)
                        else:
                            pressure_error_prev = pressure_error
                            m_dot_prev_HP = m_dot_in_HP
                            if pressure_error > 0.0:
                                m_dot_in_HP = max(m_dot_in_HP - 5.0, m_dot_min_pipe)
                            else:
                                m_dot_in_HP = min(m_dot_in_HP - 5.0, m_dot_max_pipe)
                    else:
                        pressure_error_prev = pressure_error
                        m_dot_prev_HP = m_dot_in_HP
                        if pressure_error > 0.0:
                            m_dot_in_HP = max(m_dot_in_HP - 5.0, m_dot_min_pipe)
                        else:
                            m_dot_in_HP = min(m_dot_in_HP - 5.0, m_dot_max_pipe)
                else:
                    m_dot_prev_HP = m_dot_in_HP
                    pressure_error_prev = pressure_error
                    if pressure_error > 0.0:
                        m_dot_in_HP = max(m_dot_in_HP - 5.0, m_dot_min_pipe)
                    else:
                        m_dot_in_HP = min(m_dot_in_HP + 5.0, m_dot_max_pipe)

                # Outlet mass flow and initial outlet enthalpy (constant during RK4)
                m_dot_out_HP = m_dot_HP_warmup + m_dot_HP_drain + m_dot_HPT_in
                if m_dot_out_HP > 0.0:
                    h_out_HP = ((m_dot_HP_warmup * h_HPmain_prev + m_dot_HP_drain * h_HP_drain
                                 + m_dot_HPT_in * h_HPmain_prev) / m_dot_out_HP)
                else:
                    h_out_HP = 0.0
                m_dot_ave_HP = (m_dot_in_HP + m_dot_out_HP) / 2.0

                # ---- RK4 step aa (initial state: P_HPmain_prev, h_HPmain_prev) ----
                _rho = fp.density("water", P=P_HPmain_prev, h=h_HPmain_prev)      # CONVERTED-NEEDS UNITS CHECK: eeslib uses SI; inputs confirmed in Pa/J/kg
                _u   = fp.int_energy("water", P=P_HPmain_prev, h=h_HPmain_prev)   # CONVERTED-NEEDS UNITS CHECK: eeslib uses SI; inputs confirmed in Pa/J/kg
                _drhodh = drhodhcp(P_HPmain_prev, h_HPmain_prev, dh)
                _drhodp = drhodpch(P_HPmain_prev, h_HPmain_prev, dP)
                _dudh   = dudhcp(P_HPmain_prev, h_HPmain_prev, dh)
                _dudp   = dudpch(P_HPmain_prev, h_HPmain_prev, dP)
                _h_bar  = convection_dynamicpipe(
                    P_HPmain_prev, h_HPmain_prev, self.D_HPmain.v, m_dot_ave_HP, ff_HP, T_HP_pipe)
                _Q_max = self.mc_HPmain_pipe.v * (T_hpmain_prev - T_HP_pipe) / t_crit
                _Q_pipe = _h_bar * math.pi * self.D_HPmain.v * self.Length_HPmain.v * (T_hpmain_prev - T_HP_pipe)
                Q_dot_aa = min(_Q_pipe, _Q_max) if _Q_max > 0.0 else max(_Q_pipe, _Q_max)
                _denom = vol_pipe_HP * _rho * (_drhodh * _dudp - _drhodp * _dudh)
                dhdt_aa = (-(h_in_HP * _drhodp * m_dot_in_HP - h_out_HP * _drhodp * m_dot_out_HP
                              - _u * _drhodp * m_dot_in_HP + _u * _drhodp * m_dot_out_HP
                              - _dudp * m_dot_in_HP * _rho + _dudp * m_dot_out_HP * _rho
                              - Q_dot_aa * _drhodp) / _denom)
                dPdt_aa = ((m_dot_in_HP - m_dot_out_HP) / vol_pipe_HP - dhdt_aa * _drhodh) / _drhodp
                P_aa = P_HPmain_prev + dPdt_aa * t_crit / 2.0
                h_aa = h_HPmain_prev + dhdt_aa * t_crit / 2.0
                T_pipe_aa = (self.mc_HPmain_pipe.v * T_HP_pipe + Q_dot_aa * t_crit / 2.0) / self.mc_HPmain_pipe.v

                # ---- RK4 step bb (midpoint: P_aa, h_aa) ----
                _rho  = fp.density("water", P=P_aa, h=h_aa)        # CONVERTED-NEEDS UNITS CHECK: eeslib uses SI; inputs confirmed in Pa/J/kg
                _u    = fp.int_energy("water", P=P_aa, h=h_aa)     # CONVERTED-NEEDS UNITS CHECK: eeslib uses SI; inputs confirmed in Pa/J/kg
                T_aa  = fp.temperature("water", P=P_aa, h=h_aa)    # CONVERTED-NEEDS UNITS CHECK: eeslib uses SI; inputs confirmed in Pa/J/kg
                _h_sf = fp.enthalpy("water", P=P_aa, Q=0.0)        # CONVERTED-NEEDS UNITS CHECK: eeslib uses SI; inputs confirmed in Pa/J/kg
                _drhodh = drhodhcp(P_aa, h_aa, dh)
                _drhodp = drhodpch(P_aa, h_aa, dP)
                _dudh   = dudhcp(P_aa, h_aa, dh)
                _dudp   = dudpch(P_aa, h_aa, dP)
                if m_dot_out_HP > 0.0:
                    _h_out_step = ((m_dot_HP_warmup * h_aa + m_dot_HP_drain * _h_sf
                                    + m_dot_HPT_in * h_aa) / m_dot_out_HP)
                else:
                    _h_out_step = 0.0
                _h_bar  = convection_dynamicpipe(P_aa, h_aa, self.D_HPmain.v, m_dot_ave_HP, ff_HP, T_aa)
                _Q_max  = self.mc_HPmain_pipe.v * (T_aa - T_pipe_aa) / t_crit
                _Q_pipe = _h_bar * math.pi * self.D_HPmain.v * self.Length_HPmain.v * (T_aa - T_pipe_aa)
                Q_dot_bb = min(_Q_pipe, _Q_max) if _Q_max > 0.0 else max(_Q_pipe, _Q_max)
                _denom = vol_pipe_HP * _rho * (_drhodh * _dudp - _drhodp * _dudh)
                dhdt_bb = (-(h_in_HP * _drhodp * m_dot_in_HP - _h_out_step * _drhodp * m_dot_out_HP
                              - _u * _drhodp * m_dot_in_HP + _u * _drhodp * m_dot_out_HP
                              - _dudp * m_dot_in_HP * _rho + _dudp * m_dot_out_HP * _rho
                              - Q_dot_bb * _drhodp) / _denom)
                dPdt_bb = ((m_dot_in_HP - m_dot_out_HP) / vol_pipe_HP - dhdt_bb * _drhodh) / _drhodp
                P_bb = P_HPmain_prev + dPdt_bb * t_crit / 2.0
                h_bb = h_HPmain_prev + dhdt_bb * t_crit / 2.0
                T_pipe_bb = (self.mc_HPmain_pipe.v * T_HP_pipe + Q_dot_bb * t_crit / 2.0) / self.mc_HPmain_pipe.v

                # ---- RK4 step cc (midpoint: P_bb, h_bb) ----
                _rho  = fp.density("water", P=P_bb, h=h_bb)        # CONVERTED-NEEDS UNITS CHECK: eeslib uses SI; inputs confirmed in Pa/J/kg
                _u    = fp.int_energy("water", P=P_bb, h=h_bb)     # CONVERTED-NEEDS UNITS CHECK: eeslib uses SI; inputs confirmed in Pa/J/kg
                T_bb  = fp.temperature("water", P=P_bb, h=h_bb)    # CONVERTED-NEEDS UNITS CHECK: eeslib uses SI; inputs confirmed in Pa/J/kg
                _h_sf = fp.enthalpy("water", P=P_bb, Q=0.0)        # CONVERTED-NEEDS UNITS CHECK: eeslib uses SI; inputs confirmed in Pa/J/kg
                _drhodh = drhodhcp(P_bb, h_bb, dh)
                _drhodp = drhodpch(P_bb, h_bb, dP)
                _dudh   = dudhcp(P_bb, h_bb, dh)
                _dudp   = dudpch(P_bb, h_bb, dP)
                if m_dot_out_HP > 0.0:
                    _h_out_step = ((m_dot_HP_warmup * h_bb + m_dot_HP_drain * _h_sf
                                    + m_dot_HPT_in * h_bb) / m_dot_out_HP)
                else:
                    _h_out_step = 0.0
                _h_bar  = convection_dynamicpipe(P_bb, h_bb, self.D_HPmain.v, m_dot_ave_HP, ff_HP, T_bb)
                _Q_max  = self.mc_HPmain_pipe.v * (T_bb - T_pipe_bb) / t_crit
                _Q_pipe = _h_bar * math.pi * self.D_HPmain.v * self.Length_HPmain.v * (T_bb - T_pipe_bb)
                Q_dot_cc = min(_Q_pipe, _Q_max) if _Q_max > 0.0 else max(_Q_pipe, _Q_max)
                _denom = vol_pipe_HP * _rho * (_drhodh * _dudp - _drhodp * _dudh)
                dhdt_cc = (-(h_in_HP * _drhodp * m_dot_in_HP - _h_out_step * _drhodp * m_dot_out_HP
                              - _u * _drhodp * m_dot_in_HP + _u * _drhodp * m_dot_out_HP
                              - _dudp * m_dot_in_HP * _rho + _dudp * m_dot_out_HP * _rho
                              - Q_dot_cc * _drhodp) / _denom)
                dPdt_cc = ((m_dot_in_HP - m_dot_out_HP) / vol_pipe_HP - dhdt_cc * _drhodh) / _drhodp
                P_cc = P_HPmain_prev + dPdt_cc * t_crit
                h_cc = h_HPmain_prev + dhdt_cc * t_crit
                T_pipe_cc = (self.mc_HPmain_pipe.v * T_HP_pipe + Q_dot_cc * t_crit) / self.mc_HPmain_pipe.v

                # ---- RK4 step dd (full step: P_cc, h_cc) ----
                _rho  = fp.density("water", P=P_cc, h=h_cc)        # CONVERTED-NEEDS UNITS CHECK: eeslib uses SI; inputs confirmed in Pa/J/kg
                _u    = fp.int_energy("water", P=P_cc, h=h_cc)     # CONVERTED-NEEDS UNITS CHECK: eeslib uses SI; inputs confirmed in Pa/J/kg
                T_cc  = fp.temperature("water", P=P_cc, h=h_cc)    # CONVERTED-NEEDS UNITS CHECK: eeslib uses SI; inputs confirmed in Pa/J/kg
                _h_sf = fp.enthalpy("water", P=P_cc, Q=0.0)        # CONVERTED-NEEDS UNITS CHECK: eeslib uses SI; inputs confirmed in Pa/J/kg
                _drhodh = drhodhcp(P_cc, h_cc, dh)
                _drhodp = drhodpch(P_cc, h_cc, dP)
                _dudh   = dudhcp(P_cc, h_cc, dh)
                _dudp   = dudpch(P_cc, h_cc, dP)
                if m_dot_out_HP > 0.0:
                    _h_out_step = ((m_dot_HP_warmup * h_cc + m_dot_HP_drain * _h_sf
                                    + m_dot_HPT_in * h_cc) / m_dot_out_HP)
                else:
                    _h_out_step = 0.0
                _h_bar  = convection_dynamicpipe(P_cc, h_cc, self.D_HPmain.v, m_dot_ave_HP, ff_HP, T_cc)
                _Q_max  = self.mc_HPmain_pipe.v * (T_cc - T_pipe_cc) / t_crit
                _Q_pipe = _h_bar * math.pi * self.D_HPmain.v * self.Length_HPmain.v * (T_cc - T_pipe_cc)
                Q_dot_dd = min(_Q_pipe, _Q_max) if _Q_max > 0.0 else max(_Q_pipe, _Q_max)
                _denom = vol_pipe_HP * _rho * (_drhodh * _dudp - _drhodp * _dudh)
                dhdt_dd = (-(h_in_HP * _drhodp * m_dot_in_HP - _h_out_step * _drhodp * m_dot_out_HP
                              - _u * _drhodp * m_dot_in_HP + _u * _drhodp * m_dot_out_HP
                              - _dudp * m_dot_in_HP * _rho + _dudp * m_dot_out_HP * _rho
                              - Q_dot_cc * _drhodp) / _denom)  # NOTE: Fortran uses Q_dot_cc (not dd) — replicating bug
                dPdt_dd = ((m_dot_in_HP - m_dot_out_HP) / vol_pipe_HP - dhdt_dd * _drhodh) / _drhodp

                # RK4 weighted sum → converge on new HP-main state
                P_HPmain_new = (P_HPmain_prev
                                + (dPdt_aa + 2.0 * dPdt_bb + 2.0 * dPdt_cc + dPdt_dd) * t_crit / 6.0)
                h_HPmain_new = (h_HPmain_prev
                                + (dhdt_aa + 2.0 * dhdt_bb + 2.0 * dhdt_cc + dhdt_dd) * t_crit / 6.0)
                pressure_error = P_HPmain_new - P_SGT_in

            # label 50 equivalent: update HP piping state regardless of how we got here
            P_HPmain_prev = P_HPmain_new
            h_HPmain_prev = h_HPmain_new
            m_HPmain_prev = m_HPmain_prev + (m_dot_in_HP - m_dot_out_HP) * t_crit
            Q_dot_hpmain_avg = (Q_dot_aa + Q_dot_bb + Q_dot_cc + Q_dot_dd) / 4.0
            T_HP_pipe = (self.mc_HPmain_pipe.v * T_HP_pipe + Q_dot_hpmain_avg * t_crit) / self.mc_HPmain_pipe.v
            T_HPmain = fp.temperature("water", P=P_HPmain_prev, h=h_HPmain_prev)  # CONVERTED-NEEDS UNITS CHECK: eeslib uses SI; inputs confirmed in Pa/J/kg
            x_HPmain = fp.quality("water", P=P_HPmain_prev, h=h_HPmain_prev)      # CONVERTED-NEEDS UNITS CHECK: eeslib uses SI; inputs confirmed in Pa/J/kg
            _h_sat_g_HP = fp.enthalpy("water", P=P_HPmain_prev, Q=1.0)            # CONVERTED-NEEDS UNITS CHECK: eeslib uses SI; inputs confirmed in Pa/J/kg
            if h_HPmain_prev > _h_sat_g_HP:
                x_HPmain = 100.0

            m_dot_SGT_tot += m_dot_HP_bypass * t_crit + m_dot_HP_AUX * t_crit + m_dot_in_HP * t_crit

            # Write per-iteration HP-pipe outputs
            self.P_HPmain.v = P_HPmain_prev
            self.T_HPmain.v = T_HPmain
            self.x_HPmain.v = x_HPmain
            self.T_HP_pipe.v = T_HP_pipe
            self.ff_HPmain.v = ff_HP

            # ------------------------------------------------------------------
            # LP MAIN PIPING — RK4 INTEGRATION
            # (replaces Fortran lines 2218–2390; no explicit GOTO but same pattern)
            # ------------------------------------------------------------------
            m_dot_in_LP = m_dot_SS_steam
            h_in_LP = h_HX_out
            m_dot_out_LP = (m_dot_LP_bypass + m_dot_LP_warmup + m_dot_LP_drain
                            + m_dot_LP_AUX + m_dot_LPT_in)
            vol_pipe_LP = self.D_LPmain.v ** 2.0 * math.pi / 4.0 * self.Length_LPmain.v

            # Initial state properties for the aa step (Fortran: FIT_PH + FIT_PQ before Q_dot_aa)
            _rho_LP = fp.density("water", P=P_LPmain_prev, h=h_LPmain_prev)      # CONVERTED-NEEDS UNITS CHECK: eeslib uses SI; inputs confirmed in Pa/J/kg
            _u_LP   = fp.int_energy("water", P=P_LPmain_prev, h=h_LPmain_prev)   # CONVERTED-NEEDS UNITS CHECK: eeslib uses SI; inputs confirmed in Pa/J/kg
            T_lpmain = fp.temperature("water", P=P_LPmain_prev, h=h_LPmain_prev) # CONVERTED-NEEDS UNITS CHECK: eeslib uses SI; inputs confirmed in Pa/J/kg
            _h_sf0  = fp.enthalpy("water", P=P_LPmain_prev, Q=0.0)              # CONVERTED-NEEDS UNITS CHECK: eeslib uses SI; inputs confirmed in Pa/J/kg

            # Initial outlet enthalpy weighted by outflow type (drain uses sat-liquid)
            if m_dot_out_LP > 0.0:
                h_out_LP = ((h_LP_bypass * m_dot_LP_bypass
                             + m_dot_LP_warmup * h_LPmain_prev
                             + m_dot_LP_AUX * h_LPmain_prev
                             + m_dot_LP_drain * h_LP_drain
                             + m_dot_LPT_in * h_LPmain_prev) / m_dot_out_LP)
            else:
                h_out_LP = 0.0
            m_dot_ave_LP = (m_dot_in_LP + m_dot_out_LP) / 2.0

            # NOTE: Fortran reads ff_guess=getOutputValue(93) but then uses HP ff in call — replicate with ff_LP
            # ---- RK4 step aa (initial state: P_LPmain_prev, h_LPmain_prev) ----
            _h_bar  = convection_dynamicpipe(
                P_LPmain_prev, h_LPmain_prev, self.D_LPmain.v, m_dot_ave_LP, ff_LP, T_LP_pipe)
            _Q_max  = self.mc_LPmain_pipe.v * (T_lpmain - T_LP_pipe) / t_crit
            _Q_pipe = _h_bar * math.pi * self.D_LPmain.v * self.Length_LPmain.v * (T_lpmain - T_LP_pipe)
            Q_dot_aa_LP = min(_Q_pipe, _Q_max) if _Q_max > 0.0 else max(_Q_pipe, _Q_max)
            _drhodh = drhodhcp(P_LPmain_prev, h_LPmain_prev, dh)
            _drhodp = drhodpch(P_LPmain_prev, h_LPmain_prev, dP)
            _dudh   = dudhcp(P_LPmain_prev, h_LPmain_prev, dh)
            _dudp   = dudpch(P_LPmain_prev, h_LPmain_prev, dP)
            _denom_LP = _rho_LP * vol_pipe_LP * (_drhodh * _dudp - _drhodp * _dudh)
            dPdt_aa_LP = ((h_in_LP * _drhodh * m_dot_in_LP - h_out_LP * _drhodh * m_dot_out_LP
                           - _drhodh * m_dot_in_LP * _u_LP + _drhodh * m_dot_out_LP * _u_LP
                           - _dudh * m_dot_in_LP * _rho_LP + _dudh * m_dot_out_LP * _rho_LP
                           - Q_dot_aa_LP * _drhodh) / _denom_LP)
            dhdt_aa_LP = ((m_dot_in_LP - m_dot_out_LP) / vol_pipe_LP - _drhodp * dPdt_aa_LP) / _drhodh
            P_aa_LP = P_LPmain_prev + dPdt_aa_LP * t_crit / 2.0
            h_aa_LP = h_LPmain_prev + dhdt_aa_LP * t_crit / 2.0
            T_pipe_aa_LP = (self.mc_LPmain_pipe.v * T_LP_pipe + Q_dot_aa_LP * t_crit / 2.0) / self.mc_LPmain_pipe.v

            # ---- RK4 step bb (midpoint: P_aa_LP, h_aa_LP) ----
            _rho_LP = fp.density("water", P=P_aa_LP, h=h_aa_LP)        # CONVERTED-NEEDS UNITS CHECK: eeslib uses SI; inputs confirmed in Pa/J/kg
            _u_LP   = fp.int_energy("water", P=P_aa_LP, h=h_aa_LP)     # CONVERTED-NEEDS UNITS CHECK: eeslib uses SI; inputs confirmed in Pa/J/kg
            T_aa_LP = fp.temperature("water", P=P_aa_LP, h=h_aa_LP)    # CONVERTED-NEEDS UNITS CHECK: eeslib uses SI; inputs confirmed in Pa/J/kg
            _h_sf   = fp.enthalpy("water", P=P_aa_LP, Q=0.0)           # CONVERTED-NEEDS UNITS CHECK: eeslib uses SI; inputs confirmed in Pa/J/kg
            if m_dot_out_LP > 0.0:
                h_out_LP = ((h_aa_LP * m_dot_LP_bypass + m_dot_LP_warmup * h_aa_LP
                             + m_dot_LP_AUX * h_aa_LP + m_dot_LP_drain * _h_sf
                             + m_dot_LPT_in * h_aa_LP) / m_dot_out_LP)
            _drhodh = drhodhcp(P_aa_LP, h_aa_LP, dh)
            _drhodp = drhodpch(P_aa_LP, h_aa_LP, dP)
            _dudh   = dudhcp(P_aa_LP, h_aa_LP, dh)
            _dudp   = dudpch(P_aa_LP, h_aa_LP, dP)
            _h_bar  = convection_dynamicpipe(P_aa_LP, h_aa_LP, self.D_LPmain.v, m_dot_ave_LP, ff_LP, T_pipe_aa_LP)
            _Q_max  = self.mc_LPmain_pipe.v * (T_aa_LP - T_pipe_aa_LP) / t_crit
            _Q_pipe = _h_bar * math.pi * self.D_LPmain.v * self.Length_LPmain.v * (T_aa_LP - T_pipe_aa_LP)
            Q_dot_bb_LP = min(_Q_pipe, _Q_max) if _Q_max > 0.0 else max(_Q_pipe, _Q_max)
            _denom_LP = _rho_LP * vol_pipe_LP * (_drhodh * _dudp - _drhodp * _dudh)
            dPdt_bb_LP = ((h_in_LP * _drhodh * m_dot_in_LP - h_out_LP * _drhodh * m_dot_out_LP
                           - _drhodh * m_dot_in_LP * _u_LP + _drhodh * m_dot_out_LP * _u_LP
                           - _dudh * m_dot_in_LP * _rho_LP + _dudh * m_dot_out_LP * _rho_LP
                           - Q_dot_bb_LP * _drhodh) / _denom_LP)
            dhdt_bb_LP = ((m_dot_in_LP - m_dot_out_LP) / vol_pipe_LP - _drhodp * dPdt_bb_LP) / _drhodh
            P_bb_LP = P_LPmain_prev + dPdt_bb_LP * t_crit / 2.0
            h_bb_LP = h_LPmain_prev + dhdt_bb_LP * t_crit / 2.0
            T_pipe_bb_LP = (self.mc_LPmain_pipe.v * T_LP_pipe + Q_dot_bb_LP * t_crit / 2.0) / self.mc_LPmain_pipe.v

            # ---- RK4 step cc (midpoint: P_bb_LP, h_bb_LP) ----
            _rho_LP = fp.density("water", P=P_bb_LP, h=h_bb_LP)        # CONVERTED-NEEDS UNITS CHECK: eeslib uses SI; inputs confirmed in Pa/J/kg
            _u_LP   = fp.int_energy("water", P=P_bb_LP, h=h_bb_LP)     # CONVERTED-NEEDS UNITS CHECK: eeslib uses SI; inputs confirmed in Pa/J/kg
            T_bb_LP = fp.temperature("water", P=P_bb_LP, h=h_bb_LP)    # CONVERTED-NEEDS UNITS CHECK: eeslib uses SI; inputs confirmed in Pa/J/kg
            _h_sf   = fp.enthalpy("water", P=P_bb_LP, Q=0.0)           # CONVERTED-NEEDS UNITS CHECK: eeslib uses SI; inputs confirmed in Pa/J/kg
            if m_dot_out_LP > 0.0:
                h_out_LP = ((h_bb_LP * m_dot_LP_bypass + m_dot_LP_warmup * h_bb_LP
                             + m_dot_LP_AUX * h_bb_LP + m_dot_LP_drain * _h_sf
                             + m_dot_LPT_in * h_bb_LP) / m_dot_out_LP)
            _drhodh = drhodhcp(P_bb_LP, h_bb_LP, dh)
            _drhodp = drhodpch(P_bb_LP, h_bb_LP, dP)
            _dudh   = dudhcp(P_bb_LP, h_bb_LP, dh)
            _dudp   = dudpch(P_bb_LP, h_bb_LP, dP)
            _h_bar  = convection_dynamicpipe(P_bb_LP, h_bb_LP, self.D_LPmain.v, m_dot_ave_LP, ff_LP, T_pipe_bb_LP)
            _Q_max  = self.mc_LPmain_pipe.v * (T_bb_LP - T_pipe_bb_LP) / t_crit
            _Q_pipe = _h_bar * math.pi * self.D_LPmain.v * self.Length_LPmain.v * (T_bb_LP - T_pipe_bb_LP)
            Q_dot_cc_LP = min(_Q_pipe, _Q_max) if _Q_max > 0.0 else max(_Q_pipe, _Q_max)
            _denom_LP = _rho_LP * vol_pipe_LP * (_drhodh * _dudp - _drhodp * _dudh)
            dPdt_cc_LP = ((h_in_LP * _drhodh * m_dot_in_LP - h_out_LP * _drhodh * m_dot_out_LP
                           - _drhodh * m_dot_in_LP * _u_LP + _drhodh * m_dot_out_LP * _u_LP
                           - _dudh * m_dot_in_LP * _rho_LP + _dudh * m_dot_out_LP * _rho_LP
                           - Q_dot_cc_LP * _drhodh) / _denom_LP)
            dhdt_cc_LP = ((m_dot_in_LP - m_dot_out_LP) / vol_pipe_LP - _drhodp * dPdt_cc_LP) / _drhodh
            P_cc_LP = P_LPmain_prev + dPdt_cc_LP * t_crit
            h_cc_LP = h_LPmain_prev + dhdt_cc_LP * t_crit
            T_pipe_cc_LP = (self.mc_LPmain_pipe.v * T_LP_pipe + Q_dot_cc_LP * t_crit) / self.mc_LPmain_pipe.v

            # ---- RK4 step dd (full step: P_cc_LP, h_cc_LP) ----
            _rho_LP = fp.density("water", P=P_cc_LP, h=h_cc_LP)        # CONVERTED-NEEDS UNITS CHECK: eeslib uses SI; inputs confirmed in Pa/J/kg
            _u_LP   = fp.int_energy("water", P=P_cc_LP, h=h_cc_LP)     # CONVERTED-NEEDS UNITS CHECK: eeslib uses SI; inputs confirmed in Pa/J/kg
            T_cc_LP = fp.temperature("water", P=P_cc_LP, h=h_cc_LP)    # CONVERTED-NEEDS UNITS CHECK: eeslib uses SI; inputs confirmed in Pa/J/kg
            _h_sf   = fp.enthalpy("water", P=P_cc_LP, Q=0.0)           # CONVERTED-NEEDS UNITS CHECK: eeslib uses SI; inputs confirmed in Pa/J/kg
            if m_dot_out_LP > 0.0:
                h_out_LP = ((h_cc_LP * m_dot_LP_bypass + m_dot_LP_warmup * h_cc_LP
                             + m_dot_LP_AUX * h_cc_LP + m_dot_LP_drain * _h_sf
                             + m_dot_LPT_in * h_cc_LP) / m_dot_out_LP)
            _drhodh = drhodhcp(P_cc_LP, h_cc_LP, dh)
            _drhodp = drhodpch(P_cc_LP, h_cc_LP, dP)
            _dudh   = dudhcp(P_cc_LP, h_cc_LP, dh)
            _dudp   = dudpch(P_cc_LP, h_cc_LP, dP)
            _h_bar  = convection_dynamicpipe(P_cc_LP, h_cc_LP, self.D_LPmain.v, m_dot_ave_LP, ff_LP, T_pipe_cc_LP)
            _Q_max  = self.mc_LPmain_pipe.v * (T_cc_LP - T_pipe_cc_LP) / t_crit
            _Q_pipe = _h_bar * math.pi * self.D_LPmain.v * self.Length_LPmain.v * (T_cc_LP - T_pipe_cc_LP)
            Q_dot_dd_LP = min(_Q_pipe, _Q_max) if _Q_max > 0.0 else max(_Q_pipe, _Q_max)
            _denom_LP = _rho_LP * vol_pipe_LP * (_drhodh * _dudp - _drhodp * _dudh)
            dPdt_dd_LP = ((h_in_LP * _drhodh * m_dot_in_LP - h_out_LP * _drhodh * m_dot_out_LP
                           - _drhodh * m_dot_in_LP * _u_LP + _drhodh * m_dot_out_LP * _u_LP
                           - _dudh * m_dot_in_LP * _rho_LP + _dudh * m_dot_out_LP * _rho_LP
                           - Q_dot_dd_LP * _drhodh) / _denom_LP)
            dhdt_dd_LP = ((m_dot_in_LP - m_dot_out_LP) / vol_pipe_LP - _drhodp * dPdt_dd_LP) / _drhodh

            # RK4 weighted sum → update LP main state
            P_LPmain_prev = (P_LPmain_prev
                             + (dPdt_aa_LP + 2.0 * dPdt_bb_LP + 2.0 * dPdt_cc_LP + dPdt_dd_LP) * t_crit / 6.0)
            h_LPmain_prev = (h_LPmain_prev
                             + (dhdt_aa_LP + 2.0 * dhdt_bb_LP + 2.0 * dhdt_cc_LP + dhdt_dd_LP) * t_crit / 6.0)
            Q_dot_lpmain_avg = (Q_dot_aa_LP + Q_dot_bb_LP + Q_dot_cc_LP + Q_dot_dd_LP) / 4.0
            m_LPmain_prev = m_LPmain_prev + (m_dot_in_LP - m_dot_out_LP) * t_crit
            T_LP_pipe = (self.mc_LPmain_pipe.v * T_LP_pipe + Q_dot_lpmain_avg * t_crit) / self.mc_LPmain_pipe.v
            T_LPmain   = fp.temperature("water", P=P_LPmain_prev, h=h_LPmain_prev)  # CONVERTED-NEEDS UNITS CHECK: eeslib uses SI; inputs confirmed in Pa/J/kg
            x_LPmain   = fp.quality("water", P=P_LPmain_prev, h=h_LPmain_prev)      # CONVERTED-NEEDS UNITS CHECK: eeslib uses SI; inputs confirmed in Pa/J/kg
            _h_sat_g_LP = fp.enthalpy("water", P=P_LPmain_prev, Q=1.0)              # CONVERTED-NEEDS UNITS CHECK: eeslib uses SI; inputs confirmed in Pa/J/kg
            if h_LPmain_prev > _h_sat_g_LP:
                x_LPmain = 100.0

            self.P_LPmain.v = P_LPmain_prev
            self.T_LPmain.v = T_LPmain
            self.x_LPmain.v = x_LPmain
            self.T_LP_pipe.v = T_LP_pipe
            self.ff_LPmain.v = ff_LP

            # ------------------------------------------------------------------
            # AUX PIPING — RK4 INTEGRATION
            # (replaces Fortran lines 2392–2539; no explicit GOTO)
            # ------------------------------------------------------------------
            m_dot_in_AUX = m_dot_LP_AUX + m_dot_HP_AUX
            h_in_AUX = (
                (m_dot_LP_AUX * h_LP_AUX + m_dot_HP_AUX * h_HP_AUX) / m_dot_in_AUX
                if m_dot_in_AUX > 0.0
                else 0.0
            )
            m_dot_out_AUX = m_dot_ts + m_dot_AUX_DA
            vol_pipe_AUX = math.pi / 4.0 * self.D_AuxLine.v ** 2.0 * self.Length_AuxLine.v

            # Initial state properties for AUX aa step
            _rho_AUX = fp.density("water", P=P_Aux_prev, h=h_Aux_prev)      # CONVERTED-NEEDS UNITS CHECK: eeslib uses SI; inputs confirmed in Pa/J/kg
            _u_AUX   = fp.int_energy("water", P=P_Aux_prev, h=h_Aux_prev)   # CONVERTED-NEEDS UNITS CHECK: eeslib uses SI; inputs confirmed in Pa/J/kg
            T_aux    = fp.temperature("water", P=P_Aux_prev, h=h_Aux_prev)  # CONVERTED-NEEDS UNITS CHECK: eeslib uses SI; inputs confirmed in Pa/J/kg

            # AUX h_out = current state enthalpy (no split drain correction like LP)
            h_out_AUX = h_Aux_prev
            m_dot_ave_AUX = (m_dot_in_AUX + m_dot_out_AUX) / 2.0

            # ---- RK4 step aa ----
            _h_bar   = convection_dynamicpipe(
                P_Aux_prev, h_Aux_prev, self.D_AuxLine.v, m_dot_ave_AUX, ff_AUX, T_Aux_pipe)
            _Q_max   = self.mc_AUX_pipe.v * (T_aux - T_Aux_pipe) / t_crit
            _Q_pipe  = _h_bar * math.pi * self.D_AuxLine.v * self.Length_AuxLine.v * (T_aux - T_Aux_pipe)
            Q_dot_aa_AUX = min(_Q_pipe, _Q_max) if _Q_max > 0.0 else max(_Q_pipe, _Q_max)
            _drhodh  = drhodhcp(P_Aux_prev, h_Aux_prev, dh)
            _drhodp  = drhodpch(P_Aux_prev, h_Aux_prev, dP)
            _dudh    = dudhcp(P_Aux_prev, h_Aux_prev, dh)
            _dudp    = dudpch(P_Aux_prev, h_Aux_prev, dP)
            _denom_AUX = vol_pipe_AUX * _rho_AUX * (_drhodh * _dudp - _drhodp * _dudh)
            dhdt_aa_AUX = (-(h_in_AUX * _drhodp * m_dot_in_AUX - h_out_AUX * _drhodp * m_dot_out_AUX
                              - _u_AUX * _drhodp * m_dot_in_AUX + _u_AUX * _drhodp * m_dot_out_AUX
                              - _dudp * m_dot_in_AUX * _rho_AUX + _dudp * m_dot_out_AUX * _rho_AUX
                              - Q_dot_aa_AUX * _drhodp) / _denom_AUX)
            dPdt_aa_AUX = ((m_dot_in_AUX - m_dot_out_AUX) / vol_pipe_AUX - dhdt_aa_AUX * _drhodh) / _drhodp
            P_aa_AUX = P_Aux_prev + dPdt_aa_AUX * t_crit / 2.0
            h_aa_AUX = h_Aux_prev + dhdt_aa_AUX * t_crit / 2.0
            T_pipe_aa_AUX = (self.mc_AUX_pipe.v * T_Aux_pipe + Q_dot_aa_AUX * t_crit / 2.0) / self.mc_AUX_pipe.v

            # ---- RK4 step bb ----
            _rho_AUX = fp.density("water", P=P_aa_AUX, h=h_aa_AUX)        # CONVERTED-NEEDS UNITS CHECK: eeslib uses SI; inputs confirmed in Pa/J/kg
            _u_AUX   = fp.int_energy("water", P=P_aa_AUX, h=h_aa_AUX)     # CONVERTED-NEEDS UNITS CHECK: eeslib uses SI; inputs confirmed in Pa/J/kg
            T_aa_AUX = fp.temperature("water", P=P_aa_AUX, h=h_aa_AUX)    # CONVERTED-NEEDS UNITS CHECK: eeslib uses SI; inputs confirmed in Pa/J/kg
            h_out_AUX = h_aa_AUX  # AUX h_out updates to current state enthalpy
            _drhodh  = drhodhcp(P_aa_AUX, h_aa_AUX, dh)
            _drhodp  = drhodpch(P_aa_AUX, h_aa_AUX, dP)
            _dudh    = dudhcp(P_aa_AUX, h_aa_AUX, dh)
            _dudp    = dudpch(P_aa_AUX, h_aa_AUX, dP)
            _h_bar   = convection_dynamicpipe(P_aa_AUX, h_aa_AUX, self.D_AuxLine.v, m_dot_ave_AUX, ff_AUX, T_pipe_aa_AUX)
            _Q_max   = self.mc_AUX_pipe.v * (T_aa_AUX - T_pipe_aa_AUX) / t_crit
            _Q_pipe  = _h_bar * math.pi * self.D_AuxLine.v * self.Length_AuxLine.v * (T_aa_AUX - T_pipe_aa_AUX)
            Q_dot_bb_AUX = min(_Q_pipe, _Q_max) if _Q_max > 0.0 else max(_Q_pipe, _Q_max)
            _denom_AUX = vol_pipe_AUX * _rho_AUX * (_drhodh * _dudp - _drhodp * _dudh)
            dhdt_bb_AUX = (-(h_in_AUX * _drhodp * m_dot_in_AUX - h_out_AUX * _drhodp * m_dot_out_AUX
                              - _u_AUX * _drhodp * m_dot_in_AUX + _u_AUX * _drhodp * m_dot_out_AUX
                              - _dudp * m_dot_in_AUX * _rho_AUX + _dudp * m_dot_out_AUX * _rho_AUX
                              - Q_dot_bb_AUX * _drhodp) / _denom_AUX)
            dPdt_bb_AUX = ((m_dot_in_AUX - m_dot_out_AUX) / vol_pipe_AUX - dhdt_bb_AUX * _drhodh) / _drhodp
            P_bb_AUX = P_Aux_prev + dPdt_bb_AUX * t_crit / 2.0
            h_bb_AUX = h_Aux_prev + dhdt_bb_AUX * t_crit / 2.0
            T_pipe_bb_AUX = (self.mc_AUX_pipe.v * T_Aux_pipe + Q_dot_bb_AUX * t_crit / 2.0) / self.mc_AUX_pipe.v

            # ---- RK4 step cc ----
            _rho_AUX = fp.density("water", P=P_bb_AUX, h=h_bb_AUX)        # CONVERTED-NEEDS UNITS CHECK: eeslib uses SI; inputs confirmed in Pa/J/kg
            _u_AUX   = fp.int_energy("water", P=P_bb_AUX, h=h_bb_AUX)     # CONVERTED-NEEDS UNITS CHECK: eeslib uses SI; inputs confirmed in Pa/J/kg
            T_bb_AUX = fp.temperature("water", P=P_bb_AUX, h=h_bb_AUX)    # CONVERTED-NEEDS UNITS CHECK: eeslib uses SI; inputs confirmed in Pa/J/kg
            h_out_AUX = h_bb_AUX
            _drhodh  = drhodhcp(P_bb_AUX, h_bb_AUX, dh)
            _drhodp  = drhodpch(P_bb_AUX, h_bb_AUX, dP)
            _dudh    = dudhcp(P_bb_AUX, h_bb_AUX, dh)
            _dudp    = dudpch(P_bb_AUX, h_bb_AUX, dP)
            _h_bar   = convection_dynamicpipe(P_bb_AUX, h_bb_AUX, self.D_AuxLine.v, m_dot_ave_AUX, ff_AUX, T_pipe_bb_AUX)
            _Q_max   = self.mc_AUX_pipe.v * (T_bb_AUX - T_pipe_bb_AUX) / t_crit
            _Q_pipe  = _h_bar * math.pi * self.D_AuxLine.v * self.Length_AuxLine.v * (T_bb_AUX - T_pipe_bb_AUX)
            Q_dot_cc_AUX = min(_Q_pipe, _Q_max) if _Q_max > 0.0 else max(_Q_pipe, _Q_max)
            _denom_AUX = vol_pipe_AUX * _rho_AUX * (_drhodh * _dudp - _drhodp * _dudh)
            dhdt_cc_AUX = (-(h_in_AUX * _drhodp * m_dot_in_AUX - h_out_AUX * _drhodp * m_dot_out_AUX
                              - _u_AUX * _drhodp * m_dot_in_AUX + _u_AUX * _drhodp * m_dot_out_AUX
                              - _dudp * m_dot_in_AUX * _rho_AUX + _dudp * m_dot_out_AUX * _rho_AUX
                              - Q_dot_cc_AUX * _drhodp) / _denom_AUX)
            dPdt_cc_AUX = ((m_dot_in_AUX - m_dot_out_AUX) / vol_pipe_AUX - dhdt_cc_AUX * _drhodh) / _drhodp
            P_cc_AUX = P_Aux_prev + dPdt_cc_AUX * t_crit
            h_cc_AUX = h_Aux_prev + dhdt_cc_AUX * t_crit
            T_pipe_cc_AUX = (self.mc_AUX_pipe.v * T_Aux_pipe + Q_dot_cc_AUX * t_crit) / self.mc_AUX_pipe.v

            # ---- RK4 step dd ----
            _rho_AUX = fp.density("water", P=P_cc_AUX, h=h_cc_AUX)        # CONVERTED-NEEDS UNITS CHECK: eeslib uses SI; inputs confirmed in Pa/J/kg
            _u_AUX   = fp.int_energy("water", P=P_cc_AUX, h=h_cc_AUX)     # CONVERTED-NEEDS UNITS CHECK: eeslib uses SI; inputs confirmed in Pa/J/kg
            T_cc_AUX = fp.temperature("water", P=P_cc_AUX, h=h_cc_AUX)    # CONVERTED-NEEDS UNITS CHECK: eeslib uses SI; inputs confirmed in Pa/J/kg
            h_out_AUX = h_cc_AUX
            _drhodh  = drhodhcp(P_cc_AUX, h_cc_AUX, dh)
            _drhodp  = drhodpch(P_cc_AUX, h_cc_AUX, dP)
            _dudh    = dudhcp(P_cc_AUX, h_cc_AUX, dh)
            _dudp    = dudpch(P_cc_AUX, h_cc_AUX, dP)
            _h_bar   = convection_dynamicpipe(P_cc_AUX, h_cc_AUX, self.D_AuxLine.v, m_dot_ave_AUX, ff_AUX, T_pipe_cc_AUX)
            _Q_max   = self.mc_AUX_pipe.v * (T_cc_AUX - T_pipe_cc_AUX) / t_crit
            _Q_pipe  = _h_bar * math.pi * self.D_AuxLine.v * self.Length_AuxLine.v * (T_cc_AUX - T_pipe_cc_AUX)
            Q_dot_dd_AUX = min(_Q_pipe, _Q_max) if _Q_max > 0.0 else max(_Q_pipe, _Q_max)
            _denom_AUX = vol_pipe_AUX * _rho_AUX * (_drhodh * _dudp - _drhodp * _dudh)
            dhdt_dd_AUX = (-(h_in_AUX * _drhodp * m_dot_in_AUX - h_out_AUX * _drhodp * m_dot_out_AUX
                              - _u_AUX * _drhodp * m_dot_in_AUX + _u_AUX * _drhodp * m_dot_out_AUX
                              - _dudp * m_dot_in_AUX * _rho_AUX + _dudp * m_dot_out_AUX * _rho_AUX
                              - Q_dot_dd_AUX * _drhodp) / _denom_AUX)
            dPdt_dd_AUX = ((m_dot_in_AUX - m_dot_out_AUX) / vol_pipe_AUX - dhdt_dd_AUX * _drhodh) / _drhodp

            # RK4 weighted sum → update AUX state
            P_Aux_prev = (P_Aux_prev
                          + (dPdt_aa_AUX + 2.0 * dPdt_bb_AUX + 2.0 * dPdt_cc_AUX + dPdt_dd_AUX) * t_crit / 6.0)
            h_Aux_prev = (h_Aux_prev
                          + (dhdt_aa_AUX + 2.0 * dhdt_bb_AUX + 2.0 * dhdt_cc_AUX + dhdt_dd_AUX) * t_crit / 6.0)
            Q_dot_aux_avg = (Q_dot_aa_AUX + Q_dot_bb_AUX + Q_dot_cc_AUX + Q_dot_dd_AUX) / 4.0
            m_Aux_prev = m_Aux_prev + (m_dot_in_AUX - m_dot_out_AUX) * t_crit
            T_Aux_pipe = (self.mc_AUX_pipe.v * T_Aux_pipe + Q_dot_aux_avg * t_crit) / self.mc_AUX_pipe.v
            T_AUX  = fp.temperature("water", P=P_Aux_prev, h=h_Aux_prev)   # CONVERTED-NEEDS UNITS CHECK: eeslib uses SI; inputs confirmed in Pa/J/kg
            x_AUX  = fp.quality("water", P=P_Aux_prev, h=h_Aux_prev)       # CONVERTED-NEEDS UNITS CHECK: eeslib uses SI; inputs confirmed in Pa/J/kg
            _h_sat_g_AUX = fp.enthalpy("water", P=P_Aux_prev, Q=1.0)       # CONVERTED-NEEDS UNITS CHECK: eeslib uses SI; inputs confirmed in Pa/J/kg
            if h_Aux_prev > _h_sat_g_AUX:
                x_LPmain = 100.0  # NOTE: Fortran bug — sets x_LPmain instead of x_AUX; replicating

            self.P_AUX.v = P_Aux_prev
            self.T_AUX.v = T_AUX
            self.x_AUX.v = x_AUX
            self.T_AUX_pipe.v = T_Aux_pipe
            self.ff_AUX.v = ff_AUX

        # =========================================================================
        # END OF SUB-TIMESTEP LOOP — normalise accumulators and set final outputs
        # (replaces Fortran lines 2541–2639)
        # =========================================================================
        m_dot_HPT_exh_tot /= ts
        W_dot_HPT_tot /= ts
        m_dot_LPT_exh_tot /= ts
        W_dot_LPT_tot /= ts

        # Enthalpy-weighted average for LPT exhaust
        if m_dot_LPT_exh_tot > 0.0:
            h_LPT_exh_tot = h_LPT_exh_tot / m_dot_LPT_exh_tot
        else:
            h_LPT_exh_tot = 0.0

        W_dot_tot = W_dot_HPT_tot + W_dot_LPT_tot

        m_dot_HP_warmup_tot /= ts
        if m_dot_HP_warmup_tot > 0.0:
            h_HP_warmup_tot = h_HP_warmup_tot / m_dot_HP_warmup_tot
        else:
            h_HP_warmup_tot = 0.0

        m_dot_LP_bypass_tot /= ts
        if m_dot_LP_bypass_tot > 0.0:
            h_LP_bypass_tot = h_LP_bypass_tot / m_dot_LP_bypass_tot
        else:
            h_LP_bypass_tot = 0.0

        m_dot_LP_warmup_tot /= ts
        if m_dot_LP_warmup_tot > 0.0:
            h_LP_warmup_tot = h_LP_warmup_tot / m_dot_LP_warmup_tot
        else:
            h_LP_warmup_tot = 0.0

        m_dot_AUX_DA_tot /= ts
        if m_dot_AUX_DA_tot > 0.0:
            h_AUX_DA_tot = h_AUX_DA_tot / m_dot_AUX_DA_tot
        else:
            h_AUX_DA_tot = 0.0

        m_dot_SS_drain_tot /= ts
        if m_dot_SS_drain_tot > 0.0:
            h_SS_drain_tot = h_SS_drain_tot / m_dot_SS_drain_tot
        else:
            h_SS_drain_tot = 0.0

        # SGT request update with convergence damping (replaces Fortran lines 2590–2602)
        m_dot_SGT_prev = self.m_dot_SGT.v  # previous value sent to steam drum
        m_dot_SGT_tot /= ts
        if abs(m_dot_SGT_prev - m_dot_SGT_tot) < 0.1:
            m_dot_SGT_tot = m_dot_SGT_prev  # no change — help convergence
        else:
            if m_dot_SGT_prev > m_dot_SGT_tot:
                m_dot_SGT_tot = m_dot_SGT_prev - abs(m_dot_SGT_prev - m_dot_SGT_tot) * 0.6
            else:
                m_dot_SGT_tot = m_dot_SGT_prev + abs(m_dot_SGT_prev - m_dot_SGT_tot) * 0.6

        # Condenser combined flow (rounded to 3 dp to aid convergence)
        m_dot_cond = (m_dot_LPT_exh_tot + m_dot_LP_bypass_tot
                      + m_dot_LP_warmup_tot + m_dot_HP_warmup_tot)
        m_dot_cond = round(m_dot_cond, 3)  # replaces Fortran float(nint(…*1000))/1000
        if m_dot_cond > 0.0:
            h_cond = ((m_dot_LPT_exh_tot * h_LPT_exh_tot
                       + m_dot_LP_bypass_tot * h_LP_bypass_tot
                       + m_dot_LP_warmup_tot * h_LP_warmup_tot
                       + m_dot_HP_warmup_tot * h_HP_warmup_tot)
                      / m_dot_cond)
            h_cond = float(round(h_cond))  # replaces Fortran float(nint(h_cond))
        else:
            m_dot_cond = 0.0
            h_cond = 0.0

        # Deaerator combined flow
        m_dot_DA = m_dot_AUX_DA_tot + m_dot_SS_drain_tot
        if m_dot_DA > 0.0:
            h_DA = ((m_dot_AUX_DA_tot * h_AUX_DA_tot
                     + m_dot_SS_drain_tot * h_SS_drain_tot) / m_dot_DA)
        else:
            h_DA = 0.0

        # ---- Final output assignments ----
        self.Turbine_ON_out.v = Turbine_ON
        self.m_dot_SGT.v = m_dot_SGT_tot
        self.m_dot_cond.v = m_dot_cond
        self.h_cond.v = h_cond
        self.m_dot_DA.v = m_dot_DA
        self.h_DA.v = h_DA
        self.W_dot_tot.v = W_dot_tot
        self.TRIP_TS.v = TRIP_TS
        self.m_dot_HP_bypass.v = m_dot_HP_bypass_tot
        self.m_dot_HP_AUX.v = m_dot_HP_AUX_tot
        self.m_dot_HP_drain.v = m_dot_HP_drain_tot
        self.m_dot_HP_warmup.v = m_dot_HP_warmup_tot
        self.m_dot_LP_AUX.v = m_dot_LP_AUX_tot
        self.m_dot_LP_bypass.v = m_dot_LP_bypass_tot
        self.m_dot_LP_drain.v = m_dot_LP_drain_tot
        self.m_dot_LP_warmup.v = m_dot_LP_warmup_tot
        self.m_dot_AUX_DA.v = m_dot_AUX_DA_tot

        # ---- Update dynamic piping state for is_converged block ----
        self._m_HPmain_prev = m_HPmain_prev
        self._P_HPmain_prev = P_HPmain_prev
        self._h_HPmain_prev = h_HPmain_prev
        self._T_HP_pipe = T_HP_pipe
        self._m_LPmain_prev = m_LPmain_prev
        self._P_LPmain_prev = P_LPmain_prev
        self._h_LPmain_prev = h_LPmain_prev
        self._T_LP_pipe = T_LP_pipe
        self._m_Aux_prev = m_Aux_prev
        self._P_Aux_prev = P_Aux_prev
        self._h_Aux_prev = h_Aux_prev
        self._T_Aux_pipe = T_Aux_pipe
