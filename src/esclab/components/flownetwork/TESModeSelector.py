"""TES mode selector component model (Type 6034)."""

import warnings

from esclab.simulate import Component


class TESModeSelector(Component):
    """
    Object: ESOL6034-TESModes
    Simulation Studio Model: TESModes

    Date: January 15, 2025
    last modified: January 15, 2025

    Determines the TES operating mode (Charging=1, Discharging=2, Inactive=0)
    based on pump power commands and valve positions. Guards against simultaneous
    activation of conflicting modes.

    Mode definitions:
        0 = TES not running
        1 = Charging mode
        2 = Discharging mode

    Inputs
    ------
    CP_Power   : float  Power Input for TES Charging Pump
    DP_Power   : float  Power Input for TES Discharging Pump
    HTF_CVP_i  : float  HTF Charging Valve Position, input requested by operators
    HTF_CVP_a  : float  HTF Charging Valve Position, current based on valve speed
    HTF_DVP_i  : float  HTF Discharging Valve, input requested by operators
    HTF_DVP_a  : float  HTF Discharging Valve, current valve position based on valve speed

    Outputs
    -------
    mode         : float  Current TES operating mode (Charging=1, Discharging=2, None=0)
    CP_Power_out : float  ChargingPump-Power - Input to Charging Pump
    DP_Power_out : float  DischargingPump-Power - Input to Discharging Pump
    HTF_DVP_out  : float  HTFDischargingValve-Position
    HTF_CVP_out  : float  HTFChargingValve-Position
    """

    #    VARIABLES
    # (no parameters declared in the Fortran; all logic is input/output-driven)

    #    INPUTS
    CP_Power  = Component.Input()   # Power Input for TES Charging Pump
    DP_Power  = Component.Input()   # Power Input for TES Discharging Pump
    HTF_CVP_i = Component.Input()   # HTF Charging Valve Position, input requested by operators
    HTF_CVP_a = Component.Input()   # HTF Charging Valve Position, current based on valve speed
    HTF_DVP_i = Component.Input()   # HTF Discharging Valve, input requested by operators
    HTF_DVP_a = Component.Input()   # HTF Discharging Valve, current valve position based on valve speed

    #    OUTPUTS
    mode         = Component.Output()  # Mode (Charging = 1, Discharging = 2, None = 0)
    CP_Power_out = Component.Output()  # ChargingPump-Power
    DP_Power_out = Component.Output()  # DischargingPump-Power
    # TODO-NEEDS CONVERSION REVIEW: At start-time the Fortran assigns output 4=HTF_CVP_i and
    # output 5=HTF_DVP_i, but in the main timestep body output 4=HTF_DVP_i and output 5=HTF_CVP_i.
    # The main-body assignments (below) are used as the canonical mapping.
    HTF_DVP_out  = Component.Output()  # HTFDischargingValve-Position (output 4 per main body)
    HTF_CVP_out  = Component.Output()  # HTFChargingValve-Position    (output 5 per main body)

    def presim_setup(self, **kwargs):
        # Set initial outputs to zero before start-time logic runs
        self.mode.v = 0.0
        self.CP_Power_out.v = 0.0
        self.DP_Power_out.v = 0.0
        self.HTF_DVP_out.v = 0.0
        self.HTF_CVP_out.v = 0.0

    def calculate(self):
        if self.model.is_first_step:
            # Do All of the First Timestep Manipulations Here
            CP_Power  = self.CP_Power.v
            DP_Power  = self.DP_Power.v
            HTF_CVP_i = self.HTF_CVP_i.v
            HTF_DVP_i = self.HTF_DVP_i.v

            # Select mode that is starting off
            # mode 0 = nothing running, mode 1 = charging, mode 2 = discharging
            if HTF_CVP_i > 0.0:  # charging mode established
                mode = 1.0
                # ensure discharging pump and HTF discharging valve are closed
                if DP_Power > 0.0:
                    warnings.warn(
                        'Multiple TES modes input at beginning of simulation, '
                        'charging mode defaulted to. Discharging Pump Turned off'
                    )
                    DP_Power = 0.0
                if HTF_DVP_i > 0.0:
                    warnings.warn(
                        'Multiple TES modes input at beginning of simulation, '
                        'charging mode defaulted to. HTF Discharging valve closed'
                    )
                    HTF_DVP_i = 0.0

            elif HTF_DVP_i > 0.0:  # discharging mode established
                mode = 2.0
                # ensure charging pump and HTF charging valve are closed
                if CP_Power > 0.0:
                    warnings.warn(
                        'Multiple TES modes input at beginning of simulation, '
                        'discharging mode defaulted to. Charging Pump turned off'
                    )
                    CP_Power = 0.0

            else:  # No HTF flow through TES pathway, check pumps for modes
                if CP_Power > 0.0:
                    mode = 1.0  # charging mode established
                    if DP_Power > 0.0:
                        warnings.warn(
                            'Multiple TES modes input at begining of simulation, '
                            'charging mode selected. Discharging pump turned off'
                        )
                        DP_Power = 0.0
                elif DP_Power > 0.0:
                    mode = 2.0  # discharging mode established
                else:
                    mode = 0.0  # TES is not running

            # Set the Initial Values of the Outputs (#,Value)
            # TODO-NEEDS CONVERSION REVIEW: At start-time the Fortran sets output 4=HTF_CVP_i and
            # output 5=HTF_DVP_i (swapped from main body). Matching start-time Fortran assignment here.
            self.mode.v         = mode       # Mode (Charging = 1)
            self.CP_Power_out.v = CP_Power   # ChargingPump-Power
            self.DP_Power_out.v = DP_Power   # DischargingPump-Power
            self.HTF_DVP_out.v  = HTF_CVP_i  # output 4 at start time = HTF_CVP_i
            self.HTF_CVP_out.v  = HTF_DVP_i  # output 5 at start time = HTF_DVP_i
            return

        # Modes can only switch at first timestep of iteration
        if self.model.timestep_iteration == 0:
            # Current Mode Running (1=Charging, 2=Discharging, 0=Inactive)
            mode     = self.mode.v
            CP_Power  = self.CP_Power.v
            DP_Power  = self.DP_Power.v
            HTF_CVP_i = self.HTF_CVP_i.v
            HTF_CVP_a = self.HTF_CVP_a.v
            HTF_DVP_i = self.HTF_DVP_i.v
            HTF_DVP_a = self.HTF_DVP_a.v

            if mode == 1.0:  # Current in Charging mode
                # check that charging pump and HTF charging valve position are still on/open
                if (CP_Power == 1.0) or (HTF_CVP_a > 0.0):
                    mode = 1.0  # Still in charging mode, ensure operator does not turn on discharging pump or HTF discharging valve
                    if DP_Power == 1.0:
                        warnings.warn(
                            'Cannot turn on Discharging Pump until charging pump is turned off '
                            'and HTF charging valve is fully closed'
                        )
                        DP_Power = 0.0
                    if HTF_DVP_i > 0.0:
                        warnings.warn(
                            'Cannot open HTF discharging valve while in charging mode, '
                            'turn off charging pump and fully close HTF charging valve'
                        )
                        HTF_DVP_i = 0.0
                else:  # charging pump and HTF CVP are now off, charging mode no longer activated, check if discharging pump/HTF discharging CV is open
                    if (HTF_DVP_i > 0.0) or (DP_Power == 1.0):  # now in discharging mode
                        mode = 2.0
                    else:  # TES is not in use
                        mode = 0.0

            elif mode == 2.0:  # Current mode is Discharging mode
                if (HTF_DVP_a > 0.0) or (DP_Power == 1.0):  # Still in Discharging Mode
                    mode = 2.0
                    # ensure operator does not turn on charging pump or HTF Charging valve while in discharging mode
                    if CP_Power == 1.0:
                        warnings.warn(
                            'Cannot turn on TES charging pump while in discharging mode, '
                            'fully close HTF discharging valve and turn off Discharging Pump'
                        )
                        CP_Power = 0.0
                    if HTF_CVP_i > 0.0:
                        warnings.warn(
                            'Cannot open HTF Charging Valve while in discharging mode, '
                            'fully close HTF discharging valve and turn off Discharging pump'
                        )
                        HTF_CVP_i = 0.0
                else:  # discharging pump is off and HTF discharging valve is fully closed, OK to switch modes
                    if (HTF_CVP_i > 0.0) or (CP_Power == 1.0):
                        mode = 1.0  # Switch to charging mode
                    else:
                        mode = 0.0  # TES is not running

            else:  # TES is currently not in use check if mode is turned on
                if (CP_Power == 1.0) or (HTF_CVP_i > 0.0):  # charging mode activated
                    mode = 1.0
                    # check that discharging valve or discharging pump were not also turned on
                    if DP_Power == 1.0:
                        warnings.warn(
                            'Tried to turn on two modes at once, charging is defaulted to. '
                            'Discharging Pump is turned off'
                        )
                        DP_Power = 0.0
                    if HTF_DVP_i > 0.0:
                        warnings.warn(
                            'Tried to turn on two modes at once, charging is defaulted to. '
                            'HTF discharging valve is closed'
                        )
                        HTF_DVP_i = 0.0
                elif (DP_Power == 1.0) or (HTF_DVP_i > 0.0):  # discharging mode activated
                    mode = 2.0
                else:  # all pumps and valves are still off/closed, do not change mode
                    mode = 0.0

            # Call NEW Outputs
            self.mode.v         = mode       # Mode (Charging = 1, Discharging = 2, None = 0)
            self.CP_Power_out.v = CP_Power   # ChargingPump-Power - Input to Charging Pump
            self.DP_Power_out.v = DP_Power   # DischargingPump-Power - Input to Discharging Pump
            self.HTF_DVP_out.v  = HTF_DVP_i  # HTFDischargingValve-Position  - Input to Charging HTF Valve
            self.HTF_CVP_out.v  = HTF_CVP_i  # HTFChargingValve-Position    - Input to Discharging HTF Valve
