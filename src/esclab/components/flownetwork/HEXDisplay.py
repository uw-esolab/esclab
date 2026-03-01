"""HEX Display component model (Type 4102)."""

# Object: ESOL4102-HEX-Display
# Simulation Studio Model: ESOL4102-HEX-Display
#
# Author: Nick Edwards
# Date: 2025/05/13

import numpy as np

from esclab.simulate import Component


class HEXDisplay(Component):
    """
    TRNSYS Type 4102: ESOL4102-HEX-Display.

    Resolves which HEX ports to show on the top and bottom sides of a display
    diagram based on TES mode (Charging = 1, Discharging = 2, None = 0), and
    computes time derivatives of the four display temperatures using a
    4th-order backward finite-difference stencil (Fornberg, 1988).

    Parameters
    ----------
    (none)

    Inputs
    ------
    input_tes_mode    : float
        tes mode from ESOL 6034 (Charging = 1, Discharging = 2, None = 0)
    m_dot_shell_c_in  : float
        Mass flow rate of Molten salt entering the shell in the charging direction
    T_shell_c_in      : float
        Temperature of the Molten salt entering the shell side in the charging direction
    P_shell_c_in      : float
        Pressure of the Molten Salt entering the shell side in the charging direction
    m_dot_shell_d_in  : float
        Mass flow rate of Molten salt entering the shell side of the heat exchanger in the discharging direction
    T_shell_d_in      : float
        Temperature of the Molten salt entering the shell side in the discharging direction                               ^
    P_shell_d_in      : float
        Pressure of the Molten Salt entering the shell side of the heat exchanger in the discharging direction            | shell
    m_dot_tube_c_in   : float
        Mass Flow rate of HTF entering the tube side of the heat exchanger in the charging direction                      | tube
    T_tube_c_in       : float
        Temperature of the HTF entering the tube side of the HX in the charging direction                                 v
    P_tube_c_in       : float
        Pressure of the HTF entering the tube side of the HX in the charging direction
    m_dot_tube_d_in   : float
        Mass Flow rate of HTF entering the tube side of the HX in the discharging direction
    T_tube_d_in       : float
        Temperature of HTF entering the tube side of the HX in the discharging direction                           ^
    P_tube_d_in       : float
        Pressure of the HTF entering the tube side of the HX in the discharging direction                          | in
    m_dot_shell_c_out : float
        Mass flow rate leaving the shell side of the HX in the charging direction                                  | out
    T_shell_c_out     : float
        Temperature leaving the shell side of the HX in the charging direction                                     v
    P_shell_c_out     : float
        Pressure leaving the shell side of the HX in the charging direction
    m_dot_shell_d_out : float
        Mass flow rate leaving the shell side of the HX in the discharging direction
    T_shell_d_out     : float
        Temperature leaving the shell side of the HX in the discharging direction                                         ^
    P_shell_d_out     : float
        Pressure leaving the shell side of the HX in the discharging direction                                            | shell
    m_dot_tube_c_out  : float
        Mass flow rate leaving the tube side of the HX in the charging direction                                          | tube
    T_tube_c_out      : float
        Temperature leaving the tube side of the HX in the charging direction                                             v
    P_tube_c_out      : float
        Pressure leaving the tube side of the HX in the charging direction
    m_dot_tube_d_out  : float
        Mass flow rate leaving the tube side of the HX in the discharging direction
    T_tube_d_out      : float
        Temperature leaving the tube side of the HX in the discharging direction
    P_tube_d_out      : float
        Pressure leaving the tube side of the HX in the discharging direction

    Outputs
    -------
    m_dot_shell_top    : float
        mass flow of salt of shell on display's top side
    T_shell_top        : float
        temperature of salt of shell on display's top side
    P_shell_top        : float
        pressure of salt of shell on display's top side
    m_dot_tube_top     : float
        mass flow of HTF of tube on display's top side
    T_tube_top         : float
        temperature of HTF of tube on display's top side
    P_tube_top         : float
        pressure of HTF of tube on display's top side
    m_dot_shell_bottom : float
        mass flow of salt of shell on display's bottom side
    T_shell_bottom     : float
        temperature of salt of shell on display's bottom side
    P_shell_bottom     : float
        pressure of salt of shell on display's bottom side
    m_dot_tube_bottom  : float
        mass flow of HTF of tube on display's bottom side
    T_tube_bottom      : float
        temperature of HTF of tube on display's bottom side
    P_tube_bottom      : float
        pressure of HTF of tube on display's bottom side
    dT_dt_shell_top    : float
        time derivative of temperature of salt of shell on display's top side
    dT_dt_tube_top     : float
        time derivative of temperature of HTF of tube on display's top side
    dT_dt_shell_bottom : float
        time derivative of temperature of salt of shell on display's bottom side
    dT_dt_tube_bottom  : float
        time derivative of temperature of HTF of tube on display's bottom side
    """

    # *** Model Inputs ***
    # input_tes_mode: tes mode from ESOL 6034 (Charging = 1, Discharging = 2, None = 0)
    input_tes_mode    = Component.Input()
    # m_dot_shell_c_in: Mass flow rate of Molten salt entering the shell in the charging direction
    m_dot_shell_c_in  = Component.Input()
    # T_shell_c_in: Temperature of the Molten salt entering the shell side in the charging direction
    T_shell_c_in      = Component.Input()
    # P_shell_c_in: Pressure of the Molten Salt entering the shell side in the charging direction
    P_shell_c_in      = Component.Input()
    # m_dot_shell_d_in: Mass flow rate of Molten salt entering the shell side of the HX in the discharging direction
    m_dot_shell_d_in  = Component.Input()
    # T_shell_d_in: Temperature of the Molten salt entering the shell side in the discharging direction
    T_shell_d_in      = Component.Input()
    # P_shell_d_in: Pressure of the Molten Salt entering the shell side of the HX in the discharging direction
    P_shell_d_in      = Component.Input()
    # m_dot_tube_c_in: Mass Flow rate of HTF entering the tube side of the HX in the charging direction
    m_dot_tube_c_in   = Component.Input()
    # T_tube_c_in: Temperature of the HTF entering the tube side of the HX in the charging direction
    T_tube_c_in       = Component.Input()
    # P_tube_c_in: Pressure of the HTF entering the tube side of the HX in the charging direction
    P_tube_c_in       = Component.Input()
    # m_dot_tube_d_in: Mass Flow rate of HTF entering the tube side of the HX in the discharging direction
    m_dot_tube_d_in   = Component.Input()
    # T_tube_d_in: Temperature of HTF entering the tube side of the HX in the discharging direction
    T_tube_d_in       = Component.Input()
    # P_tube_d_in: Pressure of the HTF entering the tube side of the HX in the discharging direction
    P_tube_d_in       = Component.Input()
    # m_dot_shell_c_out: Mass flow rate leaving the shell side of the HX in the charging direction
    m_dot_shell_c_out = Component.Input()
    # T_shell_c_out: Temperature leaving the shell side of the HX in the charging direction
    T_shell_c_out     = Component.Input()
    # P_shell_c_out: Pressure leaving the shell side of the HX in the charging direction
    P_shell_c_out     = Component.Input()
    # m_dot_shell_d_out: Mass flow rate leaving the shell side of the HX in the discharging direction
    m_dot_shell_d_out = Component.Input()
    # T_shell_d_out: Temperature leaving the shell side of the HX in the discharging direction
    T_shell_d_out     = Component.Input()
    # P_shell_d_out: Pressure leaving the shell side of the HX in the discharging direction
    P_shell_d_out     = Component.Input()
    # m_dot_tube_c_out: Mass flow rate leaving the tube side of the HX in the charging direction
    m_dot_tube_c_out  = Component.Input()
    # T_tube_c_out: Temperature leaving the tube side of the HX in the charging direction
    T_tube_c_out      = Component.Input()
    # P_tube_c_out: Pressure leaving the tube side of the HX in the charging direction
    P_tube_c_out      = Component.Input()
    # m_dot_tube_d_out: Mass flow rate leaving the tube side of the HX in the discharging direction
    m_dot_tube_d_out  = Component.Input()
    # T_tube_d_out: Temperature leaving the tube side of the HX in the discharging direction
    T_tube_d_out      = Component.Input()
    # P_tube_d_out: Pressure leaving the tube side of the HX in the discharging direction
    P_tube_d_out      = Component.Input()

    # *** Model Outputs ***
    m_dot_shell_top    = Component.Output()  # mass flow of salt of shell on display's top side
    T_shell_top        = Component.Output()  # temperature of salt of shell on display's top side
    P_shell_top        = Component.Output()  # pressure of salt of shell on display's top side
    m_dot_tube_top     = Component.Output()  # mass flow of HTF of tube on display's top side
    T_tube_top         = Component.Output()  # temperature of HTF of tube on display's top side
    P_tube_top         = Component.Output()  # pressure of HTF of tube on display's top side
    m_dot_shell_bottom = Component.Output()  # mass flow of salt of shell on display's bottom side
    T_shell_bottom     = Component.Output()  # temperature of salt of shell on display's bottom side
    P_shell_bottom     = Component.Output()  # pressure of salt of shell on display's bottom side
    m_dot_tube_bottom  = Component.Output()  # mass flow of HTF of tube on display's bottom side
    T_tube_bottom      = Component.Output()  # temperature of HTF of tube on display's bottom side
    P_tube_bottom      = Component.Output()  # pressure of HTF of tube on display's bottom side
    dT_dt_shell_top    = Component.Output()  # time derivative of temperature of salt of shell on display's top side
    dT_dt_tube_top     = Component.Output()  # time derivative of temperature of HTF of tube on display's top side
    dT_dt_shell_bottom = Component.Output()  # time derivative of temperature of salt of shell on display's bottom side
    dT_dt_tube_bottom  = Component.Output()  # time derivative of temperature of HTF of tube on display's bottom side

    # TODO-NEEDS CONVERSION REVIEW: dynamic array storage
    # SetNumberStoredVariables(0, 21) - dynamic array of 21 elements (0-indexed here).
    # Fortran 1-based index -> Python 0-based index mapping:
    #  0  (index  1) - tes_mode
    #  1  (index  2) - T_shell_top -1 (most recent previous timestep)
    #  2  (index  3) - T_shell_top -2
    #  3  (index  4) - T_shell_top -3
    #  4  (index  5) - T_shell_top -4
    #  5  (index  6) - T_tube_top -1
    #  6  (index  7) - T_tube_top -2
    #  7  (index  8) - T_tube_top -3
    #  8  (index  9) - T_tube_top -4
    #  9  (index 10) - T_shell_bottom -1
    # 10  (index 11) - T_shell_bottom -2
    # 11  (index 12) - T_shell_bottom -3
    # 12  (index 13) - T_shell_bottom -4
    # 13  (index 14) - T_tube_bottom -1
    # 14  (index 15) - T_tube_bottom -2
    # 15  (index 16) - T_tube_bottom -3
    # 16  (index 17) - T_tube_bottom -4
    # 17-20 (indices 18-21) - initialized to 0, unused
    # for 2 second timesteps, need 4 previous values along with current value
    _dynamic: np.ndarray  # shape (21,)

    def initialize(self):
        # Set temperature storage to 0
        # SetDynamicArrayInitialValue for indices 2..21 -> 0.0
        self._dynamic = np.zeros(21)

    def calculate(self):
        # -----------------------------------------------------------------------------------------------------------------------
        # Do All of the First Timestep Manipulations Here - There Are No Iterations at the Initial Time
        if self.model.is_first_step:
            # if tes mode is neutral, 0, then assume charging, 1
            if self.input_tes_mode.v < 2.0:
                tes_mode = 1.0
                self.m_dot_shell_top.v    = self.m_dot_shell_c_in.v
                self.T_shell_top.v        = self.T_shell_c_in.v
                self.P_shell_top.v        = self.P_shell_c_in.v
                self.m_dot_tube_top.v     = self.m_dot_tube_c_out.v
                self.T_tube_top.v         = self.T_tube_c_out.v
                self.P_tube_top.v         = self.P_tube_c_out.v
                self.m_dot_shell_bottom.v = self.m_dot_shell_c_out.v
                self.T_shell_bottom.v     = self.T_shell_c_out.v
                self.P_shell_bottom.v     = self.P_shell_c_out.v
                self.m_dot_tube_bottom.v  = self.m_dot_tube_c_in.v
                self.T_tube_bottom.v      = self.T_tube_c_in.v
                self.P_tube_bottom.v      = self.P_tube_c_in.v
            else:
                tes_mode = 2.0
                self.m_dot_shell_top.v    = self.m_dot_shell_d_out.v
                self.T_shell_top.v        = self.T_shell_d_out.v
                self.P_shell_top.v        = self.P_shell_d_out.v
                self.m_dot_tube_top.v     = self.m_dot_tube_d_in.v
                self.T_tube_top.v         = self.T_tube_d_in.v
                self.P_tube_top.v         = self.P_tube_d_in.v
                self.m_dot_shell_bottom.v = self.m_dot_shell_d_in.v
                self.T_shell_bottom.v     = self.T_shell_d_in.v
                self.P_shell_bottom.v     = self.P_shell_d_in.v
                self.m_dot_tube_bottom.v  = self.m_dot_tube_d_out.v
                self.T_tube_bottom.v      = self.T_tube_d_out.v
                self.P_tube_bottom.v      = self.P_tube_d_out.v

            # Set the Initial Values of the Outputs (#,Value)
            self.dT_dt_shell_top.v    = 0.0  # time derivative of temperature of salt of shell on display's top side
            self.dT_dt_tube_top.v     = 0.0  # time derivative of temperature of HTF of tube on display's top side
            self.dT_dt_shell_bottom.v = 0.0  # time derivative of temperature of salt of shell on display's bottom side
            self.dT_dt_tube_bottom.v  = 0.0  # time derivative of temperature of HTF of tube on display's bottom side

            # Set Initial Values of Dynamic Storage
            # TODO-NEEDS CONVERSION REVIEW: dynamic array storage
            self._dynamic[0] = tes_mode  # SetDynamicArrayValueThisIteration(1, tes_mode) - Store tes mode

            return

        # -----------------------------------------------------------------------------------------------------------------------
        # Read the Inputs (already available as self.<input>.v)

        # if tes is neutral then display previous mode
        if self.input_tes_mode.v == 0.0:
            # TODO-NEEDS CONVERSION REVIEW: dynamic array storage
            tes_mode = self._dynamic[0]  # getDynamicArrayValueLastTimestep(1)
        else:
            tes_mode = self.input_tes_mode.v

        if tes_mode == 1.0:
            self.m_dot_shell_top.v    = self.m_dot_shell_c_in.v
            self.T_shell_top.v        = self.T_shell_c_in.v
            self.P_shell_top.v        = self.P_shell_c_in.v
            self.m_dot_tube_top.v     = self.m_dot_tube_c_out.v
            self.T_tube_top.v         = self.T_tube_c_out.v
            self.P_tube_top.v         = self.P_tube_c_out.v
            self.m_dot_shell_bottom.v = self.m_dot_shell_c_out.v
            self.T_shell_bottom.v     = self.T_shell_c_out.v
            self.P_shell_bottom.v     = self.P_shell_c_out.v
            self.m_dot_tube_bottom.v  = self.m_dot_tube_c_in.v
            self.T_tube_bottom.v      = self.T_tube_c_in.v
            self.P_tube_bottom.v      = self.P_tube_c_in.v
        else:
            self.m_dot_shell_top.v    = self.m_dot_shell_d_out.v
            self.T_shell_top.v        = self.T_shell_d_out.v
            self.P_shell_top.v        = self.P_shell_d_out.v
            self.m_dot_tube_top.v     = self.m_dot_tube_d_in.v
            self.T_tube_top.v         = self.T_tube_d_in.v
            self.P_tube_top.v         = self.P_tube_d_in.v
            self.m_dot_shell_bottom.v = self.m_dot_shell_d_in.v
            self.T_shell_bottom.v     = self.T_shell_d_in.v
            self.P_shell_bottom.v     = self.P_shell_d_in.v
            self.m_dot_tube_bottom.v  = self.m_dot_tube_d_out.v
            self.T_tube_bottom.v      = self.T_tube_d_out.v
            self.P_tube_bottom.v      = self.P_tube_d_out.v

        # derivative coefficients, from Table 3 in below (free online)
        # B. Fornberg, "Generation of Finite Difference Formulas on Arbitrarily Spaced Grids," Mathematics of Computation, vol. 51, no. 184, pp. 699-706, 1988, doi: 10.2307/2008770.
        # TODO-NEEDS CONVERSION REVIEW: getTimeStepNumber() mapped to a step_number computed from model.time and model.settings; verify start_time offset is correct
        step_number = int((self.model.time - self.model.settings.start_time) / self.model.settings.timestep) + 1
        if step_number > 4:
            # TODO-NEEDS CONVERSION REVIEW: dynamic array storage
            T_1 = self._dynamic[1]   # getDynamicArrayValueLastTimestep(2)  - T_shell_top -1
            T_2 = self._dynamic[2]   # getDynamicArrayValueLastTimestep(3)  - T_shell_top -2
            T_3 = self._dynamic[3]   # getDynamicArrayValueLastTimestep(4)  - T_shell_top -3
            T_4 = self._dynamic[4]   # getDynamicArrayValueLastTimestep(5)  - T_shell_top -4

            self.dT_dt_shell_top.v = (25.0 * self.T_shell_top.v / 12.0 - 4.0 * T_1 + 3.0 * T_2 - 4.0 * T_3 / 3.0 + T_4 / 4.0) / (3600.0 * self.model.settings.timestep)

            T_1 = self._dynamic[5]   # getDynamicArrayValueLastTimestep(6)  - T_tube_top -1
            T_2 = self._dynamic[6]   # getDynamicArrayValueLastTimestep(7)  - T_tube_top -2
            T_3 = self._dynamic[7]   # getDynamicArrayValueLastTimestep(8)  - T_tube_top -3
            T_4 = self._dynamic[8]   # getDynamicArrayValueLastTimestep(9)  - T_tube_top -4

            self.dT_dt_tube_top.v = (25.0 * self.T_tube_top.v / 12.0 - 4.0 * T_1 + 3.0 * T_2 - 4.0 * T_3 / 3.0 + T_4 / 4.0) / (3600.0 * self.model.settings.timestep)

            T_1 = self._dynamic[9]    # getDynamicArrayValueLastTimestep(10) - T_shell_bottom -1
            T_2 = self._dynamic[10]   # getDynamicArrayValueLastTimestep(11) - T_shell_bottom -2
            T_3 = self._dynamic[11]   # getDynamicArrayValueLastTimestep(12) - T_shell_bottom -3
            T_4 = self._dynamic[12]   # getDynamicArrayValueLastTimestep(13) - T_shell_bottom -4

            self.dT_dt_shell_bottom.v = (25.0 * self.T_shell_bottom.v / 12.0 - 4.0 * T_1 + 3.0 * T_2 - 4.0 * T_3 / 3.0 + T_4 / 4.0) / (3600.0 * self.model.settings.timestep)

            T_1 = self._dynamic[13]   # getDynamicArrayValueLastTimestep(14) - T_tube_bottom -1
            T_2 = self._dynamic[14]   # getDynamicArrayValueLastTimestep(15) - T_tube_bottom -2
            T_3 = self._dynamic[15]   # getDynamicArrayValueLastTimestep(16) - T_tube_bottom -3
            T_4 = self._dynamic[16]   # getDynamicArrayValueLastTimestep(17) - T_tube_bottom -4

            self.dT_dt_tube_bottom.v = (25.0 * self.T_tube_bottom.v / 12.0 - 4.0 * T_1 + 3.0 * T_2 - 4.0 * T_3 / 3.0 + T_4 / 4.0) / (3600.0 * self.model.settings.timestep)
        else:
            self.dT_dt_shell_top.v    = 0.0
            self.dT_dt_tube_top.v     = 0.0
            self.dT_dt_shell_bottom.v = 0.0
            self.dT_dt_tube_bottom.v  = 0.0

        # -----------------------------------------------------------------------------------------------------------------------
        # Update dynamic storage at end of timestep
        if self.model.is_converged:
            # Check if dynamic storage needs to be updated
            # TODO-NEEDS CONVERSION REVIEW: dynamic array storage
            if tes_mode != self._dynamic[0]:  # getDynamicArrayValueLastTimestep(1)
                self._dynamic[0] = tes_mode  # SetDynamicArrayValueThisIteration(1, tes_mode)

            # shuffle temperature storage

            #  5 - T_shell_top -4
            self._dynamic[4] = self._dynamic[3]    # T_shuffle = getDynamicArrayValueLastTimestep(4); SetDynamicArrayValueThisIteration(5, T_shuffle)
            #  4 - T_shell_top -3
            self._dynamic[3] = self._dynamic[2]    # T_shuffle = getDynamicArrayValueLastTimestep(3); SetDynamicArrayValueThisIteration(4, T_shuffle)
            #  3 - T_shell_top -2
            self._dynamic[2] = self._dynamic[1]    # T_shuffle = getDynamicArrayValueLastTimestep(2); SetDynamicArrayValueThisIteration(3, T_shuffle)
            #  2 - T_shell_top -1
            self._dynamic[1] = self.T_shell_top.v  # SetDynamicArrayValueThisIteration(2, T_shell_top)

            #  9 - T_tube_top -4
            self._dynamic[8] = self._dynamic[7]    # T_shuffle = getDynamicArrayValueLastTimestep(8); SetDynamicArrayValueThisIteration(9, T_shuffle)
            #  8 - T_tube_top -3
            self._dynamic[7] = self._dynamic[6]    # T_shuffle = getDynamicArrayValueLastTimestep(7); SetDynamicArrayValueThisIteration(8, T_shuffle)
            #  7 - T_tube_top -2
            self._dynamic[6] = self._dynamic[5]    # T_shuffle = getDynamicArrayValueLastTimestep(6); SetDynamicArrayValueThisIteration(7, T_shuffle)
            #  6 - T_tube_top -1
            self._dynamic[5] = self.T_tube_top.v   # SetDynamicArrayValueThisIteration(6, T_tube_top)

            # 13 - T_shell_bottom -4
            self._dynamic[12] = self._dynamic[11]     # T_shuffle = getDynamicArrayValueLastTimestep(12); SetDynamicArrayValueThisIteration(13, T_shuffle)
            # 12 - T_shell_bottom -3
            self._dynamic[11] = self._dynamic[10]     # T_shuffle = getDynamicArrayValueLastTimestep(11); SetDynamicArrayValueThisIteration(12, T_shuffle)
            # 11 - T_shell_bottom -2
            self._dynamic[10] = self._dynamic[9]      # T_shuffle = getDynamicArrayValueLastTimestep(10); SetDynamicArrayValueThisIteration(11, T_shuffle)
            # 10 - T_shell_bottom -1
            self._dynamic[9] = self.T_shell_bottom.v  # SetDynamicArrayValueThisIteration(10, T_shell_bottom)

            # 17 - T_tube_bottom -4
            self._dynamic[16] = self._dynamic[15]     # T_shuffle = getDynamicArrayValueLastTimestep(16); SetDynamicArrayValueThisIteration(17, T_shuffle)
            # 16 - T_tube_bottom -3
            self._dynamic[15] = self._dynamic[14]     # T_shuffle = getDynamicArrayValueLastTimestep(15); SetDynamicArrayValueThisIteration(16, T_shuffle)
            # 15 - T_tube_bottom -2
            self._dynamic[14] = self._dynamic[13]     # T_shuffle = getDynamicArrayValueLastTimestep(14); SetDynamicArrayValueThisIteration(15, T_shuffle)
            # 14 - T_tube_bottom -1
            self._dynamic[13] = self.T_tube_bottom.v  # SetDynamicArrayValueThisIteration(14, T_tube_bottom)
