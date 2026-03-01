"""Shell-and-tube heat exchanger component model (Type 4101)."""

# Object: HEX-LL
# Simulation Studio Model: ESOL4101-HEX
#
# Author: Sergio Alcalde Morales
# Date:   August 05, 2024
# Last modified: August 05, 2024
# Converted by: GitHub Copilot, March 01, 2026

import numpy as np

from esclab.simulate import Component
from esclab.components.esol_properties import Incompressible as Inc


class HeatExchanger(Component):
    """
    TRNSYS Type 4101: ESOL4101-HEX
    Object: HEX-LL

    Shell-and-tube liquid-liquid heat exchanger supporting both charging and
    discharging operating modes.  A nodal temperature and enthalpy profile
    (``n_nodes = n_baffles * shell_p * n_HEX + 1``) is maintained across
    time-steps in instance arrays that mimic TRNSYS dynamic stored variables.

    At the first time-step an effectiveness-NTU (e-NTU) steady-state
    solution initialises the profile.  Every subsequent call drives a
    zero-dimensional shell/tube energy balance (``ZeroD_Eq``), a one-
    dimensional shell-side thermal-hydraulics solve (``shell``), and a
    one-dimensional tube-side incompressible solve (``tube_oneD_inc``).

    The component handles dual-direction operation: when
    ``m_dot_tube_c > m_dot_tube_d`` the HX is in *charging* mode
    (``mode_iter = 1``); otherwise it is in *discharging* mode
    (``mode_iter = 0``).  When the operating mode changes between
    time-steps the stored profiles are reversed to account for the
    reversal of flow direction.

    Parameters
    ----------
    n_baffles : int
        Number of baffles in one shell pass [-]  ``[0; +Inf]``
    n_HEX : int
        Number of heat exchangers in series [-]  ``[1; +Inf]``
    tube_p : float
        Number of tube passes in a heat exchanger [-]  ``[1; +Inf]``
    shell_p : int
        Number of shell passes in a heat exchanger [-]  ``[1; +Inf]``
    D_shell : float
        Inner diameter of the shell [m]  ``[0; +Inf]``
    L_shell : float
        Length of the shell [m]  ``[0; +Inf]``
    r_in : float
        Inner radius of the tube [m]  ``[0; +Inf]``
    r_out : float
        Outer radius of the tube [m]  ``[0; +Inf]``
    L_baffle : float
        Length of the baffle [m]  ``[0; +Inf]``
    th_baffle : float
        Thickness of the baffle relative to baffle spacing [-]  ``[0; 1]``
    th_ins : float
        Thickness of the insulation surrounding the HX [m]  ``[0; +Inf]``
    S_T : float
        Pitch between tubes in the transverse direction [m]  ``[0; +Inf]``
    S_L : float
        Pitch between tubes in the longitudinal direction [m]  ``[0; +Inf]``
    n_tubes : float
        Number of tubes bundled together [-]  ``[1; +Inf]``
    config : bool
        Heat exchanger configuration flag ``[0; 1]``
    emiss : float
        Emissivity to the environment ``[0; 1]``
    k_ins : float
        Thermal conductivity of the insulation [W/m·K]  ``[0; +Inf]``
    k_steel : float
        Thermal conductivity of the shell steel [W/m·K]  ``[0; +Inf]``
    fluid_tube : str
        Fluid identifier for the tube-side fluid ``[0; 1000]``
    fluid_shell : str
        Fluid identifier for the shell-side fluid ``[0; 1000]``

    Inputs
    ------
    v_out : float
        Velocity of the air outside of the heat exchanger [m/s]  ``[0; 300]``
    Tenv : float
        Temperature outside of the heat exchanger [K]  ``[250; 350]``
    m_dot_shell_c : float
        Mass flow rate of molten salt entering the shell
        (charging direction) [kg/s]
    T_shell_in_c : float
        Temperature of the molten salt entering the shell side
        (charging direction) [K]
    P_shell_in_c : float
        Pressure of the molten salt entering the shell side
        (charging direction) [Pa]
    m_dot_shell_d : float
        Mass flow rate of molten salt entering the shell side
        (discharging direction) [kg/s]
    T_shell_in_d : float
        Temperature of the molten salt entering the shell side
        (discharging direction) [K]
    P_shell_in_d : float
        Pressure of the molten salt entering the shell side
        (discharging direction) [Pa]
    m_dot_tube_c : float
        Mass flow rate of HTF entering the tube side
        (charging direction) [kg/s]
    T_tube_in_c : float
        Temperature of the HTF entering the tube side
        (charging direction) [K]
    P_tube_in_c : float
        Pressure of the HTF entering the tube side
        (charging direction) [Pa]
    m_dot_tube_d : float
        Mass flow rate of HTF entering the tube side
        (discharging direction) [kg/s]
    T_tube_in_d : float
        Temperature of the HTF entering the tube side
        (discharging direction) [K]
    P_tube_in_d : float
        Pressure of HTF entering the tube side
        (discharging direction) [Pa]
    mass_counter_c : float
        Mass counter for the HTF expansion system (charging direction) [kg]
    mass_counter_d : float
        Mass counter for the HTF expansion system (discharging direction) [kg]

    Outputs
    -------
    mode_iter : float
        Operating mode: 1 = charging, 0 = discharging [-]
    m_dot_shell_c_out : float
        Mass flow rate leaving the shell side (charging direction) [kg/s]
    T_shell_out_c : float
        Temperature leaving the shell side (charging direction) [K]
    P_shell_out_c : float
        Pressure leaving the shell side (charging direction) [Pa]
    m_dot_shell_d_out : float
        Mass flow rate leaving the shell side (discharging direction) [kg/s]
    T_shell_out_d : float
        Temperature leaving the shell side (discharging direction) [K]
    P_shell_out_d : float
        Pressure leaving the shell side (discharging direction) [Pa]
    m_dot_tube_c_out : float
        Mass flow rate leaving the tube side (charging direction) [kg/s]
    T_tube_out_c : float
        Temperature leaving the tube side (charging direction) [K]
    P_tube_out_c : float
        Pressure leaving the tube side (charging direction) [Pa]
    m_dot_tube_d_out : float
        Mass flow rate leaving the tube side (discharging direction) [kg/s]
    T_tube_out_d : float
        Temperature leaving the tube side (discharging direction) [K]
    P_tube_out_d : float
        Pressure leaving the tube side (discharging direction) [Pa]
    mass_counter_c_out : float
        HTF mass counter for the expansion system (charging direction) [kg]
    mass_counter_d_out : float
        HTF mass counter for the expansion system (discharging direction) [kg]
    """

    # *** Model Parameters ***
    # n_baffles: Number of baffles in one shell pass of a heat exchanger [0;+Inf]
    n_baffles = Component.Parameter()
    # n_HEX: Number of Heat Exchangers in series [1;+Inf]
    n_HEX = Component.Parameter()
    # tube_p: Number of tube passes in a heat exchanger [1;+Inf]
    tube_p = Component.Parameter()
    # shell_p: Number of shell passes in a heat exchanger [1;+Inf]
    shell_p = Component.Parameter()
    # D_shell: Inner diameter of the shell [m] [0;+Inf]
    D_shell = Component.Parameter()
    # L_shell: Length of the shell [m] [0;+Inf]
    L_shell = Component.Parameter()
    # r_in: Inner radius of the tube [m] [0;+Inf]
    r_in = Component.Parameter()
    # r_out: Outer radius of the tube [m] [0;+Inf]
    r_out = Component.Parameter()
    # L_baffle: Length of the baffle [m] [0;+Inf]
    L_baffle = Component.Parameter()
    # th_baffle: Thickness of the baffle [m] [0;1]
    th_baffle = Component.Parameter()
    # th_ins: Thickness of the insulation surrounding the HX [m] [0;+Inf]
    th_ins = Component.Parameter()
    # S_T: Pitch between the tubes in the transverse direction (against flow) [m] [0;+Inf]
    S_T = Component.Parameter()
    # S_L: Pitch between the tubes in the longitudinal direction (with flow) [m] [0;+Inf]
    S_L = Component.Parameter()
    # n_tubes: Number of tubes bundled together [-] [1;+Inf]
    n_tubes = Component.Parameter()
    # config: Heat Exchanger configuration [0;1]
    config = Component.Parameter()
    # emiss: Emissivity to the environment [0;1]
    emiss = Component.Parameter()
    # k_ins: Thermal conductivity of the insulation [W/m.K] [0;+Inf]
    k_ins = Component.Parameter()
    # k_steel: Thermal conductivity of the shell [W/m.K] [0;+Inf]
    k_steel = Component.Parameter()
    # fluid_tube: Fluid ID through the tube side of the HX (HTF = 40) [0;1000]
    fluid_tube = Component.Parameter()
    # fluid_shell: Fluid ID through the shell side of the HX (Molten Salt = 17) [0;1000]
    fluid_shell = Component.Parameter()

    # *** Model Inputs ***
    # v_out: Velocity of the air outside of the heat exchanger [m/s] [0;300]
    v_out = Component.Input()
    # Tenv: Temperature outside of the heat exchanger [K] [250;350]
    Tenv = Component.Input()
    # m_dot_shell_c: Mass flow rate of Molten salt entering the shell (charging direction) [kg/s]
    m_dot_shell_c = Component.Input()
    # T_shell_in_c: Temperature of the Molten salt entering the shell side (charging direction) [K]
    T_shell_in_c = Component.Input()
    # P_shell_in_c: Pressure of the Molten Salt entering the shell side (charging direction) [Pa]
    P_shell_in_c = Component.Input()
    # m_dot_shell_d: Mass flow rate of Molten salt entering the shell side (discharging direction) [kg/s]
    m_dot_shell_d = Component.Input()
    # T_shell_in_d: Temperature of the Molten salt entering the shell side (discharging direction) [K]
    T_shell_in_d = Component.Input()
    # P_shell_in_d: Pressure of the Molten Salt entering the shell side (discharging direction) [Pa]
    P_shell_in_d = Component.Input()
    # m_dot_tube_c: Mass Flow rate of HTF entering the tube side (charging direction) [kg/s]
    m_dot_tube_c = Component.Input()
    # T_tube_in_c: Temperature of the HTF entering the tube side (charging direction) [K]
    T_tube_in_c = Component.Input()
    # P_tube_in_c: Pressure of the HTF entering the tube side (charging direction) [Pa]
    P_tube_in_c = Component.Input()
    # m_dot_tube_d: Mass Flow rate of HTF entering the tube side (discharging direction) [kg/s]
    m_dot_tube_d = Component.Input()
    # T_tube_in_d: Temperature of the HTF entering the tube side (discharging direction) [K]
    T_tube_in_d = Component.Input()
    # P_tube_in_d: Pressure of the HTF entering the tube side (discharging direction) [Pa]
    P_tube_in_d = Component.Input()
    # mass_counter_c: Mass counter for the HTF expansion system (charging direction) [kg]
    mass_counter_c = Component.Input()
    # mass_counter_d: Mass counter for the HTF expansion system (discharging direction) [kg]
    mass_counter_d = Component.Input()

    # *** Model Outputs ***
    # mode_iter: Mode to indicate if the HX is being used in the charging (1) / discharging (0) configuration
    mode_iter = Component.Output()
    # m_dot_shell_c_out: Mass flow rate leaving the shell side of the HX in the charging direction [kg/s]
    m_dot_shell_c_out = Component.Output()
    # T_shell_out_c: Temperature leaving the shell side of the HX in the charging direction [K]
    T_shell_out_c = Component.Output()
    # P_shell_out_c: Pressure leaving the shell side of the HX in the charging direction [Pa]
    P_shell_out_c = Component.Output()
    # m_dot_shell_d_out: Mass flow rate leaving the shell side of the HX in the discharging direction [kg/s]
    m_dot_shell_d_out = Component.Output()
    # T_shell_out_d: Temperature leaving the shell side of the HX in the discharging direction [K]
    T_shell_out_d = Component.Output()
    # P_shell_out_d: Pressure leaving the shell side of the HX in the discharging direction [Pa]
    P_shell_out_d = Component.Output()
    # m_dot_tube_c_out: Mass flow rate leaving the tube side of the HX in the charging direction [kg/s]
    m_dot_tube_c_out = Component.Output()
    # T_tube_out_c: Temperature leaving the tube side of the HX in the charging direction [K]
    T_tube_out_c = Component.Output()
    # P_tube_out_c: Pressure leaving the tube side of the HX in the charging direction [Pa]
    P_tube_out_c = Component.Output()
    # m_dot_tube_d_out: Mass flow rate leaving the tube side of the HX in the discharging direction [kg/s]
    m_dot_tube_d_out = Component.Output()
    # T_tube_out_d: Temperature leaving the tube side of the HX in the discharging direction [K]
    T_tube_out_d = Component.Output()
    # P_tube_out_d: Pressure leaving the tube side of the HX in the discharging direction [Pa]
    P_tube_out_d = Component.Output()
    # mass_counter_c_out: HTF Mass Counter for the expansion system (charging direction) [kg]
    mass_counter_c_out = Component.Output()
    # mass_counter_d_out: HTF Mass Counter for the expansion system (discharging direction) [kg]
    mass_counter_d_out = Component.Output()

    # -----------------------------------------------------------------------------------------------------------------------
    def presim_setup(self, **kwargs):
        """
        Pre-simulation initialisation – analogous to the
        ``getIsFirstCallofSimulation()`` block in TRNSYS.

        Allocates the dynamic storage numpy arrays that hold the nodal
        temperature and enthalpy profiles between time-steps.  These replace
        the TRNSYS ``SetNumberStoredVariables`` / ``getDynamicArrayValueLastTimestep``
        / ``SetDynamicArrayValueThisIteration`` mechanism.
        """
        # SetNumberStoredVariables(0, n_nodes*4 + 1)
        shell_passes = int(self.shell_p.v) * int(self.n_HEX.v)
        n_nodes = int(self.n_baffles.v) * shell_passes + 1

        # Allocate dynamic storage arrays (getDynamicArrayValueLastTimestep /
        # SetDynamicArrayValueThisIteration stored as instance numpy arrays)
        self._T_shell = np.zeros(n_nodes)
        self._T_tube  = np.zeros(n_nodes)
        self._h_shell = np.zeros(n_nodes)
        self._h_tube  = np.zeros(n_nodes)
        # Scalar dynamic variable at index 4*n_nodes + 1: operating mode
        self._mode_prev = 1.0  # assume HX is initially charging to create a temperature profile

    # -----------------------------------------------------------------------------------------------------------------------
    def calculate(self):
        """Calculate heat exchanger outputs for the current simulation iteration."""

        # -----------------------------------------------------------------------------------------------------------------------
        # Calculate the number of nodes
        shell_passes = int(self.shell_p.v) * int(self.n_HEX.v)
        n_nodes = int(self.n_baffles.v) * shell_passes + 1
        d_in  = 2.0 * self.r_in.v
        d_out = 2.0 * self.r_out.v

        # -----------------------------------------------------------------------------------------------------------------------
        # Do All of the First Timestep Manipulations Here - There Are No Iterations at the Initial Time
        if self.model.is_first_step:

            # At first time step, I obtain stationary solution using the e-NTU method
            # I allocate memory for variables: h_shell, h_tube, q_htf_tube, q_htf_shell,
            # T_shell_new, T_tube_new, vel_shell, vel_tube, p_tube
            if n_nodes > 0:
                T_shell_new = np.zeros(n_nodes)
                T_tube_new  = np.zeros(n_nodes)
                h_shell_new = np.zeros(n_nodes)
                h_tube_new  = np.zeros(n_nodes)
            else:
                print("Error: n_nodes debe ser mayor que 0. Valor actual de n_nodes: ", n_nodes)
                return

            mode_iter_v = 1.0  # assume HX is initially charging, to create a temperature profile

            # Set the Initial Values of the Outputs (#,Value)
            self.mode_iter.v         = mode_iter_v              # Mode to indicate if the heat exchanger is being used in the charging/discharging configuration
            self.m_dot_shell_c_out.v = self.m_dot_shell_c.v    # Mass flow rate leaving the shell side of the HX in the charging direction
            self.T_shell_out_c.v     = self.T_shell_in_c.v     # Temperature leaving the shell side of the HX in the charging direction
            self.P_shell_out_c.v     = self.P_shell_in_c.v     # Pressure leaving the shell side of the HX in the charging direction
            self.m_dot_shell_d_out.v = self.m_dot_shell_d.v    # Mass flow rate leaving the shell side of the HX in the discharging direction
            self.T_shell_out_d.v     = self.T_shell_in_d.v     # Temperature leaving the shell side of the HX in the discharging direction
            self.P_shell_out_d.v     = self.P_shell_in_d.v     # Pressure leaving the shell side of the HX in the discharging direction
            self.m_dot_tube_c_out.v  = self.m_dot_tube_c.v     # Mass flow rate leaving the tube side of the HX in the charging direction
            self.T_tube_out_c.v      = self.T_tube_in_c.v      # Temperature leaving the tube side of the HX in the charging direction
            self.P_tube_out_c.v      = self.P_tube_in_c.v      # Pressure leaving the tube side of the HX in the charging direction
            self.m_dot_tube_d_out.v  = self.m_dot_tube_d.v     # Mass flow rate leaving the tube side of the HX in the discharging direction
            self.T_tube_out_d.v      = self.T_tube_in_d.v      # Temperature leaving the tube side of the HX in the discharging direction
            self.P_tube_out_d.v      = self.P_tube_in_d.v      # Pressure leaving the tube side of the HX in the discharging direction

            # e-NTU method for stationary solution
            # set function inputs as charging inputs
            if self.P_shell_in_c.v > 0.0:
                p_in_shell = self.P_shell_in_c.v
            else:  # input not added, set as default
                p_in_shell = 200000.0

            if self.P_tube_in_c.v > 0.0:
                p_in_tube = self.P_tube_in_c.v
            else:
                p_in_tube = 1000000.0

            if self.T_tube_in_c.v > 350.0:
                T_tube_in = self.T_tube_in_c.v
            else:
                T_tube_in = 650.0

            if self.T_shell_in_c.v > 200.0:
                T_shell_in = self.T_shell_in_c.v
            else:
                T_shell_in = 400.0

            # TODO-NEEDS LIBRARY: eNTU subroutine from SergioScripts Fortran module.
            # Computes steady-state shell and tube nodal temperature and enthalpy profiles
            # using the effectiveness-NTU method.
            # Signature (Fortran):
            #   eNTU(p_in_shell, p_in_tube, n_nodes, shell_p, n_tubes, r_out, r_in, L_shell,
            #        D_shell, T_tube_in, T_shell_in, m_tube, m_shell, fluid_shell, config, S_T, S_L,
            #        fluid_tube, n_HEX, n_baffles, T_shell_new, T_tube_new, h_shell, h_tube)
            eNTU(
                p_in_shell, p_in_tube, n_nodes,
                self.shell_p.v, self.n_tubes.v,
                self.r_out.v, self.r_in.v,
                self.L_shell.v, self.D_shell.v,
                T_tube_in, T_shell_in,
                self.m_dot_tube_c.v,    # m_tube (full bundle flow at start)
                self.m_dot_shell_c.v,   # m_shell
                self.fluid_shell.v, self.config.v,
                self.S_T.v, self.S_L.v,
                self.fluid_tube.v, self.n_HEX.v, self.n_baffles.v,
                T_shell_new, T_tube_new, h_shell_new, h_tube_new
            )

            T_shell_out = T_shell_new[n_nodes - 1]
            T_tube_out  = T_tube_new[n_nodes - 1]

            # Store initial profiles into dynamic storage
            # (SetDynamicArrayInitialValue)
            for i in range(n_nodes):
                # SetDynamicArrayInitialValue(i+1, T_shell_new(i+1))
                # - mass of water in the main piping line before the HP turbine
                self._T_shell[i] = T_shell_new[i]
                # SetDynamicArrayInitialValue(n_nodes + i+1, T_tube_new(i+1))
                self._T_tube[i]  = T_tube_new[i]
                # SetDynamicArrayInitialValue(2*n_nodes + i+1, h_shell(i+1))
                self._h_shell[i] = h_shell_new[i]
                # SetDynamicArrayInitialValue(3*n_nodes + i+1, h_tube(i+1))
                self._h_tube[i]  = h_tube_new[i]
            # SetDynamicArrayInitialValue(4*n_nodes + 1, mode_iter)
            self._mode_prev = mode_iter_v

            return

        # -----------------------------------------------------------------------------------------------------------------------
        # ReRead the Parameters if Another Unit of This Type Has Been Called Last
        # (parameter members are always current – no explicit re-read required)

        # -----------------------------------------------------------------------------------------------------------------------
        # Check if current iteration is charging or discharging mode
        mode_prev = self._mode_prev

        if self.m_dot_tube_c.v > self.m_dot_tube_d.v:  # Heat exchanger is charging
            mode_iter_v = 1.0
            if self.m_dot_shell_c.v > 0.0:
                m_shell = self.m_dot_shell_c.v
            else:
                m_shell = 1.0
            if self.T_shell_in_c.v > 0.0:
                T_shell_in = self.T_shell_in_c.v
            else:
                T_shell_in = 600.0  # default charging temperature entering
            if self.P_shell_in_c.v > 0.0:
                p_in_shell = self.P_shell_in_c.v
            else:
                p_in_shell = 101325.0
            if self.P_tube_in_c.v > 0.0:
                p_in_tube = self.P_tube_in_c.v
            else:
                p_in_tube = 101325.0
            m_tubes   = self.m_dot_tube_c.v
            T_tube_in = self.T_tube_in_c.v
            # setDynamicArrayValueThisIteration(4*n_nodes+1, mode_iter)
            # (written to self._mode_prev at end of calculate)

        else:  # Heat exchanger is discharging
            mode_iter_v = 0.0
            if self.m_dot_shell_d.v > 0.0:
                m_shell = self.m_dot_shell_d.v
            else:
                m_shell = 1.0
            if self.T_shell_in_d.v > 0.0:
                T_shell_in = self.T_shell_in_d.v
            else:
                T_shell_in = 600.0  # default charging temperature entering
            if self.P_shell_in_d.v > 0.0:
                p_in_shell = self.P_shell_in_d.v
            else:
                p_in_shell = 101325.0
            if self.P_tube_in_d.v > 0.0:
                p_in_tube = self.P_tube_in_d.v
            else:
                p_in_tube = 101325.0
            m_tubes   = self.m_dot_tube_d.v
            T_tube_in = self.T_tube_in_d.v
            # setDynamicArrayValueThisIteration(4*n_nodes+1, mode_iter)
            # (written to self._mode_prev at end of calculate)

        # Read dynamic storage from last timestep
        # (getDynamicArrayValueLastTimestep)
        if mode_prev != mode_iter_v:
            # Need to adjust temperature/enthalpy profiles on shell and tube
            # side to account for change in direction:
            #   T_shell(i) = getDynamicArrayValueLastTimestep(n_nodes - i + 1)
            # which is equivalent to reversing the stored arrays.
            T_shell = self._T_shell[::-1].copy()
            T_tube  = self._T_tube[::-1].copy()
            h_shell = self._h_shell[::-1].copy()
            h_tube  = self._h_tube[::-1].copy()
        else:
            T_shell = self._T_shell.copy()
            T_tube  = self._T_tube.copy()
            h_shell = self._h_shell.copy()
            h_tube  = self._h_tube.copy()

        # -----------------------------------------------------------------------------------------------------------------------
        gap_baffle = self.L_shell.v / (self.n_baffles.v + 1)           # Delta_x
        m_tube     = m_tubes / self.n_tubes.v                            # Mass flow at each tube of the bundle
        A_baffle   = (2.0 / 5.0) * 3.1415 * (self.D_shell.v ** 2) / 4.0  # Proportional to the transversal area of the HEX

        q_HTF_tube  = np.zeros(n_nodes)
        q_HTF_shell = np.zeros(n_nodes)
        T_shell_new = np.zeros(n_nodes)
        T_tube_new  = np.zeros(n_nodes)
        vel_shell   = np.zeros(n_nodes)
        vel_tube    = np.zeros(n_nodes)

        # Calculate the heat transfer between the shell and tube side
        # TODO-NEEDS LIBRARY: ZeroD_Eq subroutine from SergioScripts Fortran module.
        # Computes nodal heat transfer terms q_HTF_tube and q_HTF_shell from
        # current enthalpy and temperature profiles.
        # Signature (Fortran):
        #   ZeroD_Eq(n_nodes, h_tube, h_shell, T_shell, T_tube, k_steel,
        #            gap_baffle, d_in, d_out, q_HTF_tube, q_HTF_shell)
        ZeroD_Eq(
            n_nodes, h_tube, h_shell, T_shell, T_tube,
            self.k_steel.v, gap_baffle, d_in, d_out,
            q_HTF_tube, q_HTF_shell
        )

        # TODO-NEEDS LIBRARY: shell subroutine from SergioScripts Fortran module.
        # Updates shell-side nodal temperatures, enthalpy profile, velocities,
        # and outlet pressure.  h_shell is updated in-place.
        # Signature (Fortran):
        #   shell(config, shell_p, Tenv, p_in_shell, gap_baffle, n_nodes, k_steel, k_ins,
        #         S_T, S_L, m_shell, n_baffles, n_HEX, th_baffle, th_ins, T_shell, D_shell,
        #         L_baffle, A_baffle, shell_passes, L_shell, n_tubes, D_out, q_HTF_shell,
        #         v_out, timestep, fluid_shell, T_shell_new, h_shell, vel_shell, p_shell_out)
        p_shell_out = 0.0
        shell(
            self.config.v, self.shell_p.v, self.Tenv.v, p_in_shell,
            gap_baffle, n_nodes, self.k_steel.v, self.k_ins.v,
            self.S_T.v, self.S_L.v, m_shell, self.n_baffles.v,
            self.n_HEX.v, self.th_baffle.v, self.th_ins.v,
            T_shell, self.D_shell.v, self.L_baffle.v, A_baffle,
            shell_passes, self.L_shell.v, self.n_tubes.v, d_out,
            q_HTF_shell, self.v_out.v, self.model.settings.timestep,
            self.fluid_shell.v,
            T_shell_new, h_shell, vel_shell, p_shell_out
        )

        # TODO-NEEDS LIBRARY: tube_oneD_inc subroutine from SergioScripts Fortran module.
        # Updates tube-side nodal temperatures, enthalpy profile, velocities,
        # and outlet pressure.  h_tube is updated in-place.
        # Signature (Fortran):
        #   tube_oneD_inc(q_HTF_tube, T_tube, n_nodes, gap_baffle, d_in,
        #                 p_in_tube, m_tube, n_tubes, fluid_tube, timestep,
        #                 T_tube_new, h_tube, vel_tube, p_tube_out)
        p_tube_out = 0.0
        tube_oneD_inc(
            q_HTF_tube, T_tube, n_nodes, gap_baffle, d_in,
            p_in_tube, m_tube, self.n_tubes.v, self.fluid_tube.v,
            self.model.settings.timestep,
            T_tube_new, h_tube, vel_tube, p_tube_out
        )

        T_shell_out = T_shell_new[n_nodes - 1]
        T_tube_out  = T_tube_new[n_nodes - 1]

        # Save new Dynamic Storage Variables
        # (setDynamicArrayValueThisIteration)
        for i in range(n_nodes):
            # setDynamicArrayValueThisIteration(i+1, T_shell_new(i+1))
            self._T_shell[i] = T_shell_new[i]
            # setDynamicArrayValueThisIteration(n_nodes + i+1, T_tube_new(i+1))
            self._T_tube[i]  = T_tube_new[i]
            # setDynamicArrayValueThisIteration(2*n_nodes + i+1, h_shell(i+1))
            # h_shell has been updated in-place by shell()
            self._h_shell[i] = h_shell[i]
            # setDynamicArrayValueThisIteration(3*n_nodes + i+1, h_tube(i+1))
            # h_tube has been updated in-place by tube_oneD_inc()
            self._h_tube[i]  = h_tube[i]
        # setDynamicArrayValueThisIteration(4*n_nodes+1, mode_iter)
        self._mode_prev = mode_iter_v

        # -----------------------------------------------------------------------------------------------------------------------
        # Set Outputs
        # pressure could be negative while hydraulic solver is iterating over many types
        # in either mode and for both the shell and the tube pressures
        # if it is negative pass it along so that the hydraulic solver can figure itself out
        # otherwise output the calculated pressure drop
        if mode_iter_v == 1.0:
            # Charging Outputs
            T_shell_out_c_val = T_shell_out
            if self.P_shell_in_c.v < 0.0:
                P_shell_out_c_val = self.P_shell_in_c.v
            else:
                P_shell_out_c_val = p_shell_out
            T_tube_out_c_val = T_tube_out
            if self.P_tube_in_c.v < 0.0:
                P_tube_out_c_val = self.P_tube_in_c.v
            else:
                P_tube_out_c_val = p_tube_out

            # Discharging Outputs
            T_shell_out_d_val = self.T_shell_in_d.v
            P_shell_out_d_val = self.P_shell_in_d.v
            T_tube_out_d_val  = self.T_tube_in_d.v
            P_tube_out_d_val  = self.P_tube_in_d.v

        else:
            # Discharging Outputs
            T_shell_out_d_val = T_shell_out
            if self.P_shell_in_d.v < 0.0:
                P_shell_out_d_val = self.P_shell_in_d.v
            else:
                P_shell_out_d_val = p_shell_out
            T_tube_out_d_val = T_tube_out
            if self.P_tube_in_d.v < 0.0:
                P_tube_out_d_val = self.P_tube_in_d.v
            else:
                P_tube_out_d_val = p_tube_out

            # Charging Outputs
            T_shell_out_c_val = self.T_shell_in_c.v
            P_shell_out_c_val = self.P_shell_in_c.v
            T_tube_out_c_val  = self.T_tube_in_c.v
            P_tube_out_c_val  = self.P_tube_in_c.v

        # -----------------------------------------------------------------------------------------------------------------------
        # Set the Outputs from this Model (#,Value)
        # Call setOutputValue(1, 0.d0)  !Mode – commented out in original Fortran; mode_iter still tracked
        self.m_dot_shell_c_out.v = self.m_dot_shell_c.v    # Mass flow rate leaving the shell side of the HX in the charging direction
        self.T_shell_out_c.v     = T_shell_out_c_val        # Temperature leaving the shell side of the HX in the charging direction
        self.P_shell_out_c.v     = P_shell_out_c_val        # Pressure leaving the shell side of the HX in the charging direction
        self.m_dot_shell_d_out.v = self.m_dot_shell_d.v    # Mass flow rate leaving the shell side of the HX in the discharging direction
        self.T_shell_out_d.v     = T_shell_out_d_val        # Temperature leaving the shell side of the HX in the discharging direction
        self.P_shell_out_d.v     = P_shell_out_d_val        # Pressure leaving the shell side of the HX in the discharging direction
        self.m_dot_tube_c_out.v  = self.m_dot_tube_c.v     # Mass flow rate leaving the tube side of the HX in the charging direction
        self.T_tube_out_c.v      = T_tube_out_c_val         # Temperature leaving the tube side of the HX in the charging direction
        self.P_tube_out_c.v      = P_tube_out_c_val         # Pressure leaving the tube side of the HX in the charging direction
        self.m_dot_tube_d_out.v  = self.m_dot_tube_d.v     # Mass flow rate leaving the tube side of the HX in the discharging direction
        self.T_tube_out_d.v      = T_tube_out_d_val         # Temperature leaving the tube side of the HX in the discharging direction
        self.P_tube_out_d.v      = P_tube_out_d_val         # Pressure leaving the tube side of the HX in the discharging direction
        self.mass_counter_c_out.v = self.mass_counter_c.v   # HTF Mass Counter used for the expansion system in the charging direction
        self.mass_counter_d_out.v = self.mass_counter_d.v   # HTF Mass Counter used for the expansion system in the discharging direction
