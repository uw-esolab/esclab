"""Thermal Energy Storage (TES) Tank component model (Type 4100)."""

# Object: ThermalStTank
# Simulation Studio Model: ESOL4100-TESTank

# Author: Sergio Alcalde-Morales
# Editor:
# Date:    June 07, 2024
# last modified: June 07, 2024
# Ported by: GitHub Copilot, March 01, 2026

import math

import numpy as np

from esclab.simulate import Component
from esclab.components.esol_properties import Incompressible

Inc = Incompressible()


class TESTank(Component):
    """
    TRNSYS Type 4100: ESOL4100-TESTank.

    Simulates a well-mixed thermal energy storage tank with heat loss to the
    environment through insulation, radiation, and convection (natural and
    forced).  The energy balance is stepped forward in sub-timesteps whenever
    the TRNSYS timestep exceeds the critical sub-step size.

    Parameters
    ----------
    D_in : float
        Inner diameter of the tank [m].
    Height : float
        Height of the tank [m].
    Ins_Th : float
        Insulation thickness around the tank [m].
    k_iso : float
        Thermal conductivity of the insulation [W/(m·K)].
    Emiss : float
        Emissivity between the tank wall and the surroundings [-].
    T0 : float
        Initial tank temperature at the first timestep [K].
    L0 : float
        Initial tank fluid level at the first timestep [m].
    ID_Fluid : str
        Fluid identifier for HTF property lookups via
        ``esclab.components.esol_properties.Incompressible``.

    Inputs
    ------
    T_in : float
        Temperature of fluid entering the tank [K].
    m_in : float
        Mass flow rate of fluid entering the tank [kg/s].
    m_out : float
        Mass flow rate of fluid drawn from the tank [kg/s].
    T_env : float
        Temperature of the ambient air outside the tank [K].
    v_env : float
        Velocity of the air around the sides of the tank [m/s].
    v_air : float
        Velocity of the forced air underneath the tank [m/s].

    Outputs
    -------
    m_dot_out : float
        Mass flow rate [kg/s] of fluid leaving the tank (drawn by the pump).
    Vol_dot_out : float
        Volumetric flow rate [m³/s] of fluid leaving the tank.
    T_out : float
        Temperature [K] of the fluid leaving the tank.
    P_salt : float
        Pressure [Pa] of the fluid in the tank.
    T_tank : float
        Temperature of the tank [K] at this timestep.
    L_tank : float
        Level of the tank [m] at this timestep.
    m_tank : float
        Mass of fluid in the tank [kg] at this timestep.
    Tw : float
        Temperature of the tank wall [K].
    Q_loss : float
        Heat lost to the environment over the timestep [J].
    """

    # *** Model Parameters ***
    D_in = Component.Parameter()      # Inner Diameter of the tank [m]
    Height = Component.Parameter()    # Height of the tank [m]
    Ins_Th = Component.Parameter()    # Insulation thickness around the tank [m]
    k_iso = Component.Parameter()     # Thermal Conductivity of insulation
    Emiss = Component.Parameter()     # Emissivity between the tank and the surroundings
    T0 = Component.Parameter()        # Initial Tank Temperature at the first timestep [K]
    L0 = Component.Parameter()        # Initial Tank Level at the first timestep [m]
    ID_Fluid = Component.Parameter()  # Fluid ID to find fluid properties (40 = Dowtherm A)

    # *** Model Inputs ***
    T_in = Component.Input()   # Temperature of salt entering the tank [K]
    m_in = Component.Input()   # Mass of the salt entering the tank [kg/s]
    m_out = Component.Input()  # Mass of salt drawn from the tank [kg/s]
    T_env = Component.Input()  # Temperature of ambient air outside of the tank [K]
    v_env = Component.Input()  # Velocity of the air around the sides of the tank [m/s]
    v_air = Component.Input()  # Velocity of the air forced underneath the tank [m/s]

    # *** Model Outputs ***
    m_dot_out = Component.Output()   # Mass flow rate [kg/s] of molten salt leaving the tank (drawn by the pump)
    Vol_dot_out = Component.Output() # Volumetric flow rate [m^3/s] of molten salt leaving the tank (drawn by the pump)
    T_out = Component.Output()       # Temperature [K] of the molten salt leaving the tank (drawn by the pump)
    P_salt = Component.Output()      # Pressure [Pa] of the molten salt in the tank
    T_tank = Component.Output()      # Temperature of the tank [K] at this timestep
    L_tank = Component.Output()      # Level of the tank [m] at this timestep
    m_tank = Component.Output()      # Mass of the tank [kg] at this timestep
    Tw = Component.Output()          # Temperature of the tank wall [K]
    Q_loss = Component.Output()      # Heat lost to the environment [J]

    def presim_setup(self, **kwargs):
        super().presim_setup(**kwargs)

        # Dynamic storage array (2 elements):
        #   index 0 -> T_tank  (SetDynamicArrayInitialValue(1, T0))
        #   index 1 -> L_tank  (SetDynamicArrayInitialValue(2, L0))
        self._dynamic = np.zeros(2)
        self._dynamic[0] = self.T0.v  # Initial Temperature of the Tank
        self._dynamic[1] = self.L0.v  # Initial Level of the Tank

    def calculate(self):
        super().calculate()

        # Perform Any "After Convergence" Manipulations That May Be Required at the End of Each Timestep
        if self.model.is_converged:
            return

        # Do All of the First Timestep Manipulations Here - There Are No Iterations at the Initial Time
        if self.model.is_first_step:
            rho_salt = Inc.density(self.ID_Fluid.v, T=self.T0.v)  # Density evaluated at T_tank
            m_tank = rho_salt * (self.L0.v * 3.1415 * self.D_in.v**2.0) / 4.0  # Mass in the tank at the start of the sim
            Tw = 305.0  # initial wall temperature guess

            P_out = 101325.0 + 9.81 * rho_salt * self.L0.v
            T_out = self.T0.v

            # Set Dynamic Storage Values
            self._dynamic[0] = self.T0.v  # Initial Temperature of the Tank
            self._dynamic[1] = self.L0.v  # Initial Level of the Tank

            # Set the Initial Values of the Outputs (#,Value)
            self.m_dot_out.v = 0.0        # Mass flow rate [kg/s] of molten salt leaving the tank (drawn by the pump)
            self.Vol_dot_out.v = 0.0      # Volumetric flow rate [m^3/s] of molten salt leaving the tank (drawn by the pump)
            self.T_out.v = T_out          # Temperature [K] of the molten salt leaving the tank (drawn by the pump)
            self.P_salt.v = P_out         # Pressure [Pa] of the molten salt leaving the tank (drawn by the pump)
            self.T_tank.v = self.T0.v     # Temperature of the tank [K] at this timestep
            self.L_tank.v = self.L0.v     # Level of the tank [m] at this timestep
            self.m_tank.v = m_tank        # Mass of the tank [kg] at this timestep
            self.Tw.v = Tw                # Temperature of the tank wall [K]
            self.Q_loss.v = 0.0           # Heat lost to the environment [J]
            return

        # Read the Inputs
        # (inputs accessed directly via self.<input>.v throughout the calculation)

        # Constants in the code
        g = 9.81             # Gravity
        mu_air = 184.6e-07   # Viscosity of air (considered constant)
        k_air = 26.3e-03     # Thermal Conductivity of air (considered constant)
        rho_air = 1.1614     # Density of Air (considered constant)
        cp_air = 1.005       # Specific heat of air (considered constant)
        visc_air = 15.89e-06 # kinematic viscosity of air [m^2/s]
        alpha_air = 22.5e-06 # thermal diffusivity of air
        beta_air = 1.0 / self.T_env.v  # volumetric expansion coefficient
        P_salt = 101325.0    # Pressure of the salt in the tank [Pa]
        t_crit = 10.0        # Critical Timestep when step forward is split into multiple steps per iteration

        ##################### BREAK UP TIMESTEP INTO SUB-TIMESTEPS IF NEEDED #####################
        ts = self.model.timestep * 3600.0  # time step in seconds (is given in hours)
        # finaltime = int(time)  # Final time defined as an integer
        if ts > t_crit:
            sub_timestep = 10.0  # At each time step, I divide it into 10 parts to obtain a better solution
            N = math.ceil(ts / sub_timestep)  # Number of timesteps
            sub_timestep = ts / float(N)
        else:
            sub_timestep = ts
            N = 1

        #################### TANK LEVEL AND TEMP AT THE END OF LAST TIMESTEP ####################
        # GetDynamicArrayValueLastTimestep(1) -> self._dynamic[0]
        # GetDynamicArrayValueLastTimestep(2) -> self._dynamic[1]
        T_tank = self._dynamic[0]
        L_tank = self._dynamic[1]

        #################### SALT Properties based on Fluid Temperature ####################
        rho_salt = Inc.density(self.ID_Fluid.v, T=T_tank)
        m_tank = rho_salt * (L_tank * 3.1415 * self.D_in.v**2.0) / 4.0
        # CONVERTED-NEEDS UNITS CHECK: removed *1000; esol_properties.specheat already returns J/(kg·K)
        cp_salt_in = Inc.specheat(self.ID_Fluid.v, T=self.T_in.v)

        ################# TANK PROPERTIES #################
        A_total = 3.1415 * self.D_in.v * self.Height.v + 3.1415 * self.D_in.v * self.D_in.v / 4.0  # Total area of the tank (walls and bottom)
        A_ext = 3.1415 * (self.D_in.v + 2.0 * self.Ins_Th.v) * self.Height.v

        Q_loss = 0.0
        # Read current wall temperature output (getOutputValue(8) -> self.Tw.v)
        Tw = self.Tw.v

        # Here starts a "for" loop to simulate each "subtime step"
        T_new = T_tank  # ensure T_new is defined before the loop
        for i in range(1, N + 1):
            # Properties of the molten salt at previous time step
            rho_salt = Inc.density(self.ID_Fluid.v, T=T_tank)
            # CONVERTED-NEEDS UNITS CHECK: removed *1000; esol_properties.specheat already returns J/(kg·K)
            cp_salt = Inc.specheat(self.ID_Fluid.v, T=T_tank)  # Specific heat evaluated at T_tank (initial temperature at each subtime step)

            # Here starts the "while" loop
            error = 1.0
            while error > 0.000001:

                # Characteristics of the isolation
                T_iso = (Tw + T_tank) / 2.0  # Characteristic temperature of the isolation
                # k_iso = 0.04
                # Characteristics of the film air outside the tank
                T_film = (self.T_env.v + Tw) / 2.0  # Not used since T_film doesn't change enough to produce significantly changes in the air properties

                # Thermal Resistances
                # Conductivity resistance
                R_I = self.Ins_Th.v / (A_total * self.k_iso.v)
                # Radiative resistance
                sigma = 5.67e-08
                h_rad = self.Emiss.v * sigma * (self.T_env.v**2.0 + Tw**2.0) * (self.T_env.v + Tw)
                R_II = 1.0 / (A_ext * h_rad)
                # Convective resistance (natural convection)
                Pr_ext = 0.7
                if self.v_env.v < 0.1:
                    Ra_ext = (beta_air * (Tw - self.T_env.v) * g * self.Height.v**3) / visc_air / alpha_air
                    Nu_ext = 0.68 + (0.67 * Ra_ext**0.25) / (1.0 + (0.492 / Pr_ext)**(9.0 / 16.0))
                else:
                    Re_ext = rho_air * self.v_env.v * self.Height.v / mu_air
                    if Re_ext < 3000.0:
                        Nu_ext = 0.664 * Re_ext**0.5 * Pr_ext**0.333
                    else:
                        Nu_ext = 0.0296 * Re_ext**0.8 * Pr_ext**0.333

                h_ext = k_air * Nu_ext / self.Height.v
                R_III = 1.0 / (A_ext * h_ext)

                # Convective resistance (forced convection)
                D_ext = self.D_in.v + 2.0 * self.Ins_Th.v
                Re_air = rho_air * self.v_air.v * D_ext / mu_air
                Nu_air = 0.0296 * Re_air**0.8 * Pr_ext**0.333
                h_air = k_air * Nu_air / D_ext
                R_IV = 1.0 / (D_ext * h_air)

                # Parallel external resistances
                R_V = (R_II * R_III * R_IV) / (R_III * R_IV + R_II * R_II + R_II * R_IV)

                m = m_tank + (self.m_in.v - self.m_out.v) * sub_timestep  # Mass variation inside the tank
                h_in = cp_salt_in * self.T_in.v   # Inlet enthalpy
                h_out = cp_salt * T_tank           # Outlet enthalpy

                Div = m * cp_salt / sub_timestep + cp_salt * (self.m_in.v - self.m_out.v) + 1.0 / (R_I + R_V) + (cp_salt * self.m_out.v)
                T_new = ((m * cp_salt * T_tank / sub_timestep) + (h_in * self.m_in.v) + (self.T_env.v / (R_I + R_V))) / Div

                HL = (T_new - self.T_env.v) / (R_I + R_V)  # Heat losses

                Tww = Tw           # Set Previous wall temperature
                Tw = HL * R_V + self.T_env.v  # New wall temperature
                error = (Tw - Tww) * (Tw - Tww)

            # **** FIN BUCLE WHILE PARA OBTENER TW CON TNEW

            # Once "error" is lower than 1E-10, we go into the next step time
            T_tank = T_new  # Temperature at each subtime step
            m_tank = m      # Mass at each subtime step
            Q_loss = Q_loss + HL * sub_timestep

        # We go to the next subtime step

        # **** FIN BUCLE FOR DE SUBPASOS DE TIEMPO: Ya tenemos T para el siguiente paso

        # ** It is important to distinguish between time step (given by trnsys)...
        # and subtime_steps, used to obtain the solution at each time step
        T_tank = T_new  # Dynamic evolution of temperature of the tank --> Output which is an input in the next time step
        L_tank = (m_tank / rho_salt) / (3.1415 * self.D_in.v**2.0 / 4.0)  # Dynamic evolution of mass of the tank --> Output which is an input in the next time step

        # SetDynamicArrayValueThisIteration(1, T_tank) -> self._dynamic[0] = T_tank
        # SetDynamicArrayValueThisIteration(2, L_tank) -> self._dynamic[1] = L_tank
        self._dynamic[0] = T_tank
        self._dynamic[1] = L_tank

        T_out = T_tank
        Vol_dot_out = self.m_out.v / rho_salt

        self.m_dot_out.v = self.m_out.v    # Mass flow rate [kg/s] of molten salt leaving the tank (drawn by the pump)
        self.Vol_dot_out.v = Vol_dot_out   # Volumetric flow rate [m^3/s] of molten salt leaving the tank (drawn by the pump)
        self.T_out.v = T_out               # Temperature [K] of the molten salt leaving the tank (drawn by the pump)
        self.P_salt.v = P_salt             # Pressure [Pa] of the molten salt in the tank
        self.T_tank.v = T_tank             # Temperature of the tank [K] at this timestep
        self.L_tank.v = L_tank             # Level of the tank [m] at this timestep
        self.m_tank.v = m_tank             # Mass of the tank [kg] at this timestep
        self.Tw.v = Tw                     # Temperature of the tank wall [K]
        self.Q_loss.v = Q_loss             # Heat lost to the environment [J]
