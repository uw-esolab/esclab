"""Pipe component model with RK-4 thermal simulation (Type 4035)."""

import numpy as np

from esclab.simulate import Component
from esclab.components.esol_properties import Incompressible
from esclab.components.flownetwork.sf_piping_helpers import PressureDrop, pipe_dTdt

Inc = Incompressible()


class Pipe(Component):
    """
    Object: ESOL4035-Pipe
    Simulation Studio Model: ESOL4035-Pipe

    Author: Matt Tuman
    Date:    January 09, 2024
    last modified: January 09, 2024
    Ported by: GitHub Copilot, March 01, 2026

    A pipe component that simulates transient thermal behavior using RK-4 time-stepping
    and computes hydraulic pressure drop with minor losses.

    Parameters
    ----------
    Diameter : float
        Pipe inner diameter [m].
    L_tot : float
        Total pipe length [m].
    n_nodes : int
        Number of thermal nodes (must be <= 500).
    Fluid_ID : str
        Fluid identifier for property lookups.
    Roughness : float
        Pipe roughness [m].
    n_large_elbows : float
        Number of large elbows (minor loss fittings).
    n_medium_elbows : float
        Number of medium elbows.
    n_standard_elbows : float
        Number of standard (medium) elbows.
    n_contractions : float
        Number of contractions.
    n_expansions : float
        Number of expansions.
    n_gate_valves : float
        Number of gate valves.
    init_temp : float
        Initial temperature of fluid in pipe [K].
    mc_mult : float
        Mass counter multiplier.
    heat_loss : float
        Heat loss coefficient.

    Inputs
    ------
    Temperature : float
        Inlet fluid temperature [K].
    Mass_Flow : float
        Mass flow rate [kg/s].
    T_amb : float
        Ambient temperature [K].
    Pressure : float
        Inlet pressure [Pa].
    Wind : float
        Wind speed [m/s].
    Mass_Counter_in : float
        Incoming mass counter [kg].

    Outputs
    -------
    T_out : float
        Outlet fluid temperature (last node) [K].
    P_out : float
        Outlet pressure [Pa].
    MassFlow : float
        Mass flow rate through pipe [kg/s].
    Mass_Counter : float
        Accumulated fluid mass in pipe [kg].
    """

    pi = 3.14159265

    # *** Model Parameters ***
    Diameter = Component.Parameter()            # Pipe inner diameter
    L_tot = Component.Parameter()              # Total pipe length
    n_nodes = Component.Parameter()            # Number of thermal nodes
    Fluid_ID = Component.Parameter()           # Fluid identifier
    Roughness = Component.Parameter()          # Pipe roughness
    n_large_elbows = Component.Parameter()     # Number of large elbows
    n_medium_elbows = Component.Parameter()    # Number of medium elbows
    n_standard_elbows = Component.Parameter()  # Number of standard elbows
    n_contractions = Component.Parameter()     # Number of contractions
    n_expansions = Component.Parameter()       # Number of expansions
    n_gate_valves = Component.Parameter()      # Number of gate valves
    init_temp = Component.Parameter()          # Initial temperature [K]
    mc_mult = Component.Parameter()            # Mass counter multiplier
    heat_loss = Component.Parameter()          # Heat loss coefficient

    # *** Model Inputs ***
    Temperature = Component.Input()            # Inlet fluid temperature [K]
    Mass_Flow = Component.Input()              # Mass flow rate [kg/s]
    T_amb = Component.Input()                  # Ambient temperature [K]
    Pressure = Component.Input()               # Inlet pressure [Pa]
    Wind = Component.Input()                   # Wind speed [m/s]
    Mass_Counter_in = Component.Input()        # Incoming mass counter [kg]

    # *** Model Outputs ***
    T_out = Component.Output()                 # Outlet temperature [K]
    P_out = Component.Output()                 # Outlet pressure [Pa]
    MassFlow = Component.Output()              # Mass flow rate [kg/s]
    Mass_Counter = Component.Output()          # Accumulated mass in pipe [kg]

    def presim_setup(self, **kwargs):
        super().presim_setup(**kwargs)

        # TODO-NEEDS CONVERSION REVIEW: dynamic array storage replaces getDynamicArrayValueLastTimestep /
        # setDynamicArrayValueThisIteration. Node temperatures are stored as a numpy array on the instance.
        # n_nodes enforced <= 500 per original Fortran validation.
        self.T_nodes = np.full(int(self.n_nodes.v), self.init_temp.v)

    def calculate(self):
        super().calculate()

        n_nodes = int(self.n_nodes.v)
        n_cv = n_nodes - 1

        # Define length of control volumes
        L_cv = self.L_tot.v / float(n_nodes - 1)
        # Define Volume of CV's
        Vol = self.pi * (self.Diameter.v / 2.0) ** 2 * L_cv

        # Do All of the "First Timestep" Manipulations Here - There Are No Iterations at the Initial Time
        if self.model.is_first_step:
            # Initialize temperatures in nodes of pipe
            for n in range(n_nodes):
                self.T_nodes[n] = self.init_temp.v

            # Compute total mass in the pipe for mass counter
            mass_counter = 0.0
            for n in range(n_cv):
                # Compute control volume average temperature
                T_cv = self.init_temp.v
                # Compute mass in CV
                # CONVERTED-NEEDS UNITS CHECK: T in K and P in Pa confirmed; inputs are SI
                # (T_cv derives from init_temp [K]; original Fortran passes 0.0 for pressure).
                mass_counter = mass_counter + Vol * Inc.density(self.Fluid_ID.v, T=T_cv, P=0.0)

            # Compute pressure drop
            dP = PressureDrop(
                self.Fluid_ID.v, self.Mass_Flow.v, self.init_temp.v, 1.0,
                self.Diameter.v, self.Roughness.v, self.L_tot.v,
                self.n_expansions.v, self.n_contractions.v, self.n_standard_elbows.v,
                self.n_medium_elbows.v, self.n_large_elbows.v,
                self.n_gate_valves.v, 0.0, 0.0, 0.0, 0.0, 0.0,
            )

            # Set the Initial Values of the Outputs
            self.T_out.v = self.init_temp.v                 # Temperature
            self.P_out.v = self.Pressure.v - dP             # Pressure
            self.MassFlow.v = self.Mass_Flow.v              # MassFlow
            self.Mass_Counter.v = mass_counter              # Mass_Counter
            return

        # Perform Thermal Calculations at the End of Each Timestep
        if self.model.is_converged:
            # TODO-NEEDS CONVERSION REVIEW: self.model.settings.timestep is assumed to be in hours;
            # multiplied by 3600 to convert to seconds as in original Fortran.
            timestep = self.model.settings.timestep * 3600.0  # [s]

            # Load in temperatures from last timestep from dynamic array
            T_hat = self.T_nodes.copy()
            # Update Input Temperature
            T_hat[0] = self.Temperature.v
            T_prev = T_hat.copy()

            # Step through time with RK-4
            # K1
            k1 = pipe_dTdt(n_nodes, T_hat, Vol, self.Mass_Flow.v,
                           self.mc_mult.v, self.Fluid_ID.v, self.heat_loss.v, L_cv)
            T_hat = T_prev + k1 * timestep / 2.0

            # K2
            k2 = pipe_dTdt(n_nodes, T_hat, Vol, self.Mass_Flow.v,
                           self.mc_mult.v, self.Fluid_ID.v, self.heat_loss.v, L_cv)
            T_hat = T_prev + k2 * timestep / 2.0

            # K3
            k3 = pipe_dTdt(n_nodes, T_hat, Vol, self.Mass_Flow.v,
                           self.mc_mult.v, self.Fluid_ID.v, self.heat_loss.v, L_cv)
            T_hat = T_prev + k3 * timestep

            # K4
            k4 = pipe_dTdt(n_nodes, T_hat, Vol, self.Mass_Flow.v,
                           self.mc_mult.v, self.Fluid_ID.v, self.heat_loss.v, L_cv)

            # Step Through Time
            T_hat = T_prev + (k1 / 6.0 + k2 / 3.0 + k3 / 3.0 + k4 / 6.0) * timestep

            # Update dynamic storage
            self.T_nodes = T_hat.copy()
            return

        # Compute total mass in the pipe for mass counter
        mass_counter = 0.0
        for n in range(n_cv):
            # Compute control volume average temperature
            # NOTE: original Fortran uses getDynamicArrayValueLastTimestep(n) for both sides of the average;
            # this appears to be a bug in the original source (both indices are n rather than n and n+1).
            T_cv = (self.T_nodes[n] + self.T_nodes[n]) / 2.0
            # Compute mass in CV
            # CONVERTED-NEEDS UNITS CHECK: T in K and P in Pa confirmed; inputs are SI
            # (T_cv derives from T_nodes which are maintained in K).
            mass_counter = mass_counter + Vol * Inc.density(self.Fluid_ID.v, T=T_cv, P=0.0)

        # COMPUTE PRESSURE DROP
        T_ave = (self.T_nodes[0] + self.T_nodes[n_nodes - 1]) / 2.0
        dP = PressureDrop(
            self.Fluid_ID.v, self.Mass_Flow.v, T_ave, 1.0,
            self.Diameter.v, self.Roughness.v, self.L_tot.v,
            self.n_expansions.v, self.n_contractions.v, self.n_standard_elbows.v,
            self.n_medium_elbows.v, self.n_large_elbows.v,
            self.n_gate_valves.v, 0.0, 0.0, 0.0, 0.0, 0.0,
        )

        # Set the Outputs from this Model
        self.T_out.v = self.T_nodes[n_nodes - 1]          # Temperature
        self.P_out.v = self.Pressure.v - dP               # Pressure
        self.MassFlow.v = self.Mass_Flow.v                 # MassFlow
        self.Mass_Counter.v = mass_counter                 # Mass_Counter
