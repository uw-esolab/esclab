"""Variable speed pump component model (Type 4004)."""

import numpy as np

from esclab.simulate import Component
from esclab.components.esol_properties import Incompressible as Inc


class VarSpeedPump(Component):
    """
    Docstring for VarSpeedPump
    Subroutine Type4004
    Object: ESOL4004-VarPump

    Simulation Studio Model: ESOL4004-VarPump

    Author: Matt Tuman
    Editor: Mike Wagner
    Date:     January 05, 2023
    last modified: January 05, 2023
    Converted by: Mike Wagner, February 12, 2026
    """

    #     PARAMETERS
    N_pumps_parallel = Component.Parameter()         #  Number of pumps that are in parallel
    Pump_a0 = Component.Parameter()                  #  head = Pump_a0 + Pump_a1*Q_dot + Pump_a2*Q_dot^2 [-]
    Pump_a1 = Component.Parameter()                  #  head = Pump_a0 + Pump_a1*Q_dot + Pump_a2*Q_dot^2 [1/m^3/s]
    Pump_a2 = Component.Parameter()                  #  head = Pump_a0 + Pump_a1*Q_dot + Pump_a2*Q_dot^2 [1/m^3/s]
    NPSH_a0 = Component.Parameter()                  #  NPSHr = Pump_a0 + Pump_a1*Q_dot + Pump_a2*Q_dot^2 [-]
    NPSH_a1 = Component.Parameter()                  #  NPSHr = Pump_a0 + Pump_a1*Q_dot + Pump_a2*Q_dot^2 [1/m^3/s]
    NPSH_a2 = Component.Parameter()                  #  NPSHr = Pump_a0 + Pump_a1*Q_dot + Pump_a2*Q_dot^2 [1/m^3/s]
    RPM_full = Component.Parameter()                 #  RPM of pump at 100% speed
    LR = Component.Parameter()                       #  Learning rate: Close to 1 means that the mass flow can change significantly [0-1]
    D_outlet = Component.Parameter()                 #  Diameter of the pump outlet [-]
    Fluid_ID = Component.Parameter()                 #  Fluid ID
    Pump_Solver = Component.Parameter()              #  Determines if hydraulic solver is used to update mass flow

    #     INPUTS
    Mass_Flow = Component.Input()                #  Mass flow into the pump [kg/s]
    Pressure = Component.Input()                 #  Pressure of fluid at inlet of the pump [Pa]
    Temperature = Component.Input()              #  Temperature of fluid at inlet of pump [C]
    speed = Component.Input()                    #  Speed that the pump is operating at [0-1]
    error = Component.Input()                    #  Error accumulated within pressure drops throughout the system [Pa]
    mass_count = Component.Input()               #  Total mass counted before the pump [kg]

    #      Outputs
    m_dot_out = Component.Output() #  Mass Flow
    P_out = Component.Output() #  Pressure
    Temperature_out = Component.Output() #  Temperature
    mass_count_out = Component.Output()
    cavitation = Component.Output()
    g = 9.81 #  gravity

    # -----------------------------------------------------------------------------------------------------------------------
    def presim_setup(self, **kwargs):
        # Set the Initial Values of the Outputs (#,Value)
        self.m_dot_out.v = self.Mass_Flow.v
        self.P_out.v = self.Pressure.v
        self.Temperature_out.v = self.Temperature.v
        self.mass_count_out.v = self.mass_count.v
        self.cavitation.v = 0    #  Not cavitating initially

    # -----------------------------------------------------------------------------------------------------------------------
    def calculate(self):
        # -----------------------------------------------------------------------------------------------------------------------
        # Post-convergence
        if self.model.is_converged:
            # Once model has converged, check if pump is likely cavitating
            rho = Inc.density(self.Fluid_ID.v, self.Temperature.v, 0.0)
            Q_dot = self.Mass_Flow.v/rho/self.N_pumps_parallel.v
            NPSHr = self.speed.v**2 * (self.NPSH_a0.v + self.NPSH_a1.v*(Q_dot/self.speed.v) + self.NPSH_a2.v*(Q_dot/self.speed.v)**2)
            NPSH_meas = self.Pressure.v/rho/self.g - 1034000.0/rho/self.g #  NPSH in the simulation is relative to the pressure in the expansion tank (reason for subtracting 150 psi)
            if(NPSHr>NPSH_meas):
                self.cavitation.v = 1.0
            else:
                self.cavitation.v = 0.0
            return


        # -----------------------------------------------------------------------------------------------------------------------
        if (self.model.iteration == 0) or (self.model.iteration == 0 and self.model.time == self.model.settings.timestep):
            #  Set output values to the computed values from the previous timestep (don't want any manipulation without feedback)
            #  Determine pressure increase in pump
            rho = Inc.density(self.Fluid_ID.v, self.Temperature.v, 0.0)
            delta_P = self.g*self.speed.v**2*(
                rho*self.Pump_a0.v +
                self.Pump_a1.v/self.speed.v*(self.Mass_Flow.v/self.N_pumps_parallel.v) +
                self.Pump_a2.v/rho/self.speed.v**2*(self.Mass_Flow.v/self.N_pumps_parallel.v)**2
                )

            # Set the Initial Values of the Outputs (#,Value)
            self.m_dot_out.v = self.Mass_Flow.v
            self.P_out.v = self.Pressure.v+delta_P
            self.Temperature_out.v = self.Temperature.v
            self.mass_count_out.v = self.mass_count.v
            # self.cavitation.v #don't change

        else:
            if(self.Pump_Solver.v == 0.0):
                #  Compute Density of fluid
                rho = Inc.density(self.Fluid_ID.v, self.Temperature.v, 0.0)

                #  Compute head loss in system
                delta_P = self.P_out.v-(self.Pressure.v+self.error.v)
                H_L = delta_P/rho/9.81 * 3.28084 #  [ft]

                #  Compute new flow rate corresponding to head loss in system
                #  (solving a quadratic equation for the pump curve fit)
                A = self.Pump_a2.v
                B = self.Pump_a1.v*self.speed.v
                D = self.Pump_a0.v*(self.speed.v**2) - H_L

                discriminant = (B)**2 - 4*A*D
                if (discriminant>=0):
                    sol1 = (-B + np.sqrt(discriminant))/(2*A)
                    sol2 = (-B - np.sqrt(discriminant))/(2*A)
                else:
                    discriminant = 0
                    sol1 = (-B + np.sqrt(discriminant))/(2*A)
                    sol2 = (-B - np.sqrt(discriminant))/(2*A)

                Q_dot_new = max([sol1, sol2]) * 0.00006309019640343866 #  [m^3/s]
                m_dot_new = Q_dot_new*rho

                #  Update mass flow rate according to learning rate
                m_dot_adj = self.LR.v*m_dot_new + (1-self.LR.v)*self.m_dot_out.v/self.N_pumps_parallel.v

                P_out = (self.Pressure.v+delta_P)

                # Set the Outputs from this Model (#,Value)
                self.m_dot_out.v = m_dot_adj*self.N_pumps_parallel.v
                self.P_out.v = P_out
                self.Temperature_out.v = self.Temperature.v
                self.mass_count_out.v = self.mass_count.v
                # self.cavitation.v # don't change

            else:
                #  Determine pressure increase in pump
                rho = Inc.density(  self.Fluid_ID.v, self.Temperature.v, 0.0)
                delta_P = self.g*self.speed.v**2*(rho*self.Pump_a0.v +
                    self.Pump_a1.v/self.speed.v*(self.Mass_Flow.v/self.N_pumps_parallel.v) +
                    self.Pump_a2.v/rho/self.speed.v**2*(self.Mass_Flow.v/self.N_pumps_parallel.v)**2)

                # Set the Outputs from this Model (#,Value)
                self.m_dot_out.v = self.Mass_Flow.v
                self.P_out.v = self.Pressure.v + delta_P
                self.Temperature_out.v = self.Temperature.v
                self.mass_count_out.v = self.mass_count.v
                # self.cavitation.v # don't change
