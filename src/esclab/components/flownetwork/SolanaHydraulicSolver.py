"""Solana Hydraulic Solver component model (Type 6031)."""

import numpy as np

from esclab.simulate import Component
from esclab.components.esol_properties import Incompressible as Inc


class SolanaHydraulicSolver(Component):
    """
    Object: ESOL6031-SolanaHydraulicModel
    Simulation Studio Model: ESOL6031-SolanaHydraulicModel

    Author:
    Editor:
    Date:    January 02, 2025
    last modified: January 02, 2025

    Solana Hydraulic Solver: Including TES Charging/Discharging Pathways, pumps, and expansion system

              |------------------
              |                 | (0)
              | (1)             |
              |                 |
              ^ (f)             |
             / \\                |
       (2)  /   \\               |
           /     \\ (3)          |
          /       \\             |
         |         |            |
         |         |            |
         ^         ^            |
        / \\       / \\           |
       /   \\   __/   \\          |
  (4) |     | |      |          |
      |     \\ /      |          |
      |     / \\   (7)|          |
      |    |   |     |          |
      |    |   |_    |          |
      | (6)|     |   |          |
      |   _|  (5)|   |          |
      |   |      |   |          |
       \\ /        \\ /           |
        U          U            |
        |          |            |
        |          |            |
        L__________|____________|
        \\          /
         \\        /
      (8) \\      / (9)
           \\    /
            \\  /
              U
              |
              |
    """

    #    PARAMETERS
    Solver = Component.Parameter()                  #
    Fluid_ID = Component.Parameter()               # Fluid ID used to obtain fluid properties (40 = Dowtherm A)
    N_pumps = Component.Parameter()                # Number of pumps in parallel used in the solar field
    pump_speed_rate = Component.Parameter()        # Rate at which pump speed can change over a second (0.05 = 5% a second)
    Pump_Coef_A = Component.Parameter()            # Pump Curve Coefficient A - PC_A*Q**2 + PC_B*pump_speed*Q + PC_C*pump_speed**2
    Pump_Coef_B = Component.Parameter()            # Pump Curve Coefficient B - PC_A*Q**2 + PC_B*pump_speed*Q + PC_C*pump_speed**2
    Pump_Coef_C = Component.Parameter()            # Pump Curve Coefficient C - PC_A*Q**2 + PC_B*pump_speed*Q + PC_C*pump_speed**2
    Pump_Flow_guess = Component.Parameter()        # Pump Flow Guess - used to help solve the first iteration factor
    flowfrac1_guess = Component.Parameter()        # Flow fraction 1 guess
    flowfrac2_guess = Component.Parameter()        # Flow fraction 2 guess
    flowfrac3_guess = Component.Parameter()        # flow fraction 3 guess

    #    INPUTS
    # Note: Fortran SetNumberofInputs(28); inputs 1-11 are unused in this type body
    # (m_dot state is carried via outputs, not external inputs)
    P_pump_in = Component.Input()     # input 12
    P0 = Component.Input()            # input 13
    P1 = Component.Input()            # input 14
    P2 = Component.Input()            # input 15
    P3 = Component.Input()            # input 16
    P4 = Component.Input()            # input 17
    P5 = Component.Input()            # input 18
    P6 = Component.Input()            # input 19
    P7 = Component.Input()            # input 20
    P8 = Component.Input()            # input 21
    P9 = Component.Input()            # input 22
    P10 = Component.Input()           # input 23
    P_exp = Component.Input()         # input 24
    P_pump = Component.Input()        # input 25
    pump_speed_i = Component.Input()  # input 26
    T_1 = Component.Input()           # input 27
    turbine_on = Component.Input()    # input 28

    #    OUTPUTS
    m_dot_pump = Component.Output()   # 1  m_dot_pump
    m_dot_0 = Component.Output()      # 2  m_dot_0
    m_dot_1 = Component.Output()      # 3  m_dot_1
    m_dot_2 = Component.Output()      # 4  m_dot_2
    m_dot_3 = Component.Output()      # 5  m_dot_3
    m_dot_4 = Component.Output()      # 6  m_dot_4
    m_dot_5 = Component.Output()      # 7  m_dot_5
    m_dot_6 = Component.Output()      # 8  m_dot_6
    m_dot_7 = Component.Output()      # 9  m_dot_7
    m_dot_8 = Component.Output()      # 10 m_dot_8
    m_dot_9 = Component.Output()      # 11 m_dot_9
    pump_speed = Component.Output()   # 12 pump_speed (s)

    def presim_setup(self, **kwargs):
        super().presim_setup(**kwargs)

        # Tell the TRNSYS Engine How This Type Works
        # SetNumberofParameters(11)
        # SetNumberofInputs(28)
        # SetNumberofDerivatives(0)
        # SetNumberofOutputs(12)
        # SetIterationMode(1)
        # SetNumberStoredVariables(0,0)
        # SetNumberofDiscreteControls(0)

        # flow fraction 4 guess (hardcoded, not a parameter in Fortran)
        flowfrac4_guess = 0.5

        self.m_dot_1.v = self.flowfrac1_guess.v * self.Pump_Flow_guess.v
        m_dot_0 = (1.0 - self.flowfrac1_guess.v) * self.Pump_Flow_guess.v
        m_dot_2 = self.flowfrac2_guess.v * self.m_dot_1.v
        m_dot_3 = (1.0 - self.flowfrac2_guess.v) * self.m_dot_1.v
        m_dot_4 = self.flowfrac3_guess.v * m_dot_2
        m_dot_5 = (1.0 - self.flowfrac3_guess.v) * m_dot_2
        m_dot_6 = flowfrac4_guess * m_dot_3
        m_dot_7 = (1.0 - flowfrac4_guess) * m_dot_3
        m_dot_8 = m_dot_4 + m_dot_6 + m_dot_0
        m_dot_9 = m_dot_5 + m_dot_7

        # Set the Initial Values of the Outputs (#,Value)
        self.m_dot_pump.v = self.Pump_Flow_guess.v   # 1  m_dot_pump
        self.m_dot_0.v = m_dot_0                      # 2  m_dot_0
        # m_dot_1 already set above                   # 3  m_dot_1
        self.m_dot_2.v = m_dot_2                      # 4  m_dot_2
        self.m_dot_3.v = m_dot_3                      # 5  m_dot_3
        self.m_dot_4.v = m_dot_4                      # 6  m_dot_4
        self.m_dot_5.v = m_dot_5                      # 7  m_dot_5
        self.m_dot_6.v = m_dot_6                      # 8  m_dot_6
        self.m_dot_7.v = m_dot_7                      # 9  m_dot_7
        self.m_dot_8.v = m_dot_8                      # 10 m_dot_8
        self.m_dot_9.v = m_dot_9                      # 11 m_dot_9
        self.pump_speed.v = self.pump_speed_i.v       # 12 pump_speed

        return

    def calculate(self):
        super().calculate()

        g = 9.81
        T_1 = self.T_1.v
        if T_1 <= 0.0:
            T_1 = 400.0
        # Density from ESOL6015_myfunctions; mapped to Incompressible.density
        rho = Inc.density(self.Fluid_ID.v, T_1, self.P_pump_in.v)

        # learning rate assignment
        LR = 0.25
        # if (turbine_on != 0.0):
        #   LR = LR / 2

        if self.model.timestep_iteration == 0:
            if not self.model.is_first_step:  # Do not adjust pump speed at the beginning of the simulation
                # Convert timestep from hr to s
                ts = self.model.timestep * 3600.0
                if self.pump_speed.v == self.pump_speed_i.v:
                    pass  # s already equals self.pump_speed.v
                elif self.pump_speed_i.v > self.pump_speed.v:
                    self.pump_speed.v = min(self.pump_speed.v + self.pump_speed_rate.v * ts, self.pump_speed_i.v)
                else:
                    self.pump_speed.v = max(self.pump_speed.v - self.pump_speed_rate.v * ts, self.pump_speed_i.v)
            # else: s already equals self.pump_speed.v from previous timestep

        # s is always read from self.pump_speed.v (set above or carried from previous call)
        s = self.pump_speed.v

        if self.model.timestep_iteration > 0 or not self.model.is_first_step:

            # Compute K values based on previous iteration
            K0 = abs((self.P1.v - self.P0.v) / self.m_dot_0.v ** 2.0)
            K2 = abs((self.P1.v - self.P2.v) / self.m_dot_2.v ** 2.0)
            K3 = abs((self.P1.v - self.P3.v) / self.m_dot_3.v ** 2.0)
            K4 = abs((self.P2.v - self.P4.v) / self.m_dot_4.v ** 2.0)
            K5 = abs((self.P2.v - self.P5.v) / self.m_dot_5.v ** 2.0)
            K6 = abs((self.P3.v - self.P6.v) / self.m_dot_6.v ** 2.0)
            K7 = abs((self.P3.v - self.P7.v) / self.m_dot_7.v ** 2.0)
            K8 = abs((self.P4.v - self.P8.v) / self.m_dot_8.v ** 2.0)
            K9 = abs((self.P7.v - self.P9.v) / self.m_dot_9.v ** 2.0)
            K_fin = abs((self.P8.v - self.P10.v) / self.m_dot_pump.v ** 2.0)
            K_in = abs((self.P_exp.v - self.P_pump.v) / self.m_dot_pump.v ** 2.0)

            # Calculate K_hat
            Kp_tot = -K_in - K_fin + g / rho / self.N_pumps.v ** 2.0 * self.Pump_Coef_A.v
            Kp_hat = Kp_tot * abs(self.m_dot_pump.v) + g * self.Pump_Coef_B.v * s / self.N_pumps.v
            K0_hat = K0 * abs(self.m_dot_0.v)
            K2_hat = K2 * abs(self.m_dot_2.v)
            K3_hat = K3 * abs(self.m_dot_3.v)
            K4_hat = K4 * abs(self.m_dot_4.v)
            K5_hat = K5 * abs(self.m_dot_5.v)
            K6_hat = K6 * abs(self.m_dot_6.v)
            K7_hat = K7 * abs(self.m_dot_7.v)
            K8_hat = K8 * abs(self.m_dot_8.v)
            K9_hat = K9 * abs(self.m_dot_9.v)

            # Create A matrix
            A = np.zeros((11, 11))  # Set all values as zero

            # Mass Balances
            A[0, 0] = 1.0
            A[0, 1] = -1.0
            A[0, 2] = -1.0

            A[1, 2] = 1.0
            A[1, 3] = -1.0
            A[1, 4] = -1.0

            A[2, 3] = 1.0
            A[2, 5] = -1.0
            A[2, 6] = -1.0

            A[3, 4] = 1.0
            A[3, 7] = -1.0
            A[3, 8] = -1.0

            A[4, 6] = -1.0
            A[4, 8] = -1.0
            A[4, 10] = 1.0

            A[5, 1] = -1.0
            A[5, 5] = -1.0
            A[5, 7] = -1.0
            A[5, 9] = 1.0

            # Pressure Relationships
            A[6, 3] = K2_hat
            A[6, 4] = -K3_hat
            A[6, 5] = K4_hat
            A[6, 7] = -K6_hat

            A[7, 3] = K2_hat
            A[7, 4] = -K3_hat
            A[7, 6] = K5_hat
            A[7, 8] = -K7_hat

            A[8, 1] = -K0_hat
            A[8, 3] = K2_hat
            A[8, 5] = K4_hat

            A[9, 3] = K2_hat
            A[9, 4] = -K3_hat
            A[9, 5] = K4_hat
            A[9, 8] = -K7_hat
            A[9, 9] = K8_hat
            A[9, 10] = -K9_hat

            A[10, 0] = Kp_hat
            A[10, 3] = -K2_hat
            A[10, 5] = -K4_hat
            A[10, 9] = -K8_hat

            # Create b array
            b = np.zeros(11)
            b[10] = -g * s ** 2.0 * rho * self.Pump_Coef_C.v

            # Invert matrix and compute mass flow
            # Replaces Fortran matrixinv/matmul with numpy.linalg.solve
            # call find_66_matrix(Time, CurrentUnit, getTimestepIteration(), A)
            # call find_66_matrix(Time, -CurrentUnit, getTimestepIteration(), A_inv)
            m_dots = np.linalg.solve(A, b)

            # Update mass flow rates
            if self.model.timestep_iteration == 1 and self.model.is_first_step:
                # total mass flow will never be approaching zero
                self.m_dot_pump.v = m_dots[0] * LR + self.m_dot_pump.v * (1.0 - LR)
                self.m_dot_0.v = m_dots[1] * LR + self.m_dot_0.v * (1.0 - LR)
                self.m_dot_1.v = m_dots[2] * LR + self.m_dot_1.v * (1.0 - LR)
                self.m_dot_2.v = m_dots[3] * LR + self.m_dot_2.v * (1.0 - LR)
                self.m_dot_3.v = m_dots[4] * LR + self.m_dot_3.v * (1.0 - LR)
                self.m_dot_4.v = m_dots[5] * LR + self.m_dot_4.v * (1.0 - LR)
                self.m_dot_5.v = m_dots[6] * LR + self.m_dot_5.v * (1.0 - LR)
                self.m_dot_6.v = m_dots[7] * LR + self.m_dot_6.v * (1.0 - LR)
                self.m_dot_7.v = m_dots[8] * LR + self.m_dot_7.v * (1.0 - LR)
                self.m_dot_8.v = m_dots[9] * LR + self.m_dot_8.v * (1.0 - LR)
                self.m_dot_9.v = m_dots[10] * LR + self.m_dot_9.v * (1.0 - LR)
            else:
                self.m_dot_pump.v = m_dots[0] * LR + self.m_dot_pump.v * (1.0 - LR)

                # if no flow is entering TES piping after leaving pump
                if m_dots[1] < 0.001 and m_dots[1] > 0.0:
                    m_hold = np.log10(m_dots[1]) * LR + np.log10(self.m_dot_0.v) * (1.0 - LR)
                    self.m_dot_0.v = 10 ** m_hold
                else:
                    self.m_dot_0.v = m_dots[1] * LR + self.m_dot_0.v * (1.0 - LR)

                # if flow entering the field is shut off
                if m_dots[2] < 0.001 and m_dots[2] > 0.0:
                    m_hold = np.log10(m_dots[2]) * LR + np.log10(self.m_dot_1.v) * (1.0 - LR)
                    self.m_dot_1.v = 10 ** m_hold
                else:
                    self.m_dot_1.v = m_dots[2] * LR + self.m_dot_1.v * (1.0 - LR)

                # If flow through near field is shut off
                if m_dots[3] < 0.001 and m_dots[3] > 0.0:
                    m_hold = np.log10(m_dots[3]) * LR + np.log10(self.m_dot_2.v) * (1.0 - LR)
                    self.m_dot_2.v = 10 ** m_hold
                else:
                    self.m_dot_2.v = m_dots[3] * LR + self.m_dot_2.v * (1.0 - LR)

                # If flow through the far field is shut off
                if m_dots[4] < 0.001 and m_dots[4] > 0.0:
                    m_hold = np.log10(m_dots[4]) * LR + np.log10(self.m_dot_3.v) * (1.0 - LR)
                    self.m_dot_3.v = 10 ** m_hold
                else:
                    self.m_dot_3.v = m_dots[4] * LR + self.m_dot_3.v * (1.0 - LR)

                # IF near-field to PB line is shut off
                if m_dots[5] < 0.001 and m_dots[5] > 0.0:
                    m_hold = np.log10(m_dots[5]) * LR + np.log10(self.m_dot_4.v) * (1.0 - LR)
                    self.m_dot_4.v = 10 ** m_hold
                else:
                    self.m_dot_4.v = m_dots[5] * LR + self.m_dot_4.v * (1.0 - LR)

                # If near-field PB bypass is shut
                if m_dots[6] < 0.001 and m_dots[6] > 0.0:
                    m_hold = np.log10(m_dots[6]) * LR + np.log10(self.m_dot_5.v) * (1.0 - LR)
                    self.m_dot_5.v = 10 ** m_hold
                else:
                    self.m_dot_5.v = m_dots[6] * LR + self.m_dot_5.v * (1.0 - LR)

                # If far-field to PB line is shut
                if m_dots[7] < 0.001 and m_dots[7] > 0.0:
                    m_hold = np.log10(m_dots[7]) * LR + np.log10(self.m_dot_6.v) * (1.0 - LR)
                    self.m_dot_6.v = 10 ** m_hold
                else:
                    self.m_dot_6.v = m_dots[7] * LR + self.m_dot_6.v * (1.0 - LR)

                # If far-field PB bypass is shut
                if m_dots[8] < 0.001 and m_dots[8] > 0.0:
                    m_hold = np.log10(m_dots[8]) * LR + np.log10(self.m_dot_7.v) * (1.0 - LR)
                    self.m_dot_7.v = 10 ** m_hold
                else:
                    self.m_dot_7.v = m_dots[8] * LR + self.m_dot_7.v * (1.0 - LR)

                # If PB line is shut
                if m_dots[9] < 0.001 and m_dots[9] > 0.0:
                    m_hold = np.log10(m_dots[9]) * LR + np.log10(self.m_dot_8.v) * (1.0 - LR)
                    self.m_dot_8.v = 10 ** m_hold
                else:
                    self.m_dot_8.v = m_dots[9] * LR + self.m_dot_8.v * (1.0 - LR)

                # If PB bypass line is shut
                if m_dots[10] < 0.001 and m_dots[10] > 0.0:
                    m_hold = np.log10(m_dots[10]) * LR + np.log10(self.m_dot_9.v) * (1.0 - LR)
                    self.m_dot_9.v = 10 ** m_hold
                else:
                    self.m_dot_9.v = m_dots[10] * LR + self.m_dot_9.v * (1.0 - LR)

        # Set the Outputs from this Model (#,Value)
        # 1  self.m_dot_pump.v  - m_dot_pump
        # 2  self.m_dot_0.v     - m_dot_0
        # 3  self.m_dot_1.v     - m_dot_1
        # 4  self.m_dot_2.v     - m_dot_2
        # 5  self.m_dot_3.v     - m_dot_3
        # 6  self.m_dot_4.v     - m_dot_4
        # 7  self.m_dot_5.v     - m_dot_5
        # 8  self.m_dot_6.v     - m_dot_6
        # 9  self.m_dot_7.v     - m_dot_7
        # 10 self.m_dot_8.v     - m_dot_8
        # 11 self.m_dot_9.v     - m_dot_9
        # 12 self.pump_speed.v  - pump_speed (set in timestep_iteration == 0 block above)
        # Call SetOutputValue(12, 0.d0)  ! pump_speed  (commented out in original)
        # call find_66_outputs(Time, CurrentUnit, getTimestepIteration(), end - start,
        #   [m_dot_pump, m_dot_0, m_dot_1, m_dot_2, m_dot_3, m_dot_4, m_dot_5,
        #    m_dot_6, m_dot_7, m_dot_8, m_dot_9, kappa(A, A_inv, 11)])

        return
