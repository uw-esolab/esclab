"""HTF Tank-to-Tank Pump component model (Type 6032)."""

from esclab.simulate import Component
from esclab.components.esol_properties import Incompressible as Inc


class HTFTank2TankPump(Component):
    """
    Object: Tank2TankPump - Water
    Simulation Studio Model: ESOL6032-HTFTank2TankPump

    Author: AnnaGroeschel
    Editor:
    Date:    April 28, 2024
    last modified: April 28, 2024
    Converted from TRNSYS Type 6032
    """
    trnsys_type = "6032"

    # PARAMETERS (10 total; params 5-8 are commented out in the Fortran source)
    L_downcomer = Component.Parameter()      # Length of downcomer between tank and pump [m]
    PC_CoefA = Component.Parameter()         # Pump Curve Fit: PC_CoefA*(Q)**2 + PC_CoefB*(Q/pump_speed) + PC_CoefC/pump_speed**2
    PC_CoefB = Component.Parameter()         # Pump Curve Fit: PC_CoefA*(Q)**2 + PC_CoefB*(Q/pump_speed) + PC_CoefC/pump_speed**2
    PC_CoefC = Component.Parameter()         # Pump Curve Fit: PC_CoefA*(Q)**2 + PC_CoefB*(Q/pump_speed) + PC_CoefC/pump_speed**2
    # Eta_CoefA = Component.Parameter()      # Pump Efficiency Curve Fit: Eta_Coef_A*(Q/pump_speed)**4 + Eta_Coef_B*(Q/pump_speed)**3 + Eta_Coef_C*(Q/pump_speed)**2 + Eta_Coef_D*(Q/pump_speed)
    # Eta_CoefB = Component.Parameter()      # Pump Efficiency Curve Fit: Eta_Coef_A*(Q/pump_speed)**4 + Eta_Coef_B*(Q/pump_speed)**3 + Eta_Coef_C*(Q/pump_speed)**2 + Eta_Coef_D*(Q/pump_speed)
    # Eta_CoefC = Component.Parameter()      # Pump Efficiency Curve Fit: Eta_Coef_A*(Q/pump_speed)**4 + Eta_Coef_B*(Q/pump_speed)**3 + Eta_Coef_C*(Q/pump_speed)**2 + Eta_Coef_D*(Q/pump_speed)
    # Eta_CoefD = Component.Parameter()      # Pump Efficiency Curve Fit: Eta_Coef_A*(Q/pump_speed)**4 + Eta_Coef_B*(Q/pump_speed)**3 + Eta_Coef_C*(Q/pump_speed)**2 + Eta_Coef_D*(Q/pump_speed)
    Fluid_ID = Component.Parameter()         # Fluid ID used to find thermal properties (40 = Dowtherm A)
    pump_speed_rate = Component.Parameter()  # Rate of change of pump speed (pump ramp)

    # INPUTS (7 total)
    PUMP_ON = Component.Input()              # Operator Input of Pump Power signal (1 == ON)
    pump_speed = Component.Input()           # Operator Input of Pump Speed signal (0 - 1) [%base1]
    T_tank1 = Component.Input()             # Temperature of Tank prior to pump
    P_tank1 = Component.Input()             # Pressure of Tank prior to pump at top of level
    L_tank1 = Component.Input()             # Level of Tank prior to pump
    P_tank2 = Component.Input()             # Pressure of 2nd tank that HTF is sent to
    P_piping_losses = Component.Input()      # Pressure of the HTF prior to entering the 2nd tank to calculate system losses in piping network

    # OUTPUTS (11 total)
    m_dot_pump = Component.Output()          # mass flow rate leaving the pump [kg/s]
    Vol_dot_pump = Component.Output()        # Volumetric flow rate leaving the pump [m^3/s]
    P_pump_out = Component.Output()          # Pressure leaving the pump [Pa]
    T_pump_out = Component.Output()          # Temperature leaving the pump [K]
    Eta_pump = Component.Output()            # Pump Efficiency [% base 1]
    W_dot_pump = Component.Output()          # Power needed to run the pump [W]
    Point_1x = Component.Output()           # Previous Flow Rate - used to iterate to determine the correct flow for this timestep
    Point_1y = Component.Output()           # Error associated with previous flow rate - used to iterate to determine the correct flow for this timestep
    Point_2x = Component.Output()           # Current Flow Rate used - used to iterate to determine the correct flow for this timestep
    Pump_head = Component.Output()           # Pump Head [m]
    pump_speed_out = Component.Output()      # Current pump speed for this timestep (used to track pump speed ramp between timesteps)

    def presim_setup(self, **kwargs):
        rho = Inc.density(self.Fluid_ID.v, self.T_tank1.v, self.P_tank1.v)

        P_bottom = self.P_tank1.v + (self.L_tank1.v + self.L_downcomer.v) * rho * 9.81

        self.m_dot_pump.v = 0.0001

        # Set the Initial Values of the Outputs
        self.Vol_dot_pump.v = self.m_dot_pump.v / rho
        self.P_pump_out.v = P_bottom + self.PC_CoefC.v * self.pump_speed.v**2 * rho * 9.81
        self.T_pump_out.v = 400.0
        self.Eta_pump.v = 0.0
        self.W_dot_pump.v = 0.0
        self.Point_1x.v = 0.0
        self.Point_1y.v = 0.0
        self.Point_2x.v = 0.0
        self.Pump_head.v = self.PC_CoefC.v
        self.pump_speed_out.v = self.pump_speed.v

    def calculate(self):
        tol = 0.1
        Q_min = 0.0000000001
        LR = 0.5

        if self.PUMP_ON.v == 1.0:  # pump is on
            if self.model.iteration == 0:  # First Iteration of the timestep
                if not self.model.is_first_step:  # Not the first timestep of the simulation
                    # adjust pump speed if needed
                    pump_speed_o = self.pump_speed_out.v
                    ts = self.model.timestep * 3600.0  # Convert timestep from hr to s
                    if self.pump_speed.v == pump_speed_o:
                        pump_speed = pump_speed_o
                    elif self.pump_speed.v > pump_speed_o:
                        pump_speed = min(pump_speed_o + self.pump_speed_rate.v * ts, self.pump_speed.v)
                    else:
                        pump_speed = max(pump_speed_o - self.pump_speed_rate.v * ts, self.pump_speed.v)
                    # set current pump speed for this timestep
                    self.pump_speed_out.v = pump_speed
                    C = self.PC_CoefC.v * pump_speed**2
                    B = self.PC_CoefB.v * pump_speed
                    A = self.PC_CoefA.v
                    # Maximum flow allowed for given pump speed input
                    Q_max = max(
                        (-B + (B**2 - 4.0 * A * C)**0.5) / (2.0 * A),
                        (-B - (B**2 - 4.0 * A * C)**0.5) / (2.0 * A)
                    )
                    if Q_max <= 0.0:
                        Q_max = 0.00001
                    m_dot_pump = self.m_dot_pump.v  # Previous flow rate guessed
                    P_pump_out = self.P_pump_out.v  # Previous Pressure leaving the pump
                    rho = Inc.density(self.Fluid_ID.v, self.T_tank1.v, self.P_tank1.v)  # Density of HTF
                    Pump_head = P_pump_out / rho / 9.81 - self.P_tank1.v / rho / 9.81 - self.L_tank1.v - self.L_downcomer.v  # Previous pump head used in last iteration
                    System_headloss = (P_pump_out - self.P_piping_losses.v) / rho / 9.81  # Previous system headloss
                    Point_1x = m_dot_pump / rho
                    # Previous error
                    Point_1y = Pump_head - self.P_tank2.v / rho / 9.81 - System_headloss + self.P_tank1.v / rho / 9.81 + self.L_tank1.v + self.L_downcomer.v
                    if abs(Point_1y) <= tol:
                        Q_new = Point_1x
                    elif Point_1y > 0.0:  # increase flow rate
                        Q_new = min(Point_1x + 0.01, Q_max)  # Move to maximum flow pump can provide at given speed
                    else:  # decrease flow rate
                        Q_new = max(Point_1x - 0.01, Q_min)  # Move to minimum flow pump can provide
                    m_dot_pump = Q_new * rho
                    P_pump_in = self.P_tank1.v + rho * 9.81 * (self.L_tank1.v + self.L_downcomer.v)
                    P_pump_out = (A * Q_new**2 + B * Q_new + C) * rho * 9.81 + P_pump_in
                    pump_head = (P_pump_out - P_pump_in) / rho / 9.81
                    T_pump_out = self.T_tank1.v
                    Eta_pump = 0.7
                    # Eta_pump = max(Eta_CoefA * (Q_new)**4 + Eta_CoefB*pump_speed * (Q_new)**3 + Eta_CoefC*pump_speed**2*(Q_new)**2 + Eta_CoefD*Q_new*pump_speed**3, 0.1)
                    W_dot_pump = 0.0  # Correct Later
                    Point_2x = Q_new
                else:  # First timestep of simulation, no inputs to base off of
                    pump_speed = self.pump_speed_out.v
                    m_dot_pump = self.m_dot_pump.v
                    P_pump_out = self.P_pump_out.v
                    T_pump_out = self.T_pump_out.v
                    # TODO-NEEDS CONVERSION REVIEW: rho is used before being defined on the following line; ported faithfully from Fortran (possible bug in source)
                    pump_head = (P_pump_out - self.P_tank1.v) / rho / 9.81 - self.L_downcomer.v - self.L_tank1.v
                    rho = Inc.density(self.Fluid_ID.v, T_pump_out, P_pump_out)
                    Point_1x = m_dot_pump / rho
                    Point_1y = 0.0       # Not solved for yet
                    Point_2x = 0.0       # Not solved for yet
                    Eta_pump = 0.0       # Not Solved for yet
                    W_dot_pump = 0.0     # Not solved for yet
                    pump_head = (P_pump_out - self.P_tank1.v) / rho / 9.81 - self.L_downcomer.v - self.L_tank1.v
            elif (self.model.iteration == 1) and self.model.is_first_step:
                # (getTimestepIteration() == 1) and (Time == Timestep): first timestep, second iteration
                pump_speed = self.pump_speed_out.v
                C = self.PC_CoefC.v * pump_speed**2
                B = self.PC_CoefB.v * pump_speed
                A = self.PC_CoefA.v
                Q_max = max(
                    (-B + (B**2 - 4.0 * A * C)**0.5) / (2.0 * A),
                    (-B - (B**2 - 4.0 * A * C)**0.5) / (2.0 * A)
                )
                if Q_max <= 0.0:
                    Q_max = 0.00001
                Point_1x = self.Point_1x.v     # Flow used two iterations ago
                P_pump_out = self.P_pump_out.v  # Pump outlet pressure computed at last iteration
                pump_head = self.Pump_head.v
                rho = Inc.density(self.Fluid_ID.v, self.T_tank1.v, P_pump_out)
                # Compute head loss in system
                System_headloss = (P_pump_out - self.P_piping_losses.v) / rho / 9.81
                # New error for last flow sent through system
                Point_1y = self.P_tank1.v / rho / 9.81 + self.L_tank1.v + self.L_downcomer.v + pump_head - System_headloss - self.P_tank2.v / rho / 9.81
                if abs(Point_1y) <= tol:
                    Q_new = Point_1x
                elif Point_1y > 0.0:
                    Q_new = min(Point_1x + 0.05, Q_max)
                else:
                    Q_new = max(Point_1x - 0.05, Q_min)
                m_dot_pump = Q_new * rho
                P_pump_in = self.P_tank1.v + rho * 9.81 + (self.L_tank1.v + self.L_downcomer.v)
                pump_head = A * Q_new**2 + B * Q_new + C
                P_pump_out = P_pump_in + pump_head * rho * 9.81
                T_pump_out = self.T_tank1.v
                Eta_pump = 0.7
                # Eta_pump = max(Eta_CoefA * (Q_new)**4 + Eta_CoefB*pump_speed * (Q_new)**3 + Eta_CoefC*pump_speed**2*(Q_new)**2 + Eta_CoefD*Q_new*pump_speed**3, 0.1)
                W_dot_pump = 0.0  # Correct Later
                Point_1x = self.Point_2x.v
                # TODO-NEEDS CONVERSION REVIEW: Point_2y is referenced here but is never assigned in this branch; ported faithfully from Fortran (possible bug in source)
                Point_1y = Point_2y  # noqa: F821
                Point_2x = Q_new
            else:  # Third + iteration of the timestep of first timestep, or 2nd + iteration of any other timestep; enough info to iterate flow
                pump_speed = self.pump_speed_out.v
                C = self.PC_CoefC.v * pump_speed**2
                B = self.PC_CoefB.v * pump_speed
                A = self.PC_CoefA.v
                # Maximum flow allowed for given pump speed input
                Q_max = max(
                    (-B + (B**2 - 4.0 * A * C)**0.5) / (2.0 * A),
                    (-B - (B**2 - 4.0 * A * C)**0.5) / (2.0 * A)
                )
                if Q_max <= 0.0:
                    Q_max = 0.00001
                Point_1x = self.Point_1x.v     # Flow used two iterations ago
                Point_1y = self.Point_1y.v     # error received from Point_1x flow
                Point_2x = self.Point_2x.v     # flow used last iteration
                P_pump_out = self.P_pump_out.v  # Pump outlet pressure computed at last iteration
                pump_head = self.Pump_head.v
                rho = Inc.density(self.Fluid_ID.v, self.T_tank1.v, P_pump_out)
                # Compute head loss in system
                System_headloss = (P_pump_out - self.P_piping_losses.v) / rho / 9.81
                # New error for last flow sent through system
                Point_2y = self.P_tank1.v / rho / 9.81 + self.L_tank1.v + self.L_downcomer.v + pump_head - System_headloss - self.P_tank2.v / rho / 9.81
                # create line between points
                if Point_2x == Point_1x:
                    m = 0.0
                else:
                    m = (Point_2y - Point_1y) / (Point_2x - Point_1x)
                y_int = Point_1y - m * Point_1x
                if abs(Point_2y) <= tol:
                    Q_new = Point_2x
                elif m != 0.0:
                    if Point_2y >= 0.0:  # If Error is Positive, pump head guess is too high so flow must increase
                        Q_new = max(-y_int / m, Point_2x)   # Ensure new Q_dot does not decrease
                        Q_new = min(Q_new, Q_max)           # Ensure Q_dot_new is not more than pump can provide based on curve
                    else:
                        Q_new = min(-y_int / m, Point_2x)   # Ensure new Q_dot does not increase
                        Q_new = max(Q_new, Q_min)           # Ensure new Q_dot is not negative
                else:  # If Error is Negative, pump head guess is too low so flow must decrease
                    if Point_2y >= 0.0:  # Error is Positive and slope is equal to zero, pump head is too high so flow must increase
                        Q_new = min(Point_2x + 0.0001, Q_max)
                    else:  # Error is Negative and slope is equal to zero, pump head is too low so flow must decrease
                        Q_new = max(Point_2x - 0.0001, Q_min)
                Q_new = Point_2x + (Q_new - Point_2x) * LR
                m_dot_pump = Q_new * rho
                P_pump_in = self.P_tank1.v + rho * 9.81 + (self.L_tank1.v + self.L_downcomer.v)
                pump_head = A * Q_new**2 + B * Q_new + C
                P_pump_out = P_pump_in + pump_head * rho * 9.81
                T_pump_out = self.T_tank1.v
                Eta_pump = 0.7
                # Eta_pump = max(Eta_CoefA * (Q_new)**4 + Eta_CoefB*pump_speed * (Q_new)**3 + Eta_CoefC*pump_speed**2*(Q_new)**2 + Eta_CoefD*Q_new*pump_speed**3, 0.1)
                W_dot_pump = 0.0  # Correct Later
                Point_1x = Point_2x
                Point_1y = Point_2y
                Point_2x = Q_new
        else:  # pump is off, no flow through pump
            rho = Inc.density(self.Fluid_ID.v, self.T_tank1.v, self.P_tank1.v)
            m_dot_pump = Q_min * rho
            # TODO-NEEDS CONVERSION REVIEW: In Fortran, pump_speed is a local variable not set in the pump-off branch;
            # pump_speed_out.v (which tracks the ramped speed) is used here to match Fortran behavior
            P_pump_out = self.PC_CoefC.v * self.pump_speed_out.v**2 * 9.81 * rho
            T_pump_out = self.T_tank1.v
            W_dot_pump = 0.0
            Eta_pump = 0.0
            # Retain iteration tracking outputs from previous call (not updated when pump is off)
            Point_1x = self.Point_1x.v
            Point_1y = self.Point_1y.v
            Point_2x = self.Point_2x.v
            pump_head = self.Pump_head.v

        # Set the Outputs from this Model
        self.m_dot_pump.v = m_dot_pump          # mass flow rate leaving the pump [kg/s]
        # Vol_dot_pump is not recomputed in the main calculation loop; retains presim_setup value (ported faithfully from Fortran source)
        self.P_pump_out.v = P_pump_out          # Pressure leaving the pump [Pa]
        self.T_pump_out.v = T_pump_out          # Temperature leaving the pump [K]
        self.Eta_pump.v = Eta_pump              # Pump Efficiency [% base 1]
        self.W_dot_pump.v = W_dot_pump          # Power needed to run the pump [W]
        self.Point_1x.v = Point_1x             # Previous Flow Rate - used to iterate to determine the correct flow for this timestep
        self.Point_1y.v = Point_1y             # Error associated with previous flow rate - used to iterate to determine the correct flow for this timestep
        self.Point_2x.v = Point_2x             # Current Flow Rate used - used to iterate to determine the correct flow for this timestep
        self.Pump_head.v = pump_head            # Pump Head [m]
        # pump_speed_out is set inside the pump-on iteration logic; retained from previous timestep when pump is off
