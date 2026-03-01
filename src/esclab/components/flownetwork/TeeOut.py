"""Tee junction (splitting) component model (Type 4006)."""

from esclab.simulate import Component


class TeeOut(Component):
    """
    # Object: 4006-Tee_Out
    # Simulation Studio Model: ESOL4006-Tee_Out
    #

    # Author: Matt Tuman
    # Editor:
    # Date:     October 20, 2022
    # last modified: October 24, 2022
    """

    #    PARAMETERS
    Fluid_ID = Component.Parameter()        # Fluid ID: 40 for dowtherm A
    Solving_Method = Component.Parameter()  # Adjust solver

    #    INPUTS
    m_dot = Component.Input()         # Mass flow into tee [kg/s]
    Pressure = Component.Input()      # Pressure at inlet of the tee [Pa]
    Temperature = Component.Input()   # Temperature of fluid into tee [C]
    f = Component.Input()             # Fraction of flow of tee [0-1]
    P_1_out = Component.Input()       # Pressure at the end of loop 1 [Pa]
    P_2_out = Component.Input()       # Pressure at the end of loop 2 [Pa]
    m_dot_prev = Component.Input()    # Previous mass flow rate into the tee [kg/s]
    P_prev = Component.Input()        # Previous pressure at the inlet of the tee [Pa]
    error = Component.Input()         # Error in pressure due to a component with specified pressure
    Mass_Counter = Component.Input()  # Total amount of mass in the system up to this point

    #    OUTPUTS
    m_dot_1 = Component.Output()          # mass flow sent down branch 1
    P_1 = Component.Output()              # Pressure sent down branch 1
    Temp_1 = Component.Output()           # temperature sent down branch 1
    m_dot_2 = Component.Output()          # mass flow sent down branch 2
    P_2 = Component.Output()              # pressure sent down branch 2
    Temp_2 = Component.Output()           # temperature sent down branch 2
    f_out = Component.Output()            # previous flow fraction guessed
    m_dot_out = Component.Output()        # mass flow rate that entered the tee
    Pressure_out = Component.Output()     # pressure that entered the tee
    Mass_Counter_out = Component.Output() # mass counter output used for the expansion system type

    def presim_setup(self, **kwargs):
        pass

    def calculate(self):

        # Do All of the "Very First Call of the Simulation Manipulations" Here
        if self.model.is_first_step:
            # Use initial guess of fraction of flow to compute mass flow rates out of tee
            m_dot_1 = self.m_dot.v * self.f.v
            m_dot_2 = self.m_dot.v - m_dot_1

            # Set the Initial Values of the Outputs (#,Value)
            self.m_dot_1.v = m_dot_1        # mass flow sent down branch 1
            self.P_1.v = self.Pressure.v    # Pressure sent down branch 1
            self.Temp_1.v = self.Temperature.v  # temperature sent down branch 1
            self.m_dot_2.v = m_dot_2        # mass flow sent down branch 2
            self.P_2.v = self.Pressure.v    # pressure sent down branch 2
            self.Temp_2.v = self.Temperature.v  # temperature sent down branch 2
            self.f_out.v = self.f.v         # previous flow fraction guessed
            self.m_dot_out.v = self.m_dot.v  # mass flow rate that entered the tee
            self.Pressure_out.v = self.Pressure.v  # pressure that entered the tee
            self.Mass_Counter_out.v = self.Mass_Counter.v  # mass counter output used for the expansion system type
            return

        # Perform Any "After Convergence" Manipulations That May Be Required at the End of Each Timestep
        if self.model.is_converged:
            self.f.v = self.model.iteration
            return

        # TYPE PERFORMS CALCULATIONS HERE: -----------------------------------------------------------------------

        #----------------------------------------------------------------------------------------------------------
        #
        #-------------------- SOLVE USING LINEARIZED SOLUTION----------------------------------------------------
        #
        if self.Solving_Method.v == 0.0:

            if self.model.iteration > 0 or self.model.time > self.model.settings.timestep:
                # Compute K values
                m_dot_1 = self.f.v * self.m_dot_prev.v
                m_dot_2 = (1 - self.f.v) * self.m_dot_prev.v
                K_T1 = (self.P_prev.v - self.P_1_out.v) / m_dot_1 ** 2
                K_T2 = (self.P_prev.v - self.P_2_out.v) / m_dot_2 ** 2

                # Compute new f value
                f_new = K_T2 * m_dot_2 / (K_T2 * m_dot_2 + K_T1 * m_dot_1)

                # Update f value
                # f = (f_new + f)/2.d0
                LR = 0.5
                self.f.v = f_new * LR + self.f.v * (1 - LR)

                # Set Outputs
                self.m_dot_1.v = self.f.v * self.m_dot.v  # m_dot_1
                self.P_1.v = self.Pressure.v               # P_1
                self.Temp_1.v = self.Temperature.v         # Temp_1
                self.m_dot_2.v = (1 - self.f.v) * self.m_dot.v  # m_dot_2
                self.P_2.v = self.Pressure.v               # P_2
                self.Temp_2.v = self.Temperature.v         # Temp_2
                self.f_out.v = self.f.v
                self.m_dot_out.v = self.m_dot.v
                self.Pressure_out.v = self.Pressure.v
                self.Mass_Counter_out.v = self.Mass_Counter.v
                return

            # IF first timestep in the simulation
            else:
                self.m_dot_1.v = self.f.v * self.m_dot.v  # m_dot_1
                self.P_1.v = self.Pressure.v               # P_1
                self.Temp_1.v = self.Temperature.v         # Temp_1
                self.m_dot_2.v = (1 - self.f.v) * self.m_dot.v  # m_dot_2
                self.P_2.v = self.Pressure.v               # P_2
                self.Temp_2.v = self.Temperature.v         # Temp_2
                self.f_out.v = self.f.v
                self.m_dot_out.v = self.m_dot.v
                self.Pressure_out.v = self.Pressure.v
                self.Mass_Counter_out.v = self.Mass_Counter.v
                return

        #----------------------------------------------------------------------------------------------------------
        #
        #-------------------- SOLVE USING DIRECT QUADRATIC SOLUTION----------------------------------------------
        #
        elif self.Solving_Method.v == 1.0:
            # Set learning rate
            LR = 0.41

            # Compute k_t values using output from the previous iteration
            # TODO-NEEDS CONVERSION REVIEW: Most of Solving_Method==1 logic is commented out in the original
            # Fortran. K_T1/K_T2 are computed below but m_dot_1, m_dot_2, f_new_adj are never assigned in
            # this branch; the final SetOutputValue block at the end of the subroutine uses their uninitialized
            # values (effectively 0 in Fortran). In Python this will raise NameError.
            error_term = 3.0
            if error_term == 1.0:
                K_T1 = (self.P_prev.v - (self.P_1_out.v + self.error.v)) / self.m_dot_1.v ** 2
                K_T2 = (self.P_prev.v - self.P_2_out.v) / self.m_dot_2.v ** 2
            elif error_term == 2.0:
                K_T1 = (self.P_prev.v - self.P_1_out.v) / self.m_dot_1.v ** 2
                K_T2 = (self.P_prev.v - (self.P_2_out.v + self.error.v)) / self.m_dot_2.v ** 2
            else:
                K_T1 = (self.P_prev.v - self.P_1_out.v) / self.m_dot_1.v ** 2
                K_T2 = (self.P_prev.v - self.P_2_out.v) / self.m_dot_2.v ** 2

        #----------------------------------------------------------------------------------------------------------
        #
        #--------------------------- GRADIENT DESCENT APPROACH ---------------------------------------------------
        ############################ WOULD NOT RECOMMEND USING ###################################################

        elif self.Solving_Method.v == 2.0:
            Beta = 0.01
            # First iteration: allow to solve for new mass flow rates
            if self.model.is_first_iteration:
                m_dot_1 = self.f.v * self.m_dot.v
                m_dot_2 = (1 - self.f.v) * self.m_dot.v
                f_new = self.f.v
                # TODO-NEEDS CONVERSION REVIEW: f_new_adj is not assigned in the first-iteration path;
                # the original Fortran uses its uninitialized value in the final SetOutputValue block.

            # Second iteration or more
            else:
                # Load in the fraction, inlet pressure, and mass flow rate from the first iteration
                m_dot_1 = self.f.v * self.m_dot_prev.v
                m_dot_2 = (1 - self.f.v) * self.m_dot_prev.v

                # Compute k_t values using output from the first iteration
                if m_dot_1 <= 0:
                    K_T1 = (self.P_prev.v - self.P_1_out.v) / (0.00000001 ** 2)
                else:
                    K_T1 = (self.P_prev.v - self.P_1_out.v) / m_dot_1 ** 2

                if m_dot_2 <= 0:
                    K_T2 = (self.P_prev.v - self.P_2_out.v) / (0.000000001 ** 2)
                else:
                    K_T2 = (self.P_prev.v - self.P_2_out.v) / m_dot_2 ** 2

                # Compute the new fraction value: f_new
                dPdf = ((K_T2 - K_T1) * self.f.v - K_T2) * ((K_T2 - K_T1) * self.f.v ** 2 - 2 * K_T2 * self.f.v + K_T2)
                f_new = self.f.v - Beta * dPdf

                if self.P_1_out.v < 0:
                    f_new = 0.0000000001
                if self.P_2_out.v < 0:
                    f_new = 0.9999999999

                # TODO-NEEDS CONVERSION REVIEW: LR is not defined in Solving_Method==2; in the original Fortran
                # LR is declared as a shared local and would have whatever value was last written to it
                # (from Solving_Method==1's LR=0.41 assignment, if that branch ran previously). This is a bug
                # in the original Fortran code.
                f_new_adj = LR * f_new + (1 - LR) * self.f_out.v
                m_dot_1 = self.m_dot.v * f_new_adj
                m_dot_2 = self.m_dot.v * (1 - f_new_adj)

        # IF USING ESOL4050-Parallel-Flow-Solver
        elif self.Solving_Method.v == 3.0:
            # Read in flow value computed by ESOL4050
            # f_new_adj = f  # Added by C.Volkwein
            # Set Outputs
            self.m_dot_1.v = self.f.v * self.m_dot.v  # m_dot_1
            self.P_1.v = self.Pressure.v               # P_1
            self.Temp_1.v = self.Temperature.v         # Temp_1
            self.m_dot_2.v = (1 - self.f.v) * self.m_dot.v  # m_dot_2
            self.P_2.v = self.Pressure.v               # P_2
            self.Temp_2.v = self.Temperature.v         # Temp_2
            self.f_out.v = self.f.v
            self.m_dot_out.v = self.m_dot.v
            self.Pressure_out.v = self.Pressure.v
            self.Mass_Counter_out.v = self.Mass_Counter.v
            return

        #----------------------------------------------------------------------------------------------------------
        # Set the Outputs from this Model (#,Value)
        self.m_dot_1.v = m_dot_1            # m_dot_1
        self.P_1.v = self.Pressure.v        # P_1
        self.Temp_1.v = self.Temperature.v  # Temp_1
        self.m_dot_2.v = m_dot_2            # m_dot_2
        self.P_2.v = self.Pressure.v        # P_2
        self.Temp_2.v = self.Temperature.v  # Temp_2
        self.f_out.v = f_new_adj
        self.m_dot_out.v = self.m_dot.v
        self.Pressure_out.v = self.Pressure.v
        self.Mass_Counter_out.v = self.Mass_Counter.v
