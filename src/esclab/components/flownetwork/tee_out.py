"""Tee junction component model (Type 4006)."""

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
    Fluid_ID = Component.Parameter()           # Fluid ID: 40 for dowtherm A
    Solving_Method = Component.Parameter()     # Adjust solver

    #    INPUTS
    m_dot = Component.Input()              # Mass flow into tee [kg/s]
    Pressure = Component.Input()           # Pressure at inlet of the tee [Pa]
    Temperature = Component.Input()        # Temperature of fluid into tee [C]
    f = Component.Input()                  # Fraction of flow of tee [0-1]
    P_1_out = Component.Input()            # Pressure at the end of loop 1 [Pa]
    P_2_out = Component.Input()            # Pressure at the end of loop 2 [Pa]
    m_dot_prev = Component.Input()         # Previous mass flow rate into the tee [kg/s]
    P_prev = Component.Input()             # Previous pressure at the inlet of the tee [Pa]
    error = Component.Input()              # Error in pressure due do a component with specified pressure
    Mass_Counter = Component.Input()       # Total amount of mass in the system up to this point

    #    OUTPUTS
    m_dot_1 = Component.Output()
    P_1 = Component.Output()
    Temp_1 = Component.Output()
    m_dot_2 = Component.Output()
    P_2 = Component.Output()
    Temp_2 = Component.Output()
    f_out = Component.Output()
    m_dot_out = Component.Output()
    Pressure_out = Component.Output()
    Mass_Counter_out = Component.Output()

    def presim_setup(self, **kwargs):
        pass

    def calculate(self):

        # Initialize output values for the first time step
        if self.model.is_first_step:
            # Use initial guess of fraction of flow to compute mass flow rates out of tee
            self.m_dot_1.v = self.m_dot.v*self.f.v
            self.m_dot_2.v = self.m_dot.v - self.m_dot_1.v

            # Set the Initial Values of the Outputs (#,Value)
            self.P_1.v = self.Pressure.v                   # Pressure sent down branch 1
            self.Temp_1.v = self.Temperature.v             # temperature sent down branch 1
            self.P_2.v = self.Pressure.v                   # pressure sent down branch 2
            self.Temp_2.v = self.Temperature.v             # temperature sent down branch 2
            self.f_out.v = self.f.v                        # previous flow fraction guessed
            self.m_dot_out.v = self.m_dot.v                # mass flow rate that entered the tee
            self.Pressure_out.v = self.Pressure.v          # pressure that entered the tee
            self.Mass_Counter_out.v = self.Mass_Counter.v  # mass counter output used for the expansion system type
            return

        #-----------------------------------------------------------------------------------------------------------------------
        # -----------  iteration calculations ----------------------
        #-----------------------------------------------------------------------------------------------------------------------

        #-------------------- SOLVE USING LINEARIZED SOLUTION-------------------------------------------------------------
        #
        if (self.Solving_Method.v == 0.0):
            if (self.model.iteration > 0) or (self.model.time > self.model.settings.timestep):
                # Compute K values
                self.m_dot_1.v = self.f.v*self.m_dot_prev.v
                self.m_dot_2.v = (1-self.f.v)*self.m_dot_prev.v
                K_T1 = (self.P_prev.v - self.P_1_out.v)/(self.m_dot_1.v**2)
                K_T2 = (self.P_prev.v - self.P_2_out.v)/(self.m_dot_2.v**2)

                # Compute new f value
                f_new = K_T2*self.m_dot_2.v/(K_T2*self.m_dot_2.v + K_T1*self.m_dot_1.v)

                # Update f value
                #f = (f_new + f)/2.0
                LR = 0.5
                self.f.v = f_new*LR + self.f.v*(1-LR)

            # Set Outputs
            m_dot_1_new = self.f.v*self.m_dot.v
            m_dot_2_new = (1-self.f.v)*self.m_dot.v
            f_new = self.f.v

        #-----------------------------------------------------------------------------------------------------------------------
        #
        #-------------------- SOLVE USING DIRECT QUADRATIC SOLUTION-------------------------------------------------------------
        #
        elif (self.Solving_Method.v == 1.0):
            # Set learning rate
            LR = .41

            ### Compute k_t values using output from the previous iteration
            error_term = 3.0
            if(error_term == 1.0):
                K_T1 = (self.P_prev.v - (self.P_1_out.v+self.error.v)) / (self.m_dot_1.v**2)
                K_T2 = (self.P_prev.v - self.P_2_out.v) / (self.m_dot_2.v**2)
            elif(error_term == 2.0):
                K_T1 = (self.P_prev.v - self.P_1_out.v) / (self.m_dot_1.v**2)
                K_T2 = (self.P_prev.v - (self.P_2_out.v+self.error.v)) / (self.m_dot_2.v**2)
            else:
                K_T1 = (self.P_prev.v - self.P_1_out.v) / (self.m_dot_1.v**2)
                K_T2 = (self.P_prev.v - self.P_2_out.v) / (self.m_dot_2.v**2)

        #-----------------------------------------------------------------------------------------------------------------------
        #
        #--------------------------- GRADIENT DESCENT APPROACH -----------------------------------------------------------------
        ############################ WOULD NOT RECCOMEND USING ###############################

        elif (self.Solving_Method.v == 2.0):
            Beta = 0.01
            # First iteration: allow to solve for new mass flow rates
            if self.model.is_first_iteration:
                self.m_dot_1.v = self.f.v*self.m_dot.v
                self.m_dot_2.v = (1-self.f.v)*self.m_dot.v
                f_new = self.f.v

            else:
                # Second iteration or more
                # Load in the fraction, inlet pressure, and mass flow rate from the first iteration
                m_dot_1_new = self.f.v*self.m_dot_prev.v
                m_dot_2_new = (1-self.f.v)*self.m_dot_prev.v

                # Compute k_t values using output from the first iteration
                if (m_dot_1_new <= 0) :
                    K_T1 = (self.P_prev.v - self.P_1_out.v) / (.00000001**2)
                else:
                    K_T1 = (self.P_prev.v - self.P_1_out.v) / (m_dot_1_new**2)

                if (m_dot_2_new <= 0) :
                    K_T2 = (self.P_prev.v - self.P_2_out.v) / (.000000001**2)
                else:
                    K_T2 = (self.P_prev.v - self.P_2_out.v) / (m_dot_2_new**2)

                # Compute the new fraction value: f_adj
                dPdf = ((K_T2 - K_T1)*self.f.v - K_T2) * ((K_T2-K_T1)*self.f.v**2 - 2*K_T2*self.f.v + K_T2)
                f_adj = self.f.v - Beta*dPdf

                if (self.P_1_out.v < 0) :
                    f_adj = .0000000001
                if (self.P_2_out.v < 0) :
                    f_adj = .9999999999

                f_new = LR*f_adj + (1-LR)*self.f_out.v
                m_dot_1_new = self.m_dot.v*f_new
                m_dot_2_new = self.m_dot.v*(1-f_new)


        # IF USING ESOL4050-Parallel-Flow-Solver
        elif (self.Solving_Method.v == 3.0 ):
            # Read in flow value computed by ESOL4050
            f_new = self.f.v
            #f_new_adj = f  # Added by C.Volkwein

            m_dot_1_new = f_new*self.m_dot.v
            m_dot_2_new = (1-f_new)*self.m_dot.v

        #-----------------------------------------------------------------------------------------------------------------------
        #Set the Outputs from this Model (#,Value)

        # Set output values
        self.m_dot_1.v = m_dot_1_new
        self.P_1.v = self.Pressure.v
        self.Temp_1.v = self.Temperature.v
        self.m_dot_2.v = m_dot_2_new
        self.P_2.v = self.Pressure.v
        self.Temp_2.v = self.Temperature.v
        self.f_out.v = f_new
        self.m_dot_out.v = self.m_dot.v
        self.Pressure_out.v = self.Pressure.v
        self.Mass_Counter_out.v = self.Mass_Counter.v
        return
