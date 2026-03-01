"""HTF tee return (merging) junction component model (Type 4008)."""

import numpy as np

from eeslib import fluid_properties as fp

from esclab.simulate import Component


class TeeReturnHTF(Component):
    """
    Object: ESOL4008-Tee_Return
    Simulation Studio Model: ESOL4008-Tee_Return

    Author: Matt Tuman
    Date: October 20, 2022
    last modified: October 20, 2022

    Merges two HTF return branches into a single outlet stream using a
    mass-weighted energy balance and a pressure-averaging scheme.

    Parameters
    ----------
    Diameter_1 : float
        Diameter of branch 1 pipe [m].
    Diameter_2 : float
        Diameter of branch 2 pipe [m].
    Diameter_out : float
        Diameter of outlet pipe [m].
    Fluid_ID : str or float
        Fluid identifier for specific heat property lookups.

    Inputs
    ------
    m_dot_1 : float      Branch 1 mass flow rate [kg/s].
    T_1 : float          Branch 1 temperature [K].
    P_1 : float          Branch 1 pressure [Pa].
    m_dot_2 : float      Branch 2 mass flow rate [kg/s].
    T_2 : float          Branch 2 temperature [K].
    P_2 : float          Branch 2 pressure [Pa].
    Mass_Counter_1 : float   Branch 1 mass counter [kg].
    Mass_Counter_2 : float   Branch 2 mass counter [kg].

    Outputs
    -------
    m_dot        : float  Combined outlet mass flow rate [kg/s].
    Pressure     : float  Outlet pressure [Pa].
    Temperature  : float  Outlet mixed temperature [K].
    P_1_in       : float  Echo of branch 1 inlet pressure [Pa].
    P_2_in       : float  Echo of branch 2 inlet pressure [Pa].
    Mass_Counter : float  Combined mass counter [kg].
    """

    #    PARAMETERS
    Diameter_1   = Component.Parameter()  # Diameter of branch 1 pipe [m]
    Diameter_2   = Component.Parameter()  # Diameter of branch 2 pipe [m]
    Diameter_out = Component.Parameter()  # Diameter of outlet pipe [m]
    Fluid_ID     = Component.Parameter()  # Fluid identifier for property lookups

    #    INPUTS
    m_dot_1       = Component.Input()  # Mass flow rate in branch 1 [kg/s]
    T_1           = Component.Input()  # Temperature of branch 1 [K]
    P_1           = Component.Input()  # Pressure of branch 1 [Pa]
    m_dot_2       = Component.Input()  # Mass flow rate in branch 2 [kg/s]
    T_2           = Component.Input()  # Temperature of branch 2 [K]
    P_2           = Component.Input()  # Pressure of branch 2 [Pa]
    Mass_Counter_1 = Component.Input()  # Mass counter for branch 1 [kg]
    Mass_Counter_2 = Component.Input()  # Mass counter for branch 2 [kg]

    #    OUTPUTS
    m_dot        = Component.Output()  # m_dot
    Pressure     = Component.Output()  # Pressure
    Temperature  = Component.Output()  # Temperature
    P_1_in       = Component.Output()  # Input pressure 1
    P_2_in       = Component.Output()  # Input pressure 2
    Mass_Counter = Component.Output()  # Mass Counter

    def _specific_heat(self, T, P):
        # T_1/T_2 must be within a reasonable value before calling specheat
        t_eval = T if T > 273.0 else 300.0  # set as default value if temperature not reasonable
        return float(fp.specheat(self.Fluid_ID.v, T=t_eval, P=P))

    def _compute_outlet_pressure(self, P_1, P_2):
        # find pressure leaving tee out
        if P_1 * P_2 < 0.0:    # pressures have different signs, return positive one
            return max(P_1, P_2)
        elif P_1 > 0.0:         # both pressures are positive, return average
            return (P_1 + P_2) / 2.0
        else:                   # both pressures are negative, return 1 [Pa]
            return 1.0

    def _solve_mixed_temperature(self, m_dot_tot, P_1):
        """Iterative secant-method solve for mixed outlet temperature."""
        # Mass Balance
        if np.abs(m_dot_tot) < 1.0e-12:
            return float(self.Temperature.v) if np.isfinite(self.Temperature.v) else float(
                0.5 * (self.T_1.v + self.T_2.v))

        # Energy Balance
        # compute enthalpies
        cp_T1 = self._specific_heat(self.T_1.v, self.P_1.v)
        cp_T2 = self._specific_heat(self.T_2.v, self.P_2.v)
        cp_times_T_out = (self.m_dot_1.v * cp_T1 * self.T_1.v
                          + self.m_dot_2.v * cp_T2 * self.T_2.v) / m_dot_tot

        # find correct temperature for specific heat using quadratic formula
        T_out_guess = self.Temperature.v
        tol = 100.0         # error must be less than 100 [J/kg-K]
        LR = 0.5            # learning rate
        upperbound = 1000.0  # [K] htf will not be over 1000 K
        lowerbound = 273.0   # [K]
        error = 1000.0
        error_prev = float('nan')
        T_out_guess_prev = float('nan')

        cc = 1
        while abs(error) > tol and cc < 1000:
            cc = cc + 1
            # find specific heat based on guess t_out, P_1 put in but won't be used
            cp_guess = float(fp.specheat(self.Fluid_ID.v, T=T_out_guess, P=P_1))
            cp_times_t_out_guess = T_out_guess * cp_guess
            error = cp_times_T_out - cp_times_t_out_guess  # error = actual - guess

            if abs(error) < tol:
                break  # exit the while loop
            elif cc >= 2 and T_out_guess != T_out_guess_prev:  # 2 points for error method
                if T_out_guess != T_out_guess_prev:  # will not result in divide by zero when finding the slope
                    m = (error_prev - error) / (T_out_guess_prev - T_out_guess)  # slope
                    b = error - m * T_out_guess  # y-intercept
                    if m != 0.0:
                        T_out_guess_prev = T_out_guess  # update previous points
                        error_prev = error
                        T_out_new = max(min(-b / m, upperbound), lowerbound)  # make sure new guess is within bounds
                        T_out_guess = T_out_guess + (T_out_new - T_out_guess) * LR
                    else:  # slope is equal to zero
                        if error > 0.0:  # if error is positive, increase T_out
                            T_out_guess_prev = T_out_guess  # update previous points
                            error_prev = error
                            T_out_guess = min(T_out_guess + 1.0, upperbound)  # new guess value
                        else:  # error is negative, decrease T_out
                            T_out_guess_prev = T_out_guess  # update previous points
                            error_prev = error
                            T_out_guess = max(T_out_guess - 1.0, lowerbound)  # new guess value
                else:  # last two points guessed were the same
                    T_out_guess_prev = T_out_guess
                    error_prev = error
                    if error > 0.0:
                        T_out_guess = min(T_out_guess + 1.0, upperbound)  # new guess value
                    else:
                        T_out_guess = max(T_out_guess - 1.0, lowerbound)  # new guess value
            else:  # 1st while loop iteration, guess new temperature based on sign of error
                T_out_guess_prev = T_out_guess
                error_prev = error
                if error > 0.0:
                    T_out_guess = min(T_out_guess + 1.0, upperbound)  # new guess value
                else:
                    T_out_guess = max(T_out_guess - 1.0, lowerbound)  # new guess value

        return T_out_guess

    def presim_setup(self, **kwargs):
        # Set the Initial Values of the Outputs (#,Value)
        self.m_dot.v         = self.m_dot_1.v + self.m_dot_2.v      # m_dot
        self.Pressure.v      = self.P_1.v                            # Pressure
        self.Temperature.v   = 500.0                                 # Temperature
        self.P_1_in.v        = self.P_1.v                            # Input pressure 1
        self.P_2_in.v        = self.P_2.v                            # Input pressure 2
        self.Mass_Counter.v  = self.Mass_Counter_1.v + self.Mass_Counter_2.v  # Mass Counter

    def calculate(self):
        if self.model.is_first_step:
            # Set the Initial Values of the Outputs (#,Value)
            self.m_dot.v         = self.m_dot_1.v + self.m_dot_2.v      # m_dot
            self.Pressure.v      = self.P_1.v                            # Pressure
            self.Temperature.v   = 500.0                                 # Temperature
            self.P_1_in.v        = self.P_1.v                            # Input pressure 1
            self.P_2_in.v        = self.P_2.v                            # Input pressure 2
            self.Mass_Counter.v  = self.Mass_Counter_1.v + self.Mass_Counter_2.v  # Mass Counter
            return

        # Mass Balance
        m_dot_tot = self.m_dot_1.v + self.m_dot_2.v

        # Energy Balance - solve for mixed outlet temperature
        T_out = self._solve_mixed_temperature(m_dot_tot, self.P_1.v)

        # find pressure leaving tee out
        P_out = self._compute_outlet_pressure(self.P_1.v, self.P_2.v)

        # Set the Outputs from this Model (#,Value)
        self.m_dot.v         = m_dot_tot                                  # m_dot
        self.Pressure.v      = P_out                                      # Pressure
        self.Temperature.v   = T_out                                      # Temperature
        self.P_1_in.v        = self.P_1.v                                 # Input pressure 1
        self.P_2_in.v        = self.P_2.v                                 # Input pressure 2
        self.Mass_Counter.v  = self.Mass_Counter_1.v + self.Mass_Counter_2.v  # Mass Counter
