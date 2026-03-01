"""Subcooled water pump component model (Type 6027)."""

import numpy as np

from eeslib import fluid_properties as fp

from esclab.simulate import Component
from esclab.components.flownetwork.esol6015_helpers import PB_CV_data


class SubcooledWaterPump(Component):
    """
    Object: SCWaterPumps
    Simulation Studio Model: ESOL6027-SCWaterPump

    Subcooled water pump bank with N parallel pumps and individual discharge
    valves.  A hydraulic solver finds the operating point of the lead pump
    (highest valve position) by matching pump head to total system head loss.
    Non-lead active pumps are subsequently solved to match the common outlet
    pressure established by the lead pump.

    Number of TRNSYS inputs  : 4 + 2*N
    Number of TRNSYS outputs : 9 + 7*N

    TODO-NEEDS CONVERSION REVIEW: The number of inputs (4+2*N) and outputs
    (9+7*N) depends on parameter N (number of pumps).  Fixed Component.Input()
    / Component.Output() descriptors are declared only for the N-independent
    quantities (inputs 1-4, outputs 1-6).  Per-pump inputs and outputs are
    handled via Python lists stored as instance attributes.  Model connection
    code must populate the PUMP_ON and VP_input lists manually; see
    presim_setup() for placeholder assignments.

    Parameters
    ----------
    A : float
        Pump curve coefficient A  [m / (m³/s)²]  –  head = A·Q² + B·Q + C
    B : float
        Pump curve coefficient B  [m / (m³/s)]
    C : float
        Pump curve coefficient C  [m]  (shut-off head)
    Eta_A : float
        Efficiency curve coefficient A  (4th-order polynomial in Q)
    Eta_B : float
        Efficiency curve coefficient B
    Eta_C : float
        Efficiency curve coefficient C
    Eta_D : float
        Efficiency curve coefficient D  (zero-flow efficiency intercept)
    Valve_Diameter : float
        Discharge valve inner diameter  [m]
    Valve_Speed : float
        Maximum valve actuation speed   [deg/s]
    Valve_Type : int
        Discharge valve type integer passed to PB_CV_data
    N : int
        Number of parallel pumps  (param 11)
    """

    # ------------------------------------------------------------------
    # PARAMETERS  (param number matches Fortran getParameterValue index)
    # ------------------------------------------------------------------
    A = Component.Parameter()               # param 1  – pump curve coeff A
    B = Component.Parameter()               # param 2  – pump curve coeff B
    C = Component.Parameter()               # param 3  – pump curve coeff C
    Eta_A = Component.Parameter()           # param 4  – efficiency coeff A
    Eta_B = Component.Parameter()           # param 5  – efficiency coeff B
    Eta_C = Component.Parameter()           # param 6  – efficiency coeff C
    Eta_D = Component.Parameter()           # param 7  – efficiency coeff D
    Valve_Diameter = Component.Parameter()  # param 8  – valve diameter
    Valve_Speed = Component.Parameter()     # param 9  – valve actuation speed
    Valve_Type = Component.Parameter()      # param 10 – valve type
    N = Component.Parameter()               # param 11 – number of pumps

    # ------------------------------------------------------------------
    # INPUTS  (fixed portion – indices 1-4)
    # ------------------------------------------------------------------
    P_in = Component.Input()       # input 1 – inlet pressure            [Pa]
    h_in = Component.Input()       # input 2 – inlet enthalpy            [J/kg]
    P_system = Component.Input()   # input 3 – downstream system pressure [Pa]
    P_tank2 = Component.Input()    # input 4 – tank 2 pressure           [Pa]

    # TODO-NEEDS CONVERSION REVIEW: Variable-length pump inputs cannot be
    # represented as Component.Input() class descriptors.
    #   input (4+i)   for i=1..N : PUMP_ON flag (1=on, 0=off) for pump i
    #   input (4+N+i) for i=1..N : VP_input valve-position request [0–1] for pump i
    # In the calculate() method these are read from self._PUMP_ON and
    # self._VP_input lists, which the model must populate before each step.

    # ------------------------------------------------------------------
    # OUTPUTS  (fixed summary – indices 1-6)
    # ------------------------------------------------------------------
    m_dot_total = Component.Output()   # output 1 – total mass flow out       [kg/s]
    Vol_dot_out = Component.Output()   # output 2 – total volumetric flow out  [m³/s]
    P_out = Component.Output()         # output 3 – common outlet pressure     [Pa]
    h_out = Component.Output()         # output 4 – mixed outlet enthalpy      [J/kg]
    T_out = Component.Output()         # output 5 – outlet temperature         [K]
    Work_total = Component.Output()    # output 6 – total pump shaft work      [W]

    # TODO-NEEDS CONVERSION REVIEW: Variable-length per-pump outputs (indices
    # 7+7*(i-1) through 13+7*(i-1) for i=1..N) are stored in Python lists:
    #   self._P_pump_out[i]    – output  7+7*i  Pressure leaving pump i   [Pa]
    #   self._h_pump_out[i]    – output  8+7*i  Enthalpy leaving pump i   [J/kg]
    #   self._W_dot_pump[i]    – output  9+7*i  Shaft work for pump i     [W]
    #   self._Eta_pump[i]      – output 10+7*i  Efficiency of pump i      [–]
    #   self._m_dot_pump[i]    – output 11+7*i  Mass flow through pump i  [kg/s]
    #   self._VP[i]            – output 12+7*i  Discharge valve position  [0–1]
    #   self._DELTA_P_valve[i] – output 13+7*i  Valve pressure drop       [Pa]
    # Lead-pump iteration tracking (outputs 7+7*N, 8+7*N, 9+7*N) are stored
    # in self._Point_1x, self._Point_1y, self._Point_2x.

    # ------------------------------------------------------------------
    # Physical constants
    # ------------------------------------------------------------------
    rho_water = 1000.0   # Density of Subcooled Water  [kg/m³]
    g_acc = 9.81         # Gravity                     [m/s²]

    # -----------------------------------------------------------------------------------------------------------------------
    def presim_setup(self, **kwargs):
        """
        Initialise per-pump state lists and set first-step output values.

        Corresponds to the getIsStartTime() block in the Fortran source.
        """
        N = int(self.N.v)

        # TODO-NEEDS CONVERSION REVIEW: initialise PUMP_ON and VP_input from
        # actual inputs once variable-length input connections are wired.
        self._PUMP_ON = [0.0] * N    # placeholder – model must populate
        self._VP_input = [0.0] * N   # placeholder – model must populate

        # Per-pump state lists (mirror TRNSYS stored-output arrays)
        self._P_pump_out  = [0.0] * N
        self._h_pump_out  = [0.0] * N
        self._W_dot_pump  = [0.0] * N
        self._Eta_pump    = [0.0] * N
        self._m_dot_pump  = [0.0] * N
        self._VP          = [0.0] * N
        self._DELTA_P_valve = [0.0] * N

        # Lead-pump secant-iteration trackers (stored as outputs in Fortran)
        self._Point_1x = 0.0
        self._Point_1y = 0.0
        self._Point_2x = 0.0   # holds Q_new from the first iteration

        Flow = 1.0

        m_dot_total_init = 0.0
        h_out_sum = 0.0

        for i in range(N):
            # Set Valve Position leaving the pump
            VP_input = self._VP_input[i]
            self._VP[i] = VP_input

            # Solve for Pressure leaving the pump
            PUMP_ON = self._PUMP_ON[i]
            if PUMP_ON == 1:
                m_dot_pump = Flow * self.rho_water
                P_pump_out = self.C.v * self.rho_water * self.g_acc
                T_pump_out = 300.0
                h_pump_out = fp.enthalpy("water", T=T_pump_out, P=P_pump_out / 1000.0) * 1000.0
            else:
                m_dot_pump = 0.0
                P_pump_out = self.P_in.v
                T_pump_out = 300.0
                h_pump_out = fp.enthalpy("water", T=T_pump_out, P=P_pump_out / 1000.0) * 1000.0

            m_dot_total_init += m_dot_pump
            h_out_sum += m_dot_pump * h_pump_out

            self._P_pump_out[i] = P_pump_out
            self._h_pump_out[i] = h_pump_out
            self._m_dot_pump[i] = m_dot_pump

        if m_dot_total_init > 0.0:
            self.m_dot_total.v = m_dot_total_init
            self.Vol_dot_out.v = m_dot_total_init / self.rho_water
            # C * rho_water * g gives the initial outlet pressure estimate
            self.P_out.v = self.C.v * self.rho_water * self.g_acc
            self.h_out.v = h_out_sum / m_dot_total_init
        else:
            self.m_dot_total.v = 0.0000000001
            self.Vol_dot_out.v = 0.0
            self.P_out.v = self.P_in.v
            self.h_out.v = self.h_in.v

    # -----------------------------------------------------------------------------------------------------------------------
    def calculate(self):
        """
        Main iterative calculation executed at every simulation iteration.

        Corresponds to the main calculation block in the Fortran source
        (after all TRNSYS housekeeping guards).
        """
        N = int(self.N.v)
        tol = 0.01
        Q_min = 1e-6
        # Maximum volumetric flow from pump curve (positive root of A*Q^2+B*Q+C=0)
        discriminant = self.B.v**2 - 4.0 * self.A.v * self.C.v
        Q_max = max(
            (-self.B.v + np.sqrt(discriminant)) / (2.0 * self.A.v),
            (-self.B.v - np.sqrt(discriminant)) / (2.0 * self.A.v),
        )

        s_pump_in = fp.entropy("water", P=self.P_in.v / 1000.0, h=self.h_in.v / 1000.0) * 1000.0

        # TODO-NEEDS CONVERSION REVIEW: read variable-length inputs into lists each
        # call.  Replace placeholder accesses below with actual input connections.
        PUMP_ON = self._PUMP_ON    # list of N on/off flags (input indices 5..4+N)
        VP_input = self._VP_input  # list of N position requests (input indices 4+N+1..4+2*N)

        # When first iteration of timestep: update valve positions if input does not match output
        if self.model.timestep_iteration == 0:
            for i in range(N):
                vp_req = VP_input[i]  # User Input Valve Position
                vp_cur = self._VP[i]
                if vp_cur != vp_req:
                    ts = self.model.timestep * 3600  # Convert timestep from hr to seconds
                    if vp_cur < vp_req:  # increase valve position based on timestep
                        VP_d = vp_cur * 90.0           # Convert from % to degrees
                        VP_input_d = vp_req * 90.0
                        VP_d = min(VP_d + self.Valve_Speed.v * ts, VP_input_d)  # Update Valve Position
                        self._VP[i] = VP_d / 90.0      # Convert back to %
                    else:  # Decrease Valve Position based on Timestep
                        VP_d = vp_cur * 90.0
                        VP_input_d = vp_req * 90.0
                        VP_d = max(VP_d - self.Valve_Speed.v * ts, VP_input_d)
                        self._VP[i] = VP_d / 90.0      # Update Valve Position

        # Find the lead pump: the one that is experiencing the lowest system
        # losses based on the valve (active pump with the highest valve position)
        N_lead = -1   # 0-based index; -1 means none found
        VP_high = 0.0
        for i in range(N):
            if PUMP_ON[i] == 1.0:
                if self._VP[i] > VP_high:
                    VP_high = self._VP[i]
                    N_lead = i

        if N_lead != -1:  # At least one pump is on

            # ------------------------------------------------------------------
            # SOLVE THE LEAD PUMP FIRST
            # ------------------------------------------------------------------
            if self.model.timestep_iteration == 0:
                Q_prev = self._m_dot_pump[N_lead] / self.rho_water
                P_pump_prev = self._P_pump_out[N_lead]
                P_out_prev = self.P_out.v

                pump_head = P_pump_prev / self.rho_water / self.g_acc - self.P_in.v / self.rho_water / self.g_acc
                DELTA_P_sys = (P_out_prev - self.P_system.v) / self.rho_water / self.g_acc
                DELTA_P_valve = (P_pump_prev - P_out_prev) / self.rho_water / self.g_acc
                Point_1y = pump_head - DELTA_P_sys - DELTA_P_valve
                Point_1x = Q_prev

                if abs(Point_1y) <= tol:
                    Q_new = max(min(Point_1x, Q_max), Q_min)
                elif Point_1y > 0.0:   # increase flow rate
                    Q_new = min(Point_1x + 0.01, Q_max)  # Move to maximum flow pump can provide at given speed
                elif Point_1y < 0.0:   # decrease flow rate
                    Q_new = max(Point_1x - 0.01, Q_min)  # Move to minimum flow pump can provide

                s_pump_in = fp.entropy("water", P=self.P_in.v / 1000.0, h=self.h_in.v / 1000.0) * 1000.0
                m_dot_pump = Q_new * self.rho_water
                P_pump_out = (self.A.v * Q_new**2 + self.B.v * Q_new + self.C.v) * self.rho_water * 9.81 + self.P_in.v

                Eta_Pump = max(
                    self.Eta_A.v * Q_new**4
                    + self.Eta_B.v * Q_new**3
                    + self.Eta_C.v * Q_new * 2.0
                    + self.Eta_D.v,
                    0.01,
                )

                s_pump_out_s = s_pump_in
                h_pump_out_s = fp.enthalpy("water", P=P_pump_out / 1000.0, s=s_pump_out_s / 1000.0) * 1000.0
                W_dot_pump_s = m_dot_pump * (h_pump_out_s - self.h_in.v)
                W_dot_pump = W_dot_pump_s / Eta_Pump
                h_pump_out = (m_dot_pump * self.h_in.v + W_dot_pump) / m_dot_pump

                self._P_pump_out[N_lead] = P_pump_out   # Pressure out of pump [Pa]
                self._h_pump_out[N_lead] = h_pump_out   # Enthalpy out of pump [J/kg]
                self._W_dot_pump[N_lead] = W_dot_pump   # Work needed for pump [W]
                self._Eta_pump[N_lead]   = Eta_Pump     # Efficiency of pump [% base 1]
                self._m_dot_pump[N_lead] = m_dot_pump   # mass flow through pump [kg/s]
                self._Point_1x = Point_1x
                self._Point_1y = Point_1y
                self._Point_2x = Q_new

            else:  # second iteration of timestep or more
                Point_1x = self._Point_1x
                Point_1y = self._Point_1y
                Point_2x = self._Point_2x
                P_pump_prev = self._P_pump_out[N_lead]
                P_out_prev = self.P_out.v

                # Compute system head loss
                pump_head = (P_pump_prev - self.P_in.v) / self.rho_water / self.g_acc
                DELTA_P_sys = (P_out_prev - self.P_system.v) / self.rho_water / self.g_acc
                DELTA_P_valve = (P_pump_prev - P_out_prev) / self.rho_water / self.g_acc
                Point_2y = pump_head - DELTA_P_sys - DELTA_P_valve

                if Point_2x == Point_1x:
                    m = 0.0
                else:
                    m = (Point_2y - Point_1y) / (Point_2x - Point_1x)
                y_int = Point_1y - m * Point_1x

                if abs(Point_2y) <= tol:
                    Q_new = Point_2x
                elif m != 0.0:
                    if Point_2y >= 0.0:
                        Q_new = max(-y_int / m, Point_2x)   # Ensure new Q_dot does not decrease
                        Q_new = min(Q_new, Q_max)            # Ensure Q_dot_new is not more than pump can provide based on curve
                    else:
                        Q_new = min(-y_int / m, Point_2x)   # Ensure new Q_dot does not increase
                        Q_new = max(Q_new, Q_min)            # Ensure new Q_dot is not negative
                    Q_new = Point_2x + (Q_new - Point_2x) * 0.4
                else:  # If Error is Negative, pump head guess is too low so flow must decrease
                    if Point_2y >= 0.0:    # Error is Positive and slope is equal to zero, pump head is too high so flow must increase
                        Q_new = min(Point_2x + 0.1, Q_max)
                    else:                  # Error is Negative and slope is equal to zero, pump head is too low so flow must decrease
                        Q_new = max(Point_2x - 0.1, Q_min)

                s_pump_in = fp.entropy("water", P=self.P_in.v / 1000.0, h=self.h_in.v / 1000.0) * 1000.0
                m_dot_pump = Q_new * self.rho_water
                P_pump_out = (self.A.v * Q_new**2 + self.B.v * Q_new + self.C.v) * self.rho_water * 9.81 + self.P_in.v

                Eta_Pump = max(
                    self.Eta_A.v * Q_new**4
                    + self.Eta_B.v * Q_new**3
                    + self.Eta_C.v * Q_new * 2.0
                    + self.Eta_D.v,
                    0.01,
                )

                s_pump_out_s = s_pump_in
                h_pump_out_s = fp.enthalpy("water", P=P_pump_out / 1000.0, s=s_pump_out_s / 1000.0) * 1000.0
                W_dot_pump_s = m_dot_pump * (h_pump_out_s - self.h_in.v)
                W_dot_pump = W_dot_pump_s / Eta_Pump
                h_pump_out = (m_dot_pump * self.h_in.v + W_dot_pump) / m_dot_pump

                Point_1x = Point_2x
                Point_1y = Point_2y
                Point_2x = Q_new

                self._P_pump_out[N_lead] = P_pump_out   # Pressure leaving pump [Pa]
                self._h_pump_out[N_lead] = h_pump_out   # Enthalpy leaving pump [J/kg]
                self._W_dot_pump[N_lead] = W_dot_pump   # Work needed for pump [W]
                self._Eta_pump[N_lead]   = Eta_Pump     # Efficiency of pump [% base 1]
                self._m_dot_pump[N_lead] = m_dot_pump   # Flow rate through pump [kg/s]
                self._Point_1x = Point_1x
                self._Point_1y = Point_1y
                self._Point_2x = Point_2x

            # Solve for the lead valve pressure drop and outlet pressure
            VP = self._VP[N_lead]
            CV = PB_CV_data(int(self.Valve_Type.v), self.Valve_Diameter.v, VP)
            Vol_dot_gpm = self._m_dot_pump[N_lead] / self.rho_water * 15850.323140625002
            DELTA_P_Valve = Vol_dot_gpm**2 / CV**2 * 6894.76

            P_out = self._P_pump_out[N_lead] - DELTA_P_Valve

            # Solve for mass flow rate leaving other pumps based on P_out
            for i in range(N):
                if i != N_lead:
                    if PUMP_ON[i] == 1:
                        P_pump_out = self._P_pump_out[i]
                        Q_pump_guess = self._m_dot_pump[i] / self.rho_water
                        VP = self._VP[i]
                        CV = PB_CV_data(int(self.Valve_Type.v), self.Valve_Diameter.v, VP)
                        error = tol + 1.0
                        whileiterations = 0
                        error_prev = 0.0
                        Q_pump_prev = Q_pump_guess
                        while abs(error) > tol:
                            whileiterations += 1
                            Vol_dot_gpm = Q_pump_guess * 15850.323140625002
                            DELTA_P_Valve = Vol_dot_gpm**2 / CV**2 * 6894.76
                            error = (P_pump_out - DELTA_P_Valve - P_out) / self.rho_water / self.g_acc
                            if Q_pump_guess == Q_min:
                                if error < 0.0:
                                    error = tol / 2.0
                                    Q_pump_guess = Q_min
                            if Q_pump_guess == Q_max:
                                if error > 0.0:
                                    error = tol / 2.0
                                    Q_pump_guess = Q_max
                            if abs(error) < tol:
                                m_dot_pump = Q_pump_guess * self.rho_water
                                P_pump_out = (
                                    self.A.v * Q_pump_guess**2
                                    + self.B.v * Q_pump_guess
                                    + self.C.v
                                ) * self.rho_water * self.g_acc + self.P_in.v
                            elif whileiterations == 1:
                                error_prev = error
                                Q_pump_prev = Q_pump_guess
                                if error > 0.0:  # Increase Flow Rate
                                    Q_pump_guess = min(Q_pump_guess + 0.01, Q_max)
                                else:
                                    Q_pump_guess = max(Q_pump_guess - 0.01, Q_min)
                                P_pump_out = (
                                    self.A.v * Q_pump_guess**2
                                    + self.B.v * Q_pump_guess
                                    + self.C.v
                                ) * self.rho_water * self.g_acc + self.P_in.v
                            else:
                                if Q_pump_prev != Q_pump_guess:
                                    m_slope = (error - error_prev) / (Q_pump_guess - Q_pump_prev)
                                    y_int_inner = error - m_slope * Q_pump_guess
                                    Q_pump_prev = Q_pump_guess
                                    error_prev = error
                                    if m_slope != 0.0:
                                        if error > 0.0:
                                            Q_new_inner = min(-y_int_inner / m_slope, Q_max)
                                        else:
                                            Q_new_inner = max(-y_int_inner / m_slope, Q_min)
                                        Q_pump_guess = Q_pump_guess + (Q_new_inner - Q_pump_guess) * 0.4
                                    else:
                                        if error > 0.0:
                                            Q_pump_guess = min(Q_pump_guess + 0.01, Q_max)
                                        else:
                                            Q_pump_guess = max(Q_pump_guess - 0.01, Q_min)
                                else:
                                    if error > 0.0:
                                        Q_pump_guess = min(Q_pump_guess + 0.01, Q_max)
                                    else:
                                        Q_pump_guess = max(Q_pump_guess - 0.01, Q_min)
                                P_pump_out = (
                                    self.A.v * Q_pump_guess**2
                                    + self.B.v * Q_pump_guess
                                    + self.C.v
                                ) * self.rho_water * self.g_acc + self.P_in.v

                        # solve for enthalpy leaving the pump
                        Eta_pump_i = max(
                            self.Eta_A.v * Q_pump_guess**4
                            + self.Eta_B.v * Q_pump_guess**3
                            + self.Eta_C.v * Q_pump_guess**2
                            + self.Eta_D.v * Q_pump_guess,
                            0.01,
                        )
                        s_pump_out_s = s_pump_in
                        h_pump_out_s = fp.enthalpy("water", P=P_pump_out / 1000.0, s=s_pump_out_s / 1000.0) * 1000.0
                        W_dot_pump_s = m_dot_pump * (h_pump_out_s - self.h_in.v)
                        W_dot_pump = W_dot_pump_s / Eta_pump_i
                        h_pump_out = (m_dot_pump * self.h_in.v + W_dot_pump) / m_dot_pump

                        self._P_pump_out[i]    = P_pump_out       # Pressure leaving Pump [Pa]
                        self._h_pump_out[i]    = h_pump_out       # Enthalpy leaving Pump [J/kg]
                        self._W_dot_pump[i]    = W_dot_pump       # Work needed for Pump [W]
                        self._Eta_pump[i]      = Eta_pump_i       # Pump Efficiency [% base 1]
                        self._m_dot_pump[i]    = m_dot_pump       # Flow rate through pump [kg/s]
                        self._DELTA_P_valve[i] = DELTA_P_Valve    # pressure drop across valve [Pa]

                    else:  # PUMP IS NOT ON
                        self._P_pump_out[i]  = self.P_in.v
                        self._h_pump_out[i]  = self.h_in.v
                        self._W_dot_pump[i]  = 0.0
                        self._Eta_pump[i]    = 0.0
                        self._m_dot_pump[i]  = 0.0

            # solve for mixing the pump flows together
            Work_total_sum = 0.0
            h_out_sum = 0.0
            m_dot_total_sum = 0.0
            for i in range(N):
                Work_total_sum += self._W_dot_pump[i]
                m_dot_pump_i = self._m_dot_pump[i]
                m_dot_total_sum += m_dot_pump_i
                h_out_sum += self._h_pump_out[i] * m_dot_pump_i

            h_out_mixed = h_out_sum / m_dot_total_sum
            T_out_val = fp.temperature("water", P=P_out / 1000.0, h=h_out_mixed / 1000.0)
            Vol_dot_out_val = m_dot_total_sum / self.rho_water

            self.m_dot_total.v = m_dot_total_sum
            self.Vol_dot_out.v = Vol_dot_out_val
            self.P_out.v = P_out
            self.h_out.v = h_out_mixed
            self.T_out.v = T_out_val
            self.Work_total.v = Work_total_sum

        else:  # All pumps are off

            T_in = fp.temperature("water", P=self.P_in.v / 1000.0, h=self.h_in.v / 1000.0)
            self.m_dot_total.v = 0.00001
            self.Vol_dot_out.v = 0.0
            self.P_out.v = self.P_in.v
            self.h_out.v = self.h_in.v
            self.T_out.v = T_in
            self.Work_total.v = 0.0

            for i in range(N):
                self._P_pump_out[i]    = self.P_in.v
                self._h_pump_out[i]    = self.h_in.v
                self._W_dot_pump[i]    = 0.0
                self._Eta_pump[i]      = 0.0
                self._m_dot_pump[i]    = 0.0
                self._DELTA_P_valve[i] = 0.0
