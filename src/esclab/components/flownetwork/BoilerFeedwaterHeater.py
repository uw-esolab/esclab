"""Boiler feedwater heater component model (Type 6019)."""

from eeslib import fluid_properties as fp

from esclab.simulate import Component

class BoilerFeedwaterHeater(Component):
    """
    TRNSYS Type 6019: ESOL6019-BFWH.

    Boiler Feedwater Heater Type. Heats feedwater using turbine stage extraction
    steam and/or drain bleed from a downstream feedwater heater. Solves for the
    off-design feedwater temperature increase based on the steam-drum pressure
    ratio relative to design, then determines the required steam extraction from
    the turbine stage while respecting maximum extraction limits.

    Parameters
    ----------
    DELTA_T_design : float
        Design feedwater temperature increase [K].
    TTD_BFWH : float
        Terminal Temperature Difference for BFWH; (T_sat_extraction - T_fw_out) = TTD [K].
    eta_sc : float
        Subcooling effectiveness [-].
    P_SD_design : float
        Design Steam Drum pressure; used to calculate off-design DELTA_T [Pa].
    perc_min : float
        Minimum turbine stage pressure as a fraction of design where extractions
        can still be drawn from the turbine [%].
    percent_max_extraction : float
        Fraction of turbine stage flow that can be extracted at most [-].
    flow_tol : float
        Tolerance that controls when the type solves for a new off-design
        temperature increase based on how close the turbines are to convergence
        (turbine error) [kg/s].
    OD_Coef : float
        Scaling coefficient for off-design DELTA_T:
        DELTA_T_OD / DELTA_T_design = OD_Coef * log(P_SD / P_SD_design) + 1 [-].
    DELTA_T_tol : float
        Tolerance for changing the turbine stage temperature [K].
    m_dot_b_guess : float
        Guess steam extraction to request from turbine during first iteration [kg/s].
    DELTA_T_guess : float
        Guess feedwater temperature increase during first iteration [K].

    Inputs
    ------
    Turbine_On : float
        Turbine running signal (1 == On, 0 == Off) [-].
    m_dot_fw : float
        Feedwater mass flow rate entering the BFWH [kg/s].
    P_fw_in : float
        Feedwater inlet pressure [Pa].
    h_fw_in : float
        Feedwater inlet enthalpy [J/kg].
    m_dot_turbine_stage_max : float
        Maximum flow entering the turbine stage; used to find maximum extraction [kg/s].
    P_turbine_stage : float
        Pressure of the turbine stage [Pa].
    h_turbine_stage : float
        Enthalpy of the turbine stage [J/kg].
    m_dot_drain : float
        Additional bleed entering from the upstream feedwater heater [kg/s].
    h_drain_in : float
        Enthalpy of bleed entering from the upstream feedwater heater [J/kg].
    P_SD : float
        Current steam drum pressure; used to find off-design DELTA_T [Pa].

    Outputs
    -------
    m_dot_fw_out : float
        Mass flow rate of feedwater leaving the BFWH [kg/s].
    Vol_dot_fw : float
        Volumetric flow rate of feedwater leaving the BFWH [m³/s].
    P_fw_out : float
        Pressure of feedwater leaving the BFWH [Pa].
    h_fw_out : float
        Enthalpy of feedwater leaving the BFWH [J/kg].
    T_fw_out : float
        Temperature of feedwater leaving the BFWH [K].
    m_dot_total : float
        Total steam extraction mass flow leaving the feedwater heater
        (drain from previous + extraction from turbine combined) [kg/s].
    Vol_dot_total : float
        Total steam extraction volumetric flow leaving the feedwater heater
        (drain from previous + extraction from turbine combined) [m³/s].
    P_shell : float
        Pressure of the extraction side leaving the feedwater heater [Pa].
    h_b_out : float
        Enthalpy of the extraction side leaving the feedwater heater [J/kg].
    T_b_out : float
        Temperature of the extraction leaving the feedwater heater [K].
    Q_dot_act : float
        Total heat transfer from the turbine extraction to the feedwater [W].
    m_dot_b : float
        Requested turbine bleed from BFWH type [kg/s].
    DELTA_T_OD : float
        Off-design feedwater temperature increase wanted by the BFWH
        (not always achieved) [K].
    DELTA_T_act : float
        Actual feedwater temperature increase [K].
    """

    # *** Model Parameters ***
    # DELTA_T_design: design feedwater temperature increase
    DELTA_T_design = Component.Parameter()
    # TTD_BFWH: Terminal Temperature Difference for BFWH, (T_sat_extraction - T_fw_out) = TTD
    TTD_BFWH = Component.Parameter()
    # eta_sc: subcooling effectiveness
    eta_sc = Component.Parameter()
    # P_SD_design: Design SD Pressure, used to calculate off-design DELTA_T
    P_SD_design = Component.Parameter()
    # perc_min: finds minimum turbine stage pressures where extractions can still be drawn from turbine [%]
    perc_min = Component.Parameter()
    # percent_max_extraction: Percent_max_extraction, used to limit extraction amount from turbine if necessary
    percent_max_extraction = Component.Parameter()
    # flow_tol: Tolerance that lets type solve for new off-design temperature increase
    #           based on how close the turbines are to convergence (turbine error)
    flow_tol = Component.Parameter()
    # OD_Coef: Scaling Coefficient for OD DELTA T
    #          --> DELTA_T_OD/DELTA_T_Design = Coef * log(P_turbine_stage_OD/P_turbine_stage_design) + 1
    OD_Coef = Component.Parameter()
    # DELTA_T_tol: Tolerance for changing the turbine stage temperature
    DELTA_T_tol = Component.Parameter()
    # m_dot_b_guess: Guess steam extraction to request from turbine during first iteration
    m_dot_b_guess = Component.Parameter()
    # DELTA_T_guess: Guess Feedwater Temperature Increase during first iteration
    DELTA_T_guess = Component.Parameter()
    # FWH_sequence_number: Parameter to identify the sequence from low to high pressure of the FWH's. 
    FWH_sequence_number = Component.Parameter()

    # *** Model Inputs ***
    Turbine_On = Component.Input()
    m_dot_fw = Component.Input()
    P_fw_in = Component.Input()
    h_fw_in = Component.Input()
    m_dot_turbine_stage_max = Component.Input()
    P_turbine_stage = Component.Input()
    h_turbine_stage = Component.Input()
    m_dot_drain = Component.Input()
    h_drain_in = Component.Input()
    P_SD = Component.Input()

    # *** Model Outputs ***
    m_dot_fw_out = Component.Output()
    Vol_dot_fw = Component.Output()
    P_fw_out = Component.Output()
    h_fw_out = Component.Output()
    T_fw_out = Component.Output()
    m_dot_total = Component.Output()
    Vol_dot_total = Component.Output()
    P_shell = Component.Output()
    h_b_out = Component.Output()
    T_b_out = Component.Output()
    Q_dot_act = Component.Output()
    m_dot_b = Component.Output()
    DELTA_T_OD = Component.Output()
    DELTA_T_act = Component.Output()

    def calculate(self):

        # Do All of the First Timestep Manipulations Here - There Are No Iterations at the Initial Time
        if self.model.is_first_step:

            # CONVERTED-NEEDS UNITS CHECK: removed unit scale factors; eeslib uses SI
            T_fw_in = fp.temperature("water", P=self.P_fw_in.v, h=self.h_fw_in.v)
            rho_fw = fp.density("water", P=self.P_fw_in.v, h=self.h_fw_in.v)

            if rho_fw > 0.0:
                Vol_dot_fw = self.m_dot_fw.v / rho_fw
            else:
                Vol_dot_fw = self.m_dot_fw.v / 1000.0

            # Check the Parameters for Problems (#,ErrorType,Text)
            # Sample Code: If( PAR1 <= 0.) Call FoundBadParameter(1,'Fatal','The first parameter provided to this model is not acceptable.')

            # Set the Initial Values of the Outputs (#,Value)
            self.m_dot_fw_out.v = self.m_dot_fw.v       # Mass flow rate of feedwater leaving the boiler feedwater heater
            self.Vol_dot_fw.v = Vol_dot_fw               # Volumetric flow rate of feedwater leaving the boiler feedwater heater
            self.P_fw_out.v = self.P_fw_in.v             # Pressure of feedwater leaving
            self.h_fw_out.v = self.h_fw_in.v             # Enthalpy of feedwater leaving
            self.T_fw_out.v = T_fw_in                    # Temperature of Feedwater leaving
            self.m_dot_total.v = 0.0                     # Mass flow rate of steam extraction draining from the boiler feedwater heater
            self.Vol_dot_total.v = 0.0                   # Volumetric flow rate of steam extraction draining from boiler feedwater heater
            self.P_shell.v = 300000.0                    # Pressure of steam side leaving
            self.h_b_out.v = 0.0                         # Enthalpy of steam side leaving
            self.T_b_out.v = 0.0                         # Temperature of steam side leaving
            self.Q_dot_act.v = 0.0                       # Heat Transferred from Steam Side to BFWHs
            self.m_dot_b.v = self.m_dot_b_guess.v        # Steam Flow Requested from Turbine Stage
            self.DELTA_T_OD.v = self.DELTA_T_guess.v     # Off-design feedwater temperature increase through boiler feedwater heater
            self.DELTA_T_act.v = self.DELTA_T_guess.v

            return

        # ReRead the Parameters if Another Unit of This Type Has Been Called Last
        # (getIsReReadParameters block – no parameter re-reading required here)

        # Read the Inputs (parameters and inputs accessed directly via self members per step 5)

        # CHECK IF TURBINE IS ON
        if self.Turbine_On.v == 1:
            if self.m_dot_fw.v >= 0.1:
                # SOLVE FOR FEEDWATER TEMPERATURE INCREASE THROUGH BFWH BASED ON TURBINE STAGE PRESSURE
                DELTA_T_old = self.DELTA_T_OD.v  # getOutputValue(13)

                # Solve for new off-design FW temp Increase
                P_min = self.P_SD_design.v * self.perc_min.v  # bleed will be zero if turbine stage pressure is under this pressure
                P_ratio = self.P_SD.v / self.P_SD_design.v
                if self.P_SD.v <= P_min:
                    DELTA_T_new = 0.0
                elif self.P_SD.v >= self.P_SD_design.v:
                    DELTA_T_new = self.DELTA_T_design.v
                else:
                    # (MJW) The source fortran checks the the FWH correction curve based on deck file unit number,
                    # which is almost guaranteed to have caused problems. The most recent deck file 'predawn_solana' 
                    # has unit numbers that differ from the expected range. 
                    # The correlations were verified as being mapped to 5-3 using the BFWH Ramping.xlsx sheet in 
                    # Anna's folder.
                    if self.FWH_sequence_number.v == 5:  #BFWH5, lowest pressure 
                        T_ratio = max(min(-10.752*P_ratio**4.+30.139*P_ratio**3. -31.579*P_ratio**2. +15.49*P_ratio-2.2949, 1.),0.)
                    elif self.FWH_sequence_number.v == 4: #BFWH4, second lowest pressure
                        T_ratio = -18.391*P_ratio**4. +50.51*P_ratio**3.-50.34*P_ratio**2.+21.902*P_ratio-2.7355
                    elif self.FWH_sequence_number.v == 3: #BFWH3, third lowest pressure
                        T_ratio = -8.1877*P_ratio**4. +23.805*P_ratio**3.-24.809*P_ratio**2.+11.266*P_ratio-1.0745
                    else:  #other BFWHs
                        T_ratio = self.OD_Coef.v * __import__("math").log(P_ratio) + 1.0
                        DELTA_T_new = T_ratio * self.DELTA_T_design.v

                LR = 0.2
                if abs(DELTA_T_new - DELTA_T_old) <= self.DELTA_T_tol.v:
                    # IF NEW DELTA_T AND OLD DELTA_T ARE WITHIN XX OF EACH OTHER USE OLD DELTA_T TO HELP CONVERGENCE
                    DELTA_T_OD = DELTA_T_old
                else:
                    DELTA_T_OD = DELTA_T_old + (DELTA_T_new - DELTA_T_old) * LR

                DELTA_T_act = DELTA_T_OD

                # CONVERTED-NEEDS UNITS CHECK: removed unit scale factors; eeslib uses SI
                T_fw_in = fp.temperature("water", P=self.P_fw_in.v, h=self.h_fw_in.v)
                # CONVERTED-NEEDS UNITS CHECK: removed unit scale factors; eeslib uses SI
                T_sat_fw = fp.temperature("water", P=self.P_fw_in.v, quality=0.0)

                if T_sat_fw < T_fw_in + DELTA_T_OD:
                    # Decrease DELTA_T to avoid entering the saturation region
                    DELTA_T_act = T_sat_fw - T_fw_in

                T_fw_out = T_fw_in + DELTA_T_act
                T_sat_b = T_fw_out + self.TTD_BFWH.v
                # Solve for shell pressure based on saturated temperature of steam extraction
                # CONVERTED-NEEDS UNITS CHECK: removed unit scale factors; eeslib uses SI
                P_shell = fp.pressure("water", T=T_sat_b, quality=1.0)
                # CONVERTED-NEEDS UNITS CHECK: removed unit scale factors; eeslib uses SI
                h_b_sat_g = fp.enthalpy("water", T=T_sat_b, quality=1.0)

                if T_fw_out < T_sat_fw:
                    # CONVERTED-NEEDS UNITS CHECK: removed unit scale factors; eeslib uses SI
                    h_fw_out = fp.enthalpy("water", T=T_fw_out, P=self.P_fw_in.v)
                else:
                    # CONVERTED-NEEDS UNITS CHECK: removed unit scale factors; eeslib uses SI
                    h_fw_out = fp.enthalpy("water", P=self.P_fw_in.v, quality=0.0)
                Q_dot_fw = self.m_dot_fw.v * (h_fw_out - self.h_fw_in.v)  # total amount of heat transfer needed by the feedwater

                # Solve for heat from drain to fw
                # CONVERTED-NEEDS UNITS CHECK: removed unit scale factors; eeslib uses SI
                h_b_sat_f = fp.enthalpy("water", P=P_shell, quality=0.0)
                # CONVERTED-NEEDS UNITS CHECK: removed unit scale factors; eeslib uses SI
                T_sat_b = fp.temperature("water", P=P_shell, quality=0.0)

                if abs(T_sat_b - T_fw_in) >= 0.1:
                    # CONVERTED-NEEDS UNITS CHECK: removed unit scale factors; eeslib uses SI
                    h_b_min = fp.enthalpy("water", P=P_shell, T=T_fw_in)
                else:
                    # CONVERTED-NEEDS UNITS CHECK: removed unit scale factors; eeslib uses SI
                    h_b_min = fp.enthalpy("water", P=P_shell, T=T_fw_in - 0.1)

                # Find heat added from draining steam extraction to FW
                Q_dot_drain = (
                    self.m_dot_drain.v * (max(self.h_drain_in.v, h_b_sat_f) - h_b_sat_f)
                    + self.eta_sc.v * (self.m_dot_drain.v * (min(self.h_drain_in.v, h_b_sat_f) - h_b_min))
                )

                if Q_dot_fw >= Q_dot_drain:
                    # Calculate the amount of steam extraction from turbine needed
                    Q_dot_b = Q_dot_fw - Q_dot_drain
                    m_dot_b = Q_dot_b / (
                        (max(self.h_turbine_stage.v, h_b_sat_f) - h_b_sat_f)
                        + self.eta_sc.v * (min(self.h_turbine_stage.v, h_b_sat_f) - h_b_min)
                    )
                    check = (
                        m_dot_b * (max(self.h_turbine_stage.v, h_b_sat_f) - h_b_sat_f)
                        + m_dot_b * self.eta_sc.v * (min(self.h_turbine_stage.v, h_b_sat_f) - h_b_min)
                    )
                    m_dot_b_prev = self.m_dot_b.v  # getOutputValue(12)
                    m_dot_max = self.m_dot_turbine_stage_max.v * self.percent_max_extraction.v
                    # make sure it's not drawing more than the turbines can handle
                    if m_dot_b > m_dot_max:
                        m_dot_b = m_dot_max
                    Q_dot_act = Q_dot_fw
                    DELTA_T_act = DELTA_T_OD

                    m_dot_b_prev = self.m_dot_b.v  # getOutputValue(12)
                    if abs(m_dot_b - m_dot_b_prev) <= self.flow_tol.v:
                        # keep flow the same to bring to convergence
                        m_dot_b = m_dot_b_prev
                        Q_dot_b = (
                            m_dot_b * (max(self.h_turbine_stage.v, h_b_sat_f) - h_b_sat_f)
                            + m_dot_b * self.eta_sc.v * (min(self.h_turbine_stage.v, h_b_sat_f) - h_b_min)
                        )
                        Q_dot_act = Q_dot_drain + Q_dot_b
                        h_fw_out = (self.h_fw_in.v * self.m_dot_fw.v + Q_dot_act) / self.m_dot_fw.v
                        # CONVERTED-NEEDS UNITS CHECK: removed unit scale factors; eeslib uses SI
                        T_fw_out = fp.temperature("water", P=self.P_fw_in.v, h=h_fw_out)
                        DELTA_T_act = T_fw_out - T_fw_in

                else:
                    # No steam extraction needed, calculate the actual temperature increase
                    # of the feedwater from the steam draining
                    m_dot_b = 0.0
                    Q_dot_b = 0.0
                    Q_dot_act = Q_dot_drain
                    h_fw_out = (self.h_fw_in.v * self.m_dot_fw.v + Q_dot_drain) / self.m_dot_fw.v
                    # CONVERTED-NEEDS UNITS CHECK: removed unit scale factors; eeslib uses SI
                    T_fw_out = fp.temperature("water", P=self.P_fw_in.v, h=h_fw_out)
                    DELTA_T_act = T_fw_out - T_fw_in

                m_dot_total = self.m_dot_drain.v + m_dot_b
                if m_dot_total > 0.0:
                    h_b_out = (
                        m_dot_b * self.h_turbine_stage.v
                        + self.m_dot_drain.v * self.h_drain_in.v
                        - Q_dot_act
                    ) / m_dot_total
                    # CONVERTED-NEEDS UNITS CHECK: removed unit scale factors; eeslib uses SI
                    T_b_out = fp.temperature("water", P=P_shell, h=h_b_out)
                    rho_b = fp.density("water", P=P_shell, h=h_b_out)
                    if rho_b > 0.0:
                        Vol_dot_total = m_dot_total / rho_b
                    else:
                        Vol_dot_total = m_dot_total / 1000.0
                else:
                    h_b_out = 0.0
                    T_b_out = 0.0
                    Vol_dot_total = 0.0

            else:  # no fw entering the BFWH Type
                h_fw_out = self.h_fw_in.v
                Q_dot_act = 0.0
                m_dot_total = self.m_dot_drain.v
                m_dot_b = 0.0
                Vol_dot_total = 0.0
                h_b_out = self.h_drain_in.v
                T_fw_out = 0.0
                DELTA_T_OD = self.DELTA_T_OD.v
                DELTA_T_act = self.DELTA_T_act.v
                P_shell = self.P_shell.v
                T_b_out = self.T_b_out.v

        else:
            h_fw_out = self.h_fw_in.v  # No heat transfer when turbines are off
            Q_dot_act = 0.0
            m_dot_total = 0.0
            m_dot_b = 0.0
            Vol_dot_total = 0.0
            P_shell = 0.0
            h_b_out = 0.0
            T_b_out = 0.0
            T_fw_out = 0.0
            DELTA_T_OD = self.DELTA_T_OD.v
            DELTA_T_act = self.DELTA_T_act.v

        # Calculating Feedwater Volumetric Flow Rate
        # CONVERTED-NEEDS UNITS CHECK: removed unit scale factors; eeslib uses SI
        rho_fw = fp.density("water", P=self.P_fw_in.v, h=h_fw_out)

        if rho_fw > 0.0:
            Vol_dot_fw = self.m_dot_fw.v / rho_fw
        else:
            Vol_dot_fw = self.m_dot_fw.v / 1000.0

        # Set the Outputs from this Model (#,Value)
        self.m_dot_fw_out.v = self.m_dot_fw.v      # feedwater mass flow leaving the BFWH
        self.Vol_dot_fw.v = Vol_dot_fw              # Volumetric flow rate leaving the BFWH
        self.P_fw_out.v = self.P_fw_in.v            # Pressure of feedwater leaving the BFWH
        self.h_fw_out.v = h_fw_out                  # enthalpy of feedwater leaving the BFWH
        self.T_fw_out.v = T_fw_out                  # temperature of feedwater leaving the BFWH
        self.m_dot_total.v = m_dot_total            # total steam extraction mass flow leaving the feedwater heater (drain from previous and extraction from turbine combined)
        self.Vol_dot_total.v = Vol_dot_total        # total steam extraction volume flow leaving the feedwater heater (drain from previous and extraction from turbine combined)
        self.P_shell.v = P_shell                    # pressure of the extraction side leaving the feedwater heater
        self.h_b_out.v = h_b_out                    # enthalpy of the extraction side leaving the feedwater heater
        self.T_b_out.v = T_b_out                    # Temperature of the extraction leaving the feedwater heater
        self.Q_dot_act.v = Q_dot_act                # Total heat transfer from the turbine extraction to the feedwater
        self.m_dot_b.v = m_dot_b                    # requested turbine bleed from BFWH type
        self.DELTA_T_OD.v = DELTA_T_OD             # off-design temperature increase wanted by the BFWH - not always the case
        self.DELTA_T_act.v = DELTA_T_act            # actual feedwater temperature increase
