"""Deaerator-Pump component model (Type 6011)."""

import math

import eeslib.fluid_properties as fp

from esclab.simulate import Component
from esclab.components.flownetwork.esol6015_helpers import drhodhcp, drhodpch, dudhcp, dudpch, tank_level


class DeaeratorPump(Component):
    """
    Deaerator-Pump (esclab refactor of TRNSYS Type 6011).

    Three-pump deaerator system (one variable-speed pump, two constant-speed pumps)
    with a horizontal cylindrical tank model. Tank pressure, enthalpy, level, and mass
    are tracked via 4th-order Runge-Kutta integration. LP turbine steam extraction is
    computed for deaeration.

    Pump hydraulics are solved via the coupled-equations matrix.  Each active pump
    contributes a linearised head-curve equation; the matrix jointly solves pump
    operating points with the connected network.

    Parameters
    ----------
    pump_head_coeffs : ndarray (3, 3)
        [[A0, B0, C0], [A1, B1, C1], [A2, B2, C2]]
        Head curve per pump: dP = rho*g*(A*Q^2 + B*Q + C), Q in m^3/s, dP in Pa.
        Pump 0 is variable-speed; pumps 1-2 are constant-speed.
    pump_eta_coeffs : ndarray (3, 4)
        Efficiency curve per pump: eta = Ea*Q^4 + Eb*Q^3 + Ec*Q^2 + Ed*Q  [-]
    pump_npsh_coeffs : ndarray (3, 4)
        NPSH required curve per pump: NPSHr = Na*Q^3 + Nb*Q^2 + Nc*Q + Nd  [m]
        Coefficients for pump 0 are speed-scaled (Na*s^3, Nb*s^2, Nc*s).
    """
    trnsys_type = "6011"

    # -------------------------------------------------------------------------
    # Parameters
    # -------------------------------------------------------------------------

    pump_head_coeffs = Component.Parameter()    # (3x3) ndarray – head curve coefficients
    pump_eta_coeffs = Component.Parameter()     # (3x4) ndarray – efficiency curve coefficients
    pump_npsh_coeffs = Component.Parameter()    # (3x4) ndarray – NPSH curve coefficients

    D_tank = Component.Parameter()
    Length_tank = Component.Parameter()
    Length_tank2pump = Component.Parameter()
    P_tank_ini = Component.Parameter()
    L_tank_ini = Component.Parameter()
    m_dot_vent_frac = Component.Parameter()
    m_dot_LPB1_max = Component.Parameter()
    DA_ss_LPB1 = Component.Parameter()
    extraction_tol = Component.Parameter()
    LL_Alarm_sp = Component.Parameter()
    LL_Trip_sp = Component.Parameter()
    HL_Alarm_sp = Component.Parameter()
    HL_Trip_sp = Component.Parameter()

    # -------------------------------------------------------------------------
    # Inputs
    # -------------------------------------------------------------------------

    Turbine_ON = Component.Input()          # 1 = turbine on, 0 = off
    pump_power = Component.Input()          # (3,) array: 1/0 per pump (0=var-speed, 1,2=const-speed)
    pump_speed = Component.Input()          # fractional speed for pump 0 [-]
    m_dot_pump_in = Component.Input()       # loop-closure mass-flow from downstream [kg/s]
    m_dot_fw_in = Component.Input()
    P_fw_in = Component.Input()
    h_fw_in = Component.Input()
    m_dot_TB = Component.Input()
    P_TB = Component.Input()
    h_TB = Component.Input()
    m_dot_HPFWH = Component.Input()
    P_HPFWH = Component.Input()
    h_HPFWH = Component.Input()
    P_LPB1 = Component.Input()
    h_LPB1 = Component.Input()

    # -------------------------------------------------------------------------
    # Outputs
    # -------------------------------------------------------------------------

    # Combined pump discharge
    m_dot_pump_out = Component.Output()     # mass flow rate leaving all pumps [kg/s]
    Vol_dot_pump = Component.Output()       # volumetric flow rate leaving the pumps [m^3/s]
    P_pump_out = Component.Output()         # pressure leaving the pumps [Pa]
    h_pump_out = Component.Output()         # enthalpy leaving the pumps [J/kg]
    T_pump_out = Component.Output()         # temperature leaving the pumps [K]

    # Deaerator steam/bleed flows
    m_dot_LPB1 = Component.Output()         # extraction flow from low pressure turbine stage 1 [kg/s]
    m_dot_vent_out = Component.Output()     # mass of steam released from deaerator [kg/s]
    h_vent = Component.Output()             # enthalpy of steam released from deaerator [J/kg]

    # Tank state (written at convergence each timestep)
    m_tank = Component.Output()             # mass in the tank at the start of this timestep [kg]
    P_tank = Component.Output()             # pressure of the deaerator this timestep [Pa]
    L_tank = Component.Output()             # level of the tank this timestep [m]
    h_tank = Component.Output()             # enthalpy of the tank this timestep [J/kg]

    # Individual pump mass flows (solved by matrix)
    m_dot_P1 = Component.Output()           # mass flow rate out of pump 0 [kg/s]
    m_dot_P2 = Component.Output()           # mass flow rate out of pump 1 [kg/s]
    m_dot_P3 = Component.Output()           # mass flow rate out of pump 2 [kg/s]

    # Pump power and efficiency
    W_dot_total = Component.Output()        # total pump power [W]
    W_dot_P1 = Component.Output()           # power input to pump 0 [W]
    Eta_P1 = Component.Output()             # efficiency of pump 0 [-]
    W_dot_P2 = Component.Output()           # power input to pump 1 [W]
    Eta_P2 = Component.Output()             # efficiency of pump 1 [-]
    W_dot_P3 = Component.Output()           # power input to pump 2 [W]
    Eta_P3 = Component.Output()             # efficiency of pump 2 [-]

    # Pump cavitation trips
    P1_trip = Component.Output()            # pump 0 cavitation trip (1=cavitation, 0=ok) [-]
    P2_trip = Component.Output()            # pump 1 cavitation trip [-]
    P3_trip = Component.Output()            # pump 2 cavitation trip [-]

    # Level alarm/trip signals
    LL_Alarm_out = Component.Output()       # low level alarm signal (1=alarm, 0=ok) [-]
    LL_Trip_out = Component.Output()        # low level trip signal [-]
    HL_Alarm_out = Component.Output()       # high level alarm signal [-]
    HL_Trip_out = Component.Output()        # high level trip signal [-]

    # -------------------------------------------------------------------------

    def presim_setup(self, **kwargs):
        """Initialise instance-variable tank state and all output ports."""
        P_ini = self.P_tank_ini.v
        L_ini = self.L_tank_ini.v
        rho_g = fp.density("water", P=P_ini, Q=1.0)
        rho_f = fp.density("water", P=P_ini, Q=0.0)
        T_ini = fp.temperature("water", P=P_ini, Q=0.0)
        R_tank = self.D_tank.v / 2.0
        Vol_tank = math.pi * R_tank ** 2.0 * self.Length_tank.v
        Area_liq = (
            math.acos((R_tank - L_ini) / R_tank) * R_tank ** 2.0
            - (R_tank - L_ini) * math.sqrt(2.0 * R_tank * L_ini - L_ini ** 2.0)
        )
        m_f = Area_liq * self.Length_tank.v * rho_f
        m_g = (Vol_tank - Area_liq * self.Length_tank.v) * rho_g
        m_tot = m_f + m_g
        x_ini = m_g / m_tot
        h_ini = fp.enthalpy("water", P=P_ini, Q=x_ini)

        # Instance-variable tank state (pattern from Capacitor.U_C_prev)
        self._P_tank = P_ini
        self._h_tank = h_ini
        self._m_tank = m_tot
        self._L_tank = L_ini
        self._rho_tank_f = rho_f
        self._T_tank = T_ini

        # Minimal initial pump operating point (essentially zero flow through pump 0)
        P_bottom = P_ini + L_ini * rho_f * 9.81
        m_dot_ini = 1e-4
        Q_ini = m_dot_ini / rho_f
        h_p_in = fp.enthalpy("water", P=P_bottom, T=T_ini)
        s_p_in = fp.entropy("water", P=P_bottom, T=T_ini)
        A0, B0, C0 = self.pump_head_coeffs.v[0]
        P_out_ini = (A0 * Q_ini ** 2.0 + B0 * Q_ini + C0) * rho_f * 9.81 + P_bottom
        h_out_s = fp.enthalpy("water", P=P_out_ini, S=s_p_in)
        h_out_ini = h_p_in + (h_out_s - h_p_in) / 0.01   # minimal efficiency
        T_out_ini = fp.temperature("water", P=P_out_ini, H=h_out_ini)

        if self.Turbine_ON.v == 1.0:
            m_dot_LPB1_ini = self.m_dot_LPB1_max.v
        else:
            m_dot_LPB1_ini = 0.0

        LL_Alarm, LL_Trip, HL_Alarm, HL_Trip = self._check_alarms(L_ini)

        self.m_dot_pump_out.v = m_dot_ini
        self.Vol_dot_pump.v = Q_ini
        self.P_pump_out.v = P_out_ini
        self.h_pump_out.v = h_out_ini
        self.T_pump_out.v = T_out_ini
        self.m_dot_LPB1.v = m_dot_LPB1_ini
        self.m_dot_vent_out.v = 0.0
        self.h_vent.v = 0.0
        self.m_tank.v = m_tot
        self.P_tank.v = P_ini
        self.L_tank.v = L_ini
        self.h_tank.v = h_ini
        self.m_dot_P1.v = m_dot_ini
        self.m_dot_P2.v = 0.0
        self.m_dot_P3.v = 0.0
        self.W_dot_total.v = 0.0
        self.W_dot_P1.v = 0.0
        self.Eta_P1.v = 0.0
        self.W_dot_P2.v = 0.0
        self.Eta_P2.v = 0.0
        self.W_dot_P3.v = 0.0
        self.Eta_P3.v = 0.0
        self.P1_trip.v = 0.0
        self.P2_trip.v = 0.0
        self.P3_trip.v = 0.0
        self.LL_Alarm_out.v = LL_Alarm
        self.LL_Trip_out.v = LL_Trip
        self.HL_Alarm_out.v = HL_Alarm
        self.HL_Trip_out.v = HL_Trip

    def _check_alarms(self, level):
        """Return (LL_Alarm, LL_Trip, HL_Alarm, HL_Trip) for the given tank level."""
        if self.LL_Alarm_sp.v >= level:
            LL_Alarm = 1.0
            LL_Trip = 1.0 if self.LL_Trip_sp.v >= level else 0.0
        else:
            LL_Alarm = 0.0
            LL_Trip = 0.0
        if self.HL_Alarm_sp.v <= level:
            HL_Alarm = 1.0
            HL_Trip = 1.0 if self.HL_Trip_sp.v <= level else 0.0
        else:
            HL_Alarm = 0.0
            HL_Trip = 0.0
        return LL_Alarm, LL_Trip, HL_Alarm, HL_Trip

    def calculate(self):
        g = 9.81
        pump_outputs = [self.m_dot_P1, self.m_dot_P2, self.m_dot_P3]
        pump_speed_vals = [max(min(self.pump_speed.v, 1.0), 0.01), 1.0, 1.0]
        pump_on = self.pump_power.v
        head_coeffs = self.pump_head_coeffs.v
        P_bottom = self._P_tank + self._rho_tank_f * g * (self._L_tank + self.Length_tank2pump.v)

        # ------------------------------------------------------------------
        # Phase A: Matrix solve – linearised pump head curves
        # ------------------------------------------------------------------
        if self.coupled_eqs is not None:
            num_active = 0
            for i, (m_dot_pi, s_i) in enumerate(zip(pump_outputs, pump_speed_vals)):
                A, B, C = head_coeffs[i]
                if pump_on[i]:
                    num_active += 1
                    Q_i0 = max(m_dot_pi.v, 1e-6) / self._rho_tank_f
                    slope_i = g * (2.0 * A * Q_i0 + B * s_i)
                    intercept_i = self._rho_tank_f * g * (C * s_i ** 2.0 - A * Q_i0 ** 2.0)
                    self.coupled_eqs.add_equation(
                        {self.P_pump_out: 1.0, m_dot_pi: -slope_i},
                        rhs=P_bottom + intercept_i,
                    )
                else:
                    self.coupled_eqs.add_equation({m_dot_pi: 1.0}, rhs=0.0)

            # Flow combination: m_dot_pump_out = sum(m_dot_Pi)
            self.coupled_eqs.add_equation(
                {self.m_dot_pump_out: 1.0,
                 self.m_dot_P1: -1.0, self.m_dot_P2: -1.0, self.m_dot_P3: -1.0},
                rhs=0.0,
            )

            # 6th equation: loop closure when pumps run; static pressure when all off
            if num_active > 0:
                self.coupled_eqs.add_equation(
                    {self.m_dot_pump_out: 1.0,
                     self.coupled_eqs.source(self.m_dot_pump_in): -1.0},
                    rhs=0.0,
                )
            else:
                self.coupled_eqs.add_equation({self.P_pump_out: 1.0}, rhs=P_bottom)

            # Enthalpy: mass-averaged specific work at current operating point
            h_p_in = fp.enthalpy("water", P=P_bottom, T=self._T_tank)
            s_p_in = fp.entropy("water", P=P_bottom, T=self._T_tank)
            W_dot_sum = 0.0
            m_dot_total = sum(max(m.v, 0.0) for m in pump_outputs)
            if m_dot_total > 0.0:
                for i, (m_dot_pi, s_i) in enumerate(zip(pump_outputs, pump_speed_vals)):
                    if pump_on[i] and m_dot_pi.v > 0.0:
                        A, B, C = head_coeffs[i]
                        Q_i = m_dot_pi.v / self._rho_tank_f
                        P_i_out = (A * Q_i ** 2.0 + B * s_i * Q_i + C * s_i ** 2.0) * self._rho_tank_f * g + P_bottom
                        h_i_s = fp.enthalpy("water", P=P_i_out, S=s_p_in)
                        ea, eb, ec, ed = self.pump_eta_coeffs.v[i]
                        eta_i = max(ea * Q_i ** 4.0 + eb * Q_i ** 3.0 + ec * Q_i ** 2.0 + ed * Q_i, 0.01)
                        W_dot_sum += m_dot_pi.v * (h_i_s - h_p_in) / eta_i
            W_dot_specific = W_dot_sum / m_dot_total if m_dot_total > 0.0 else 0.0
            self.coupled_eqs.add_equation({self.h_pump_out: 1.0}, rhs=h_p_in + W_dot_specific)
            return

        # ------------------------------------------------------------------
        # Phase B: Sequential – LP extraction, vent, per-pump diagnostics
        # ------------------------------------------------------------------
        P_tank = self._P_tank
        h_tank = self._h_tank
        m_tank = self._m_tank
        L_tank = self._L_tank
        rho_f = self._rho_tank_f

        # Step 5: LP extraction
        if self.Turbine_ON.v == 1.0:
            m_dot_LPB1_prev = self.m_dot_LPB1.v
            h_sat_da = fp.enthalpy("water", P=P_tank, Q=0.05)
            m_dot_LPB1 = (
                self.m_dot_fw_in.v * (self.h_fw_in.v - h_sat_da)
                + self.m_dot_TB.v * (self.h_TB.v - h_sat_da)
                + self.m_dot_HPFWH.v * (self.h_HPFWH.v - h_sat_da)
            ) / (h_sat_da - self.h_LPB1.v)
            m_dot_LPB1 = max(min(m_dot_LPB1, self.m_dot_LPB1_max.v), 0.0)
            delta_m = abs(m_dot_LPB1 - m_dot_LPB1_prev)
            if self.model.time != self.model.timestep or self.model.iteration != 0:
                if delta_m > self.extraction_tol.v:
                    if m_dot_LPB1 > m_dot_LPB1_prev:
                        m_dot_LPB1 = min(
                            m_dot_LPB1_prev + delta_m * 0.1,
                            m_dot_LPB1_prev + self.DA_ss_LPB1.v,
                            self.m_dot_LPB1_max.v,
                        )
                    else:
                        m_dot_LPB1 = max(
                            m_dot_LPB1_prev - delta_m * 0.1,
                            m_dot_LPB1_prev - self.DA_ss_LPB1.v,
                            0.0,
                        )
                else:
                    m_dot_LPB1 = m_dot_LPB1_prev
        else:
            m_dot_LPB1 = 0.0

        # Step 6: inlet mixing and vent enthalpy
        m_dot_in = self.m_dot_fw_in.v + self.m_dot_TB.v + self.m_dot_HPFWH.v + m_dot_LPB1
        if m_dot_in > 0.0:
            h_in_mix = (
                self.m_dot_fw_in.v * self.h_fw_in.v
                + self.m_dot_TB.v * self.h_TB.v
                + self.m_dot_HPFWH.v * self.h_HPFWH.v
                + m_dot_LPB1 * self.h_LPB1.v
            ) / m_dot_in
        else:
            h_in_mix = h_tank
        h_vent = fp.enthalpy("water", P=P_tank, Q=1.0)

        # Per-pump diagnostics (m_dot_Pi already solved by the matrix)
        h_p_in = fp.enthalpy("water", P=P_bottom, T=self._T_tank)
        s_p_in = fp.entropy("water", P=P_bottom, T=self._T_tank)
        pump_W_outs = [self.W_dot_P1, self.W_dot_P2, self.W_dot_P3]
        pump_eta_outs = [self.Eta_P1, self.Eta_P2, self.Eta_P3]
        h_out_vals = []
        P_out_vals = []
        W_total = 0.0
        for i, (m_dot_pi, W_pi, eta_pi, s_i) in enumerate(
            zip(pump_outputs, pump_W_outs, pump_eta_outs, pump_speed_vals)
        ):
            if pump_on[i] and m_dot_pi.v > 0.0:
                A, B, C = head_coeffs[i]
                Q_i = m_dot_pi.v / rho_f
                P_i_out = (A * Q_i ** 2.0 + B * s_i * Q_i + C * s_i ** 2.0) * rho_f * g + P_bottom
                h_i_s = fp.enthalpy("water", P=P_i_out, S=s_p_in)
                ea, eb, ec, ed = self.pump_eta_coeffs.v[i]
                eta_i = max(ea * Q_i ** 4.0 + eb * Q_i ** 3.0 + ec * Q_i ** 2.0 + ed * Q_i, 0.01)
                W_i = m_dot_pi.v * (h_i_s - h_p_in) / eta_i
                eta_pi.v = eta_i
                W_pi.v = W_i
                W_total += W_i
                h_out_vals.append((m_dot_pi.v, h_p_in + W_i / m_dot_pi.v))
                P_out_vals.append((m_dot_pi.v, P_i_out))
            else:
                eta_pi.v = 0.0
                W_pi.v = 0.0
                h_out_vals.append((0.0, h_p_in))
                P_out_vals.append((0.0, P_bottom))

        self.W_dot_total.v = W_total
        m_dot_pump = self.m_dot_pump_out.v
        if m_dot_pump > 0.0:
            h_combined = sum(m * h for m, h in h_out_vals) / m_dot_pump
            P_combined = sum(m * P for m, P in P_out_vals) / m_dot_pump
            T_combined = fp.temperature("water", P=P_combined, H=h_combined)
            Vol_dot = m_dot_pump / rho_f
        else:
            h_combined = h_p_in
            P_combined = P_bottom
            T_combined = self._T_tank
            Vol_dot = 0.0
        self.Vol_dot_pump.v = Vol_dot
        self.T_pump_out.v = T_combined
        self.m_dot_vent_out.v = self.m_dot_vent_frac.v
        self.h_vent.v = h_vent
        self.m_dot_LPB1.v = m_dot_LPB1

        if not self.model.is_converged:
            return

        # ------------------------------------------------------------------
        # Phase C: Convergence – RK4 tank integration, NPSH, alarms, state advance
        # ------------------------------------------------------------------
        m_dot_pump_val = m_dot_pump
        m_dot_vent = self.m_dot_vent_frac.v
        ts = self.model.timestep * 3600.0
        dh_fd = 1000.0
        dP_fd = 1000.0
        Vol_tank = (
            3.14 / 4.0 * self.D_tank.v ** 2.0 * (self.Length_tank.v - self.D_tank.v)
            + 4.0 / 3.0 * 3.14 * (self.D_tank.v / 2.0) ** 3.0
        )

        def _rk4_rates(P, h):
            _drhodh = drhodhcp(P_tank=P, h_tank=h, dh=dh_fd)
            _drhodP = drhodpch(P_tank=P, h_tank=h, dP=dP_fd)
            _dudh = dudhcp(P_tank=P, h_tank=h, dh=dh_fd)
            _dudP = dudpch(P_tank=P, h_tank=h, dP=dP_fd)
            rho = m_tank / Vol_tank
            u = fp.internalenergy("water", P=P, H=h)
            h_sf = fp.enthalpy("water", P=P, Q=0.0)
            h_sv = fp.enthalpy("water", P=P, Q=1.0)
            denom = _dudh * _drhodP - _dudP * _drhodh
            denom = math.copysign(max(abs(denom), 0.000025), denom)
            _drhodh = math.copysign(max(abs(_drhodh), 0.0001), _drhodh)
            dPdt = (
                ((u - h_in_mix) * m_dot_in + (-u + h_sf) * m_dot_pump_val
                 - m_dot_vent * (u - h_sv)) * _drhodh
                + rho * _dudh * (m_dot_in - m_dot_pump_val - m_dot_vent)
            ) / (m_tank * (_dudh * _drhodP - _dudP * _drhodh))
            dhdt = ((m_dot_in - m_dot_pump_val - m_dot_vent) / Vol_tank - _drhodP * dPdt) / _drhodh
            return dPdt, dhdt

        dPdt_aa, dhdt_aa = _rk4_rates(P_tank, h_tank)
        P_aa = P_tank + dPdt_aa * ts / 2.0
        h_aa = h_tank + dhdt_aa * ts / 2.0

        dPdt_bb, dhdt_bb = _rk4_rates(P_aa, h_aa)
        P_bb = P_tank + dPdt_bb * ts / 2.0
        h_bb = h_tank + dhdt_bb * ts / 2.0

        dPdt_cc, dhdt_cc = _rk4_rates(P_bb, h_bb)
        P_cc = P_tank + dPdt_cc * ts
        h_cc = h_tank + dhdt_cc * ts

        dPdt_dd, dhdt_dd = _rk4_rates(P_cc, h_cc)

        P_new = P_tank + (dPdt_aa + 2.0 * dPdt_bb + 2.0 * dPdt_cc + dPdt_dd) * ts / 6.0
        h_new = h_tank + (dhdt_aa + 2.0 * dhdt_bb + 2.0 * dhdt_cc + dhdt_dd) * ts / 6.0

        x_new = fp.quality("water", P=P_new, H=h_new)
        T_new = fp.temperature("water", P=P_new, H=h_new)
        rho_f_new = fp.density("water", P=P_new, Q=0.0)
        dmdt = m_dot_in - m_dot_pump_val - m_dot_vent
        m_new = m_tank + dmdt * ts
        m_g_new = m_new * x_new
        m_f_new = m_new - max(m_g_new, 0.0)
        Vol_liq = m_f_new / rho_f_new
        L_new = tank_level(Vol_liq, self.D_tank.v, self.Length_tank.v, L_tank, 0.01)

        # NPSH cavitation check
        P_ref = 87726.1
        T_ref = 373.0
        D_inlet = 0.1524
        P_pump_inlet = P_tank + (L_tank + self.Length_tank2pump.v) * rho_f * g
        T_v = fp.temperature("water", P=P_tank, Q=1.0)
        lnP1P2 = 8.314 * (1.0 / T_ref - 1.0 / T_v)
        P_v = P_ref * math.exp(lnP1P2)
        trip_outs = [self.P1_trip, self.P2_trip, self.P3_trip]
        npsh_coeffs = self.pump_npsh_coeffs.v
        for i, (m_dot_pi, trip_pi, s_i) in enumerate(
            zip(pump_outputs, trip_outs, pump_speed_vals)
        ):
            if pump_on[i]:
                Q_i = m_dot_pi.v / rho_f
                vel = Q_i / (math.pi / 4.0 * D_inlet ** 2.0)
                na, nb, nc, nd = npsh_coeffs[i]
                NPSHr = (na * s_i ** 3.0) * Q_i ** 3.0 + (nb * s_i ** 2.0) * Q_i ** 2.0 + (nc * s_i) * Q_i + nd
                NPSHa = P_pump_inlet / (rho_f * g) + vel ** 2.0 / (2.0 * g)
                trip_pi.v = 0.0 if NPSHa > NPSHr else 1.0
            else:
                trip_pi.v = 0.0

        # Alarm check
        LL_Alarm, LL_Trip, HL_Alarm, HL_Trip = self._check_alarms(L_tank)
        self.LL_Alarm_out.v = LL_Alarm
        self.LL_Trip_out.v = LL_Trip
        self.HL_Alarm_out.v = HL_Alarm
        self.HL_Trip_out.v = HL_Trip

        # Advance instance-variable state and write diagnostic output ports
        self._P_tank = P_new
        self._h_tank = h_new
        self._m_tank = m_new
        self._L_tank = L_new
        self._rho_tank_f = rho_f_new
        self._T_tank = T_new

        self.m_tank.v = m_new
        self.P_tank.v = P_new
        self.L_tank.v = L_new
        self.h_tank.v = h_new













