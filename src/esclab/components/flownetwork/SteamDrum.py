"""Type 6003 steam drum converted from Fortran."""

import math
from collections import deque

from eeslib import fluid_properties as fp

from esclab.simulate import Component


class SteamDrum(Component):
    """
    TRNSYS Type 6003: Steam Drum (ESOL6003-SteamDrum).

    Parameters
    ----------
    d_tank : float
        Parameter 1, tank diameter [m].
    length_tank : float
        Parameter 2, tank length [m].
    initial_l_tank : float
        Parameter 3, initial tank level [m].
    initial_p_tank : float
        Parameter 4, initial tank pressure [Pa].
    rated_heat_transfer : float
        Parameter 5, rated heat-transfer coefficient input used for UA scaling.
    rated_htf_flow : float
        Parameter 6, rated HTF mass flow [kg/s].
    rated_exp : float
        Parameter 7, exponent for off-design UA relation.
    no_shell_passes : float
        Parameter 8, number of shell passes.
    no_tube_passes : float
        Parameter 9, number of tube passes.
    length_hx : float
        Parameter 10, evaporator tube length [m].
    tube_od : float
        Parameter 11, tube outside diameter [m].
    tube_th : float
        Parameter 12, tube wall thickness [m].
    no_tubes : float
        Parameter 13, number of tubes.
    fluid_id : float | str
        Parameter 14, HTF identifier for property calls.
    low_level_alarm_cond..high_htf_hr_trip_cond : float
        Parameters 15-28, alarm/trip thresholds.

    Inputs
    ------
    m_dot_in : float
        Input 1, feedwater mass flow into drum [kg/s].
    h_in : float
        Input 2, feedwater specific enthalpy [J/kg].
    p_in : float
        Input 3, feedwater pressure [Pa].
    htf_mass_in : float
        Input 4, HTF mass flow into evaporator [kg/s].
    htf_temp_in : float
        Input 5, HTF inlet temperature [K].
    htf_p_in : float
        Input 6, HTF inlet pressure [Pa].
    m_dot_superheat_req : float
        Input 7, requested steam flow to superheater [kg/s].

    Outputs
    -------
    m_dot_superheat..high_htf_hr_trip : float
        Outputs 1-36 mapped directly to the Type6003 `setOutputValue` ordering.

    Notes
    -----
    The original Fortran uses `drhodhcp`, `drhodpch`, `dudhcp`, `dudpch`, and
    `tank_level` from `ESOL6015_myfunctions`. This conversion preserves the
    Type6003 RK4 pressure/enthalpy structure by evaluating those derivatives
    numerically from water properties.
    """

    # Parameters (Type6003: 1-28)
    d_tank = Component.Parameter()
    length_tank = Component.Parameter()
    initial_l_tank = Component.Parameter()
    initial_p_tank = Component.Parameter()
    rated_heat_transfer = Component.Parameter()
    rated_htf_flow = Component.Parameter()
    rated_exp = Component.Parameter()
    no_shell_passes = Component.Parameter()
    no_tube_passes = Component.Parameter()
    length_hx = Component.Parameter()
    tube_od = Component.Parameter()
    tube_th = Component.Parameter()
    no_tubes = Component.Parameter()
    fluid_id = Component.Parameter()
    low_level_alarm_cond = Component.Parameter()
    low_level_trip_cond = Component.Parameter()
    high_level_alarm_cond = Component.Parameter()
    high_level_trip_cond = Component.Parameter()
    high_pressure_alarm_cond = Component.Parameter()
    high_pressure_trip_cond = Component.Parameter()
    high_delta_t_fw_alarm_cond = Component.Parameter()
    high_delta_t_fw_trip_cond = Component.Parameter()
    high_tank_hr_alarm_cond = Component.Parameter()
    high_tank_hr_trip_cond = Component.Parameter()
    high_delta_t_htf_alarm_cond = Component.Parameter()
    high_delta_t_htf_trip_cond = Component.Parameter()
    high_htf_hr_alarm_cond = Component.Parameter()
    high_htf_hr_trip_cond = Component.Parameter()

    # Inputs (Type6003: 1-7)
    m_dot_in = Component.Input()
    h_in = Component.Input()
    p_in = Component.Input()
    htf_mass_in = Component.Input()
    htf_temp_in = Component.Input()
    htf_p_in = Component.Input()
    m_dot_superheat_req = Component.Input()

    # Outputs (Type6003: 1-36)
    m_dot_superheat = Component.Output()
    vol_dot_fw = Component.Output()
    p_superheat = Component.Output()
    t_superheat = Component.Output()
    h_superheat = Component.Output()
    htf_mass_out = Component.Output()
    vol_dot_htf = Component.Output()
    htf_temp_out = Component.Output()
    htf_p_out = Component.Output()
    m_tank_end = Component.Output()
    l_tank_end = Component.Output()
    p_tank_end = Component.Output()
    t_tank_end = Component.Output()
    h_tank_end = Component.Output()
    m_tank_begin = Component.Output()
    l_tank_begin = Component.Output()
    p_tank_begin = Component.Output()
    h_tank_begin = Component.Output()
    evaporator_effectiveness = Component.Output()
    m_dot_evap = Component.Output()
    tank_heating_rate = Component.Output()
    htf_inlet_heating_rate = Component.Output()
    low_level_alarm = Component.Output()
    low_level_trip = Component.Output()
    high_level_alarm = Component.Output()
    high_level_trip = Component.Output()
    high_pressure_alarm = Component.Output()
    high_pressure_trip = Component.Output()
    high_delta_t_fw_alarm = Component.Output()
    high_delta_t_fw_trip = Component.Output()
    high_tank_hr_alarm = Component.Output()
    high_tank_hr_trip = Component.Output()
    high_delta_t_htf_alarm = Component.Output()
    high_delta_t_htf_trip = Component.Output()
    high_htf_hr_alarm = Component.Output()
    high_htf_hr_trip = Component.Output()

    @staticmethod
    def _safe(value, default=0.0):
        return value if value == value else default

    @staticmethod
    def _safe_pow(base, exponent):
        if base <= 0.0:
            return 0.0
        return base**exponent

    @staticmethod
    def _hr_from_history(history_values, current_value, dt_seconds):
        n_points = len(history_values)
        if n_points <= 0:
            return 0.0
        hr_sum = 0.0
        for index in range(1, n_points):
            hr_sum += (history_values[index] - history_values[index - 1]) / dt_seconds
        hr_sum += (current_value - history_values[-1]) / dt_seconds
        return hr_sum / n_points * 60.0

    @staticmethod
    def _tank_total_volume(d_tank, length_tank):
        return (
            math.pi * d_tank**2.0 / 4.0 * (length_tank - d_tank)
            + 4.0 / 3.0 * math.pi * (d_tank / 2.0) ** 3.0
        )

    @staticmethod
    def _liquid_volume_from_level(level, d_tank, length_tank):
        level_limited = max(0.0, min(level, d_tank))
        radius = d_tank / 2.0
        asin_arg = max(-1.0, min(1.0, 1.0 - 2.0 * level_limited / max(d_tank, 1.0e-9)))
        root_term = max(radius**2.0 - (radius - level_limited) ** 2.0, 0.0)

        area_liquid = (
            math.pi * d_tank**2.0 / 4.0 * (0.5 - (math.asin(asin_arg) / math.pi))
            - math.sqrt(root_term) * (radius - level_limited)
        )
        vol_liquid = (
            area_liquid * (length_tank - d_tank)
            + math.pi * (level_limited**2.0 * d_tank / 2.0 - level_limited**3.0 / 2.0)
        )
        return max(vol_liquid, 0.0)

    def _tank_level_from_volume(self, target_volume, d_tank, length_tank, level_guess):
        low = 0.0
        high = max(d_tank, 1.0e-9)
        level = max(min(level_guess, high), low)
        for _ in range(60):
            mid = 0.5 * (low + high)
            vol_mid = self._liquid_volume_from_level(mid, d_tank, length_tank)
            if abs(vol_mid - target_volume) < 1.0e-6:
                level = mid
                break
            if vol_mid < target_volume:
                low = mid
            else:
                high = mid
            level = mid
        return max(min(level, d_tank), 0.0)

    def _sat_props_from_pressure(self, pressure_pa):
        pressure_eval = max(pressure_pa, 1.0)

        t_sat = None
        h_f = None
        h_g = None
        rho_f = None
        rho_g = None

        for p_value in (pressure_eval, pressure_eval / 1000.0):
            try:
                t_sat = float(fp.temperature("water", P=p_value, Q=1.0))
            except Exception:
                pass
            try:
                h_f_raw = float(fp.enthalpy("water", P=p_value, Q=0.0))
                h_g_raw = float(fp.enthalpy("water", P=p_value, Q=1.0))
                h_f = h_f_raw * 1000.0 if abs(h_f_raw) < 1.0e4 else h_f_raw
                h_g = h_g_raw * 1000.0 if abs(h_g_raw) < 1.0e4 else h_g_raw
            except Exception:
                pass
            try:
                rho_f = float(fp.density("water", P=p_value, Q=0.0))
                rho_g = float(fp.density("water", P=p_value, Q=1.0))
            except Exception:
                pass

        if t_sat is None or h_f is None or h_g is None or rho_f is None or rho_g is None:
            pressure_ratio = pressure_eval / 101325.0
            t_sat = max(273.15, min(650.0, 373.15 + 42.0 * math.log(pressure_ratio + 1.0)))
            h_f = 4200.0 * (t_sat - 273.15)
            h_fg = max(2.5e6 - 1800.0 * (t_sat - 273.15), 5.0e5)
            h_g = h_f + h_fg
            rho_f = max(1000.0 - 0.35 * (t_sat - 273.15), 600.0)
            rho_g = max(pressure_eval / (461.5 * max(t_sat, 200.0)), 0.05)

        return t_sat, h_f, h_g, max(rho_f, 1.0e-6), max(rho_g, 1.0e-6)

    def _temperature_from_ph(self, pressure_pa, enthalpy_j_kg):
        pressure_eval = max(pressure_pa, 1.0)
        h_eval = enthalpy_j_kg
        for kwargs in (
            {"fluid": "water", "P": pressure_eval, "H": h_eval},
            {"fluid": "water", "P": pressure_eval / 1000.0, "H": h_eval / 1000.0},
            {"fluid": "water", "P": pressure_eval, "h": h_eval},
            {"fluid": "water", "P": pressure_eval / 1000.0, "h": h_eval / 1000.0},
        ):
            try:
                temperature = float(fp.temperature(**kwargs))
                if temperature == temperature and 200.0 <= temperature <= 2000.0:
                    return temperature
            except Exception:
                continue

        t_sat, h_f, h_g, _, _ = self._sat_props_from_pressure(pressure_eval)
        if h_eval <= h_f:
            return max(273.15, min(t_sat - 0.5, 273.15 + h_eval / 4200.0))
        if h_eval >= h_g:
            return min(1200.0, t_sat + (h_eval - h_g) / 2100.0)
        return t_sat

    def _cp_htf(self, fluid_id, temperature_k, pressure_pa):
        t_eval = max(temperature_k, 273.15)
        p_eval = max(pressure_pa, 1.0)
        try:
            cp = float(fp.specheat(fluid_id, T=t_eval, P=p_eval))
            if cp == cp and cp > 0.0:
                return cp * 1000.0 if cp < 100.0 else cp
        except Exception:
            pass
        return 2200.0

    def _rho_htf(self, fluid_id, temperature_k, pressure_pa):
        try:
            rho = float(fp.density(fluid_id, T=max(temperature_k, 273.15), P=max(pressure_pa, 1.0)))
            if rho == rho and rho > 0.0:
                return rho
        except Exception:
            pass
        return 1000.0

    def _rho_water_from_ph(self, pressure_pa, enthalpy_j_kg):
        pressure_eval = max(pressure_pa, 1.0)
        h_eval = enthalpy_j_kg

        for kwargs in (
            {"fluid": "water", "P": pressure_eval, "H": h_eval},
            {"fluid": "water", "P": pressure_eval / 1000.0, "H": h_eval / 1000.0},
            {"fluid": "water", "P": pressure_eval, "h": h_eval},
            {"fluid": "water", "P": pressure_eval / 1000.0, "h": h_eval / 1000.0},
        ):
            try:
                rho = float(fp.density(**kwargs))
                if rho == rho and rho > 0.0:
                    return rho
            except Exception:
                continue

        t_sat, h_f, h_g, rho_f, rho_g = self._sat_props_from_pressure(pressure_eval)
        if h_eval <= h_f:
            t_est = max(273.15, min(t_sat, 273.15 + h_eval / 4200.0))
            return max(600.0, 1000.0 - 0.35 * (t_est - 273.15))
        if h_eval >= h_g:
            t_est = t_sat + (h_eval - h_g) / 2100.0
            return max(0.05, pressure_eval / (461.5 * max(t_est, 200.0)))

        quality = (h_eval - h_f) / max(h_g - h_f, 1.0e-9)
        return 1.0 / max((1.0 - quality) / rho_f + quality / rho_g, 1.0e-12)

    def _u_water_from_ph(self, pressure_pa, enthalpy_j_kg):
        rho = self._rho_water_from_ph(pressure_pa, enthalpy_j_kg)
        return enthalpy_j_kg - pressure_pa / max(rho, 1.0e-9)

    def _state_derivatives(self, pressure_pa, enthalpy_j_kg, m_dot_in, m_dot_superheat, m_dot_blowdown, h_in, q_dot_actual, vol_tank):
        pressure_eval = max(pressure_pa, 1.0)
        enthalpy_eval = max(enthalpy_j_kg, 1.0)

        dp = 1000.0
        dh = 1000.0

        rho = self._rho_water_from_ph(pressure_eval, enthalpy_eval)
        u_tank = self._u_water_from_ph(pressure_eval, enthalpy_eval)
        _, h_sat_f, h_sat_g, _, _ = self._sat_props_from_pressure(pressure_eval)

        rho_h_plus = self._rho_water_from_ph(pressure_eval, enthalpy_eval + dh)
        rho_h_minus = self._rho_water_from_ph(pressure_eval, max(enthalpy_eval - dh, 1.0))
        drhodhcp = (rho_h_plus - rho_h_minus) / max((2.0 * dh), 1.0e-9)

        rho_p_plus = self._rho_water_from_ph(pressure_eval + dp, enthalpy_eval)
        rho_p_minus = self._rho_water_from_ph(max(pressure_eval - dp, 1.0), enthalpy_eval)
        drhodpch = (rho_p_plus - rho_p_minus) / max((2.0 * dp), 1.0e-9)

        u_h_plus = self._u_water_from_ph(pressure_eval, enthalpy_eval + dh)
        u_h_minus = self._u_water_from_ph(pressure_eval, max(enthalpy_eval - dh, 1.0))
        dudhcp = (u_h_plus - u_h_minus) / max((2.0 * dh), 1.0e-9)

        u_p_plus = self._u_water_from_ph(pressure_eval + dp, enthalpy_eval)
        u_p_minus = self._u_water_from_ph(max(pressure_eval - dp, 1.0), enthalpy_eval)
        dudpch = (u_p_plus - u_p_minus) / max((2.0 * dp), 1.0e-9)

        denominator = drhodhcp * dudpch - drhodpch * dudhcp
        denominator = math.copysign(max(abs(denominator), 1.0e-5), denominator if denominator != 0.0 else 1.0)
        drhodhcp_safe = math.copysign(max(abs(drhodhcp), 5.0e-5), drhodhcp if drhodhcp != 0.0 else 1.0)

        dpdt = (
            ((u_tank - h_sat_g) * m_dot_superheat + (u_tank - h_sat_f) * m_dot_blowdown + (h_in - u_tank) * m_dot_in + q_dot_actual) * drhodhcp_safe
            + dudhcp * rho * (m_dot_blowdown + m_dot_superheat - m_dot_in)
        ) / max(vol_tank * rho * denominator, 1.0e-12)

        dhdt = ((m_dot_in - m_dot_superheat - m_dot_blowdown) / max(vol_tank, 1.0e-12) - drhodpch * dpdt) / drhodhcp_safe
        return dpdt, dhdt

    def _initialize_hr_history(self, n_int, t_tank_value, htf_temp_value):
        self._n_int = max(int(n_int), 1)
        self._tank_temp_hist = deque([t_tank_value] * self._n_int, maxlen=self._n_int)
        self._htf_temp_hist = deque([htf_temp_value] * self._n_int, maxlen=self._n_int)

    def calculate(self):
        d_tank = max(self._safe(self.d_tank.v, 1.0), 1.0e-3)
        length_tank = max(self._safe(self.length_tank.v, 1.0), 1.0e-3)

        rated_heat_transfer = max(self._safe(self.rated_heat_transfer.v, 0.0), 0.0)
        rated_htf_flow = max(self._safe(self.rated_htf_flow.v, 1.0), 1.0e-9)
        rated_exp = self._safe(self.rated_exp.v, 0.8)
        no_shell_passes = max(self._safe(self.no_shell_passes.v, 1.0), 1.0)
        no_tube_passes = max(self._safe(self.no_tube_passes.v, 1.0), 1.0)
        length_hx = max(self._safe(self.length_hx.v, 1.0), 1.0e-6)
        tube_od = max(self._safe(self.tube_od.v, 0.01), 1.0e-6)
        tube_th = max(self._safe(self.tube_th.v, 0.001), 0.0)
        no_tubes = max(self._safe(self.no_tubes.v, 1.0), 1.0)
        fluid_id = self._safe(self.fluid_id.v, "water")

        m_dot_in = max(self._safe(self.m_dot_in.v), 0.0)
        h_in = self._safe(self.h_in.v)
        p_in = max(self._safe(self.p_in.v, 101325.0), 1.0)
        htf_mass_in = max(self._safe(self.htf_mass_in.v), 0.0)
        htf_temp_in = self._safe(self.htf_temp_in.v, 300.0)
        htf_p_in = self._safe(self.htf_p_in.v, p_in)
        m_dot_superheat_req = max(self._safe(self.m_dot_superheat_req.v, 60.0), 0.0)

        timestep_seconds = max(self.model.settings.timestep * 3600.0, 1.0e-9)
        n_int = max(int(math.ceil(60.0 / timestep_seconds)), 1)
        if not hasattr(self, "_n_int"):
            self._initialize_hr_history(n_int, htf_temp_in, htf_temp_in)
        elif self._n_int != n_int:
            tank_seed = self.t_tank_end.v if self.t_tank_end.v == self.t_tank_end.v else htf_temp_in
            self._initialize_hr_history(n_int, tank_seed, htf_temp_in)

        vol_tank = self._tank_total_volume(d_tank, length_tank)

        # Do All of the First Timestep Manipulations Here - There Are No Iterations at the Intial Time
        if self.model.is_first_step:
            l_tank = max(min(self._safe(self.initial_l_tank.v, 0.5 * d_tank), d_tank), 0.0)
            p_tank = max(self._safe(self.initial_p_tank.v, p_in), 1.0)
            t_tank, h_sat_f, h_out, rho_tank_f, rho_tank_g = self._sat_props_from_pressure(p_tank)

            vol_liquid = self._liquid_volume_from_level(l_tank, d_tank, length_tank)
            vol_vapor = max(vol_tank - vol_liquid, 0.0)

            m_tank_f = vol_liquid * rho_tank_f
            m_tank_g = vol_vapor * rho_tank_g
            m_tank_tot = max(m_tank_f + m_tank_g, 1.0)
            x_tank = m_tank_g / m_tank_tot
            h_tank = h_sat_f + x_tank * (h_out - h_sat_f)

            # inital guess for steam flow leaving
            m_dot_superheat = 60.0
            vol_dot_fw = m_dot_superheat / max(rho_tank_g, 1.0e-9)

            # Volumetric flow rate for HTF
            rho_htf = self._rho_htf(fluid_id, htf_temp_in, htf_p_in)
            vol_dot_htf = htf_mass_in / max(rho_htf, 1.0e-9)

            # Set the Initial Values of the Outputs (#,Value)
            self.m_dot_superheat.v = m_dot_superheat
            self.vol_dot_fw.v = vol_dot_fw
            self.p_superheat.v = p_tank
            self.t_superheat.v = t_tank
            self.h_superheat.v = h_out
            self.htf_mass_out.v = htf_mass_in
            self.vol_dot_htf.v = vol_dot_htf
            self.htf_temp_out.v = htf_temp_in
            self.htf_p_out.v = htf_p_in
            self.m_tank_end.v = m_tank_tot
            self.l_tank_end.v = l_tank
            self.p_tank_end.v = p_tank
            self.t_tank_end.v = t_tank
            self.h_tank_end.v = h_tank
            self.m_tank_begin.v = m_tank_tot
            self.l_tank_begin.v = l_tank
            self.p_tank_begin.v = p_tank
            self.h_tank_begin.v = h_tank

            self.evaporator_effectiveness.v = 0.0
            self.m_dot_evap.v = 0.0
            self.tank_heating_rate.v = 0.0
            self.htf_inlet_heating_rate.v = 0.0
            self.low_level_alarm.v = 0.0
            self.low_level_trip.v = 0.0
            self.high_level_alarm.v = 0.0
            self.high_level_trip.v = 0.0
            self.high_pressure_alarm.v = 0.0
            self.high_pressure_trip.v = 0.0
            self.high_delta_t_fw_alarm.v = 0.0
            self.high_delta_t_fw_trip.v = 0.0
            self.high_tank_hr_alarm.v = 0.0
            self.high_tank_hr_trip.v = 0.0
            self.high_delta_t_htf_alarm.v = 0.0
            self.high_delta_t_htf_trip.v = 0.0
            self.high_htf_hr_alarm.v = 0.0
            self.high_htf_hr_trip.v = 0.0

            self._initialize_hr_history(n_int, t_tank, htf_temp_in)
            return

        # Read previous-step storage outputs
        m_tank_prev = max(self._safe(self.m_tank_begin.v, 1.0), 1.0)
        p_tank_prev = max(self._safe(self.p_tank_begin.v, p_in), 1.0)
        l_tank_prev = max(self._safe(self.l_tank_begin.v, 0.0), 0.0)
        h_tank_prev = self._safe(self.h_tank_begin.v, h_in)

        t_tank_prev, h_sat_f_prev, h_sat_g_prev, _, rho_g_prev = self._sat_props_from_pressure(p_tank_prev)

        # Evaporator calculations
        htf_mass_eval = max(htf_mass_in, 1.0e-5)
        a_s = math.pi * max(tube_od - 2.0 * tube_th, 1.0e-6) * length_hx * no_tube_passes * no_tubes
        ua_rated = rated_heat_transfer * a_s
        ua_od = ua_rated * self._safe_pow(htf_mass_eval / rated_htf_flow, rated_exp)

        cp_htf_max = self._cp_htf(fluid_id, htf_temp_in, htf_p_in)
        cp_htf_min = self._cp_htf(fluid_id, t_tank_prev, htf_p_in)
        cp_htf_ave = max((cp_htf_max + cp_htf_min) / 2.0, 1.0)

        ntu_od = ua_od / max(htf_mass_eval * cp_htf_ave, 1.0e-9)
        if ntu_od < 1.0e-10:
            eta_1pass = 0.0
        else:
            exp_term = math.exp(-ntu_od)
            eta_1pass = min(
                2.0 * (1.0 + (1.0 + exp_term) / max(1.0 - exp_term, 1.0e-12)) ** -1.0,
                0.99999999,
            )

        cr = 0.0
        ratio = (1.0 - eta_1pass * cr) / max(1.0 - eta_1pass, 1.0e-12)
        eta_od = (ratio**no_shell_passes - 1.0) / max(ratio**no_shell_passes - cr, 1.0e-12)
        eta_od = max(min(eta_od, 1.0), 0.0)

        q_dot_max = htf_mass_eval * cp_htf_ave * (htf_temp_in - t_tank_prev)
        q_dot_actual = eta_od * q_dot_max

        htf_mass_out = htf_mass_in
        htf_temp_out = (
            htf_mass_eval * cp_htf_ave * htf_temp_in - q_dot_actual
        ) / max(htf_mass_eval * cp_htf_ave, 1.0e-9)
        htf_p_out = htf_p_in

        m_dot_evap = q_dot_actual / max(h_sat_g_prev - h_sat_f_prev, 1.0e-6)

        # FOR VERIFICATION ONLY - BLOWDOWN IS NOT IN MODEL
        m_dot_blowdown = 0.0

        # Tank calculations
        m_dot_superheat = m_dot_superheat_req
        m_tank_new = max(m_tank_prev + (m_dot_in - m_dot_superheat - m_dot_blowdown) * timestep_seconds, 1.0)

        # aa calculations
        dpdt_aa, dhdt_aa = self._state_derivatives(
            p_tank_prev,
            h_tank_prev,
            m_dot_in,
            m_dot_superheat,
            m_dot_blowdown,
            h_in,
            q_dot_actual,
            vol_tank,
        )

        p_aa = p_tank_prev + dpdt_aa * timestep_seconds / 2.0
        h_aa = h_tank_prev + dhdt_aa * timestep_seconds / 2.0

        # bb calculations
        dpdt_bb, dhdt_bb = self._state_derivatives(
            p_aa,
            h_aa,
            m_dot_in,
            m_dot_superheat,
            m_dot_blowdown,
            h_in,
            q_dot_actual,
            vol_tank,
        )

        p_bb = p_tank_prev + dpdt_bb * timestep_seconds / 2.0
        h_bb = h_tank_prev + dhdt_bb * timestep_seconds / 2.0

        # cc calculations
        dpdt_cc, dhdt_cc = self._state_derivatives(
            p_bb,
            h_bb,
            m_dot_in,
            m_dot_superheat,
            m_dot_blowdown,
            h_in,
            q_dot_actual,
            vol_tank,
        )

        p_cc = p_tank_prev + dpdt_cc * timestep_seconds
        h_cc = h_tank_prev + dhdt_cc * timestep_seconds

        # dd calculations
        dpdt_dd, dhdt_dd = self._state_derivatives(
            p_cc,
            h_cc,
            m_dot_in,
            m_dot_superheat,
            m_dot_blowdown,
            h_in,
            q_dot_actual,
            vol_tank,
        )

        # End of timestep Pressure and Enthlapy
        p_tank_new = p_tank_prev + (dpdt_aa + 2.0 * dpdt_bb + 2.0 * dpdt_cc + dpdt_dd) * timestep_seconds / 6.0
        h_tank_new = h_tank_prev + (dhdt_aa + 2.0 * dhdt_bb + 2.0 * dhdt_cc + dhdt_dd) * timestep_seconds / 6.0

        # Fallback if RK4 finite-difference derivatives become numerically unstable.
        energy_prev = m_tank_prev * h_tank_prev
        energy_in = m_dot_in * h_in * timestep_seconds
        energy_out = (m_dot_superheat * h_sat_g_prev + m_dot_blowdown * h_sat_f_prev) * timestep_seconds
        energy_new = max(energy_prev + energy_in - energy_out + q_dot_actual * timestep_seconds, 1.0)
        h_tank_fallback = energy_new / max(m_tank_new, 1.0)
        mass_ratio = m_tank_new / max(m_tank_prev, 1.0e-9)
        enthalpy_ratio = h_tank_fallback / max(h_tank_prev, 1.0e-9)
        p_tank_fallback = max(p_tank_prev * (0.70 + 0.20 * mass_ratio + 0.10 * enthalpy_ratio), 1.0)

        rk4_invalid = (
            p_tank_new != p_tank_new
            or h_tank_new != h_tank_new
            or p_tank_new < 0.2 * p_tank_prev
            or p_tank_new > 5.0 * p_tank_prev
            or h_tank_new < 0.1 * h_tank_prev
            or h_tank_new > 5.0 * h_tank_prev
        )
        if rk4_invalid:
            p_tank_new = p_tank_fallback
            h_tank_new = h_tank_fallback

        p_tank_new = max(p_tank_new, 1.0)
        h_tank_new = max(h_tank_new, 1.0)

        t_tank_new, h_sat_f_new, h_sat_g_new, rho_f_new, rho_g_new = self._sat_props_from_pressure(p_tank_new)
        x_tank_new = max(min((h_tank_new - h_sat_f_new) / max(h_sat_g_new - h_sat_f_new, 1.0e-9), 1.0), 0.0)
        m_tank_g_new = m_tank_new * x_tank_new
        m_tank_f_new = m_tank_new - m_tank_g_new

        vol_liquid = m_tank_f_new / max(rho_f_new, 1.0e-9)
        l_tank_new = self._tank_level_from_volume(vol_liquid, d_tank, length_tank, l_tank_prev)

        h_out = h_sat_g_new
        vol_dot_fw = m_dot_superheat / max(rho_g_new, 1.0e-9)
        rho_htf = self._rho_htf(fluid_id, htf_temp_out, htf_p_out)
        vol_dot_htf = htf_mass_out / max(rho_htf, 1.0e-9)

        # Set the Outputs from this Model (#,Value)
        self.m_dot_superheat.v = m_dot_superheat
        self.vol_dot_fw.v = vol_dot_fw
        self.p_superheat.v = p_tank_prev
        self.t_superheat.v = t_tank_prev
        self.h_superheat.v = h_out
        self.htf_mass_out.v = htf_mass_out
        self.vol_dot_htf.v = vol_dot_htf
        self.htf_temp_out.v = htf_temp_out
        self.htf_p_out.v = htf_p_out
        self.m_tank_end.v = m_tank_new
        self.l_tank_end.v = l_tank_new
        self.p_tank_end.v = p_tank_new
        self.t_tank_end.v = t_tank_new
        self.h_tank_end.v = h_tank_new
        self.evaporator_effectiveness.v = eta_od
        self.m_dot_evap.v = m_dot_evap

        # Perform Any "After Convergence" Manipulations at end of each timestep.
        if self.model.is_converged:
            # Save new tank level, pressure and enthalpy as previous values for next timestep.
            self.m_tank_begin.v = self.m_tank_end.v
            self.l_tank_begin.v = self.l_tank_end.v
            self.p_tank_begin.v = self.p_tank_end.v
            self.h_tank_begin.v = self.h_tank_end.v

            # Find Alarms and Trips - Steam Drum Alarms.
            self.low_level_alarm.v = 0.0 if l_tank_new > self._safe(self.low_level_alarm_cond.v, -1.0e12) else 1.0
            if self.low_level_alarm.v == 0.0:
                self.low_level_trip.v = 0.0
            else:
                self.low_level_trip.v = 0.0 if l_tank_new > self._safe(self.low_level_trip_cond.v, -1.0e12) else 1.0

            self.high_level_alarm.v = 0.0 if l_tank_new < self._safe(self.high_level_alarm_cond.v, 1.0e12) else 1.0
            if self.high_level_alarm.v == 0.0:
                self.high_level_trip.v = 0.0
            else:
                self.high_level_trip.v = 0.0 if l_tank_new < self._safe(self.high_level_trip_cond.v, 1.0e12) else 1.0

            self.high_pressure_alarm.v = 0.0 if p_tank_new < self._safe(self.high_pressure_alarm_cond.v, 1.0e12) else 1.0
            if self.high_pressure_alarm.v == 0.0:
                self.high_pressure_trip.v = 0.0
            else:
                self.high_pressure_trip.v = 0.0 if p_tank_new < self._safe(self.high_pressure_trip_cond.v, 1.0e12) else 1.0

            # High Temp Diff between water entering drum and saturation temp.
            t_in = self._temperature_from_ph(p_in, h_in)
            delta_t_fw = abs(t_in - t_tank_new)
            self.high_delta_t_fw_alarm.v = 0.0 if delta_t_fw < self._safe(self.high_delta_t_fw_alarm_cond.v, 1.0e12) else 1.0
            if self.high_delta_t_fw_alarm.v == 0.0:
                self.high_delta_t_fw_trip.v = 0.0
            else:
                self.high_delta_t_fw_trip.v = 0.0 if delta_t_fw < self._safe(self.high_delta_t_fw_trip_cond.v, 1.0e12) else 1.0

            # High HR in Steam Drum.
            tank_hist = list(self._tank_temp_hist)
            hr_tank = self._hr_from_history(tank_hist, t_tank_new, timestep_seconds)
            self.tank_heating_rate.v = hr_tank
            self.high_tank_hr_alarm.v = 0.0 if abs(hr_tank) < self._safe(self.high_tank_hr_alarm_cond.v, 1.0e12) else 1.0
            if self.high_tank_hr_alarm.v == 0.0:
                self.high_tank_hr_trip.v = 0.0
            else:
                self.high_tank_hr_trip.v = 0.0 if abs(hr_tank) < self._safe(self.high_tank_hr_trip_cond.v, 1.0e12) else 1.0

            # Evaporator Alarms and Trips.
            delta_t_htf = htf_temp_in - t_tank_new
            self.high_delta_t_htf_alarm.v = 0.0 if delta_t_htf < self._safe(self.high_delta_t_htf_alarm_cond.v, 1.0e12) else 1.0
            if self.high_delta_t_htf_alarm.v == 0.0:
                self.high_delta_t_htf_trip.v = 0.0
            else:
                self.high_delta_t_htf_trip.v = 0.0 if delta_t_htf < self._safe(self.high_delta_t_htf_trip_cond.v, 1.0e12) else 1.0

            htf_hist = list(self._htf_temp_hist)
            hr_htf = self._hr_from_history(htf_hist, htf_temp_in, timestep_seconds)
            self.htf_inlet_heating_rate.v = hr_htf
            self.high_htf_hr_alarm.v = 0.0 if abs(hr_htf) < self._safe(self.high_htf_hr_alarm_cond.v, 1.0e12) else 1.0
            if self.high_htf_hr_alarm.v == 0.0:
                self.high_htf_hr_trip.v = 0.0
            else:
                self.high_htf_hr_trip.v = 0.0 if abs(hr_htf) < self._safe(self.high_htf_hr_trip_cond.v, 1.0e12) else 1.0

            self._tank_temp_hist.append(t_tank_new)
            self._htf_temp_hist.append(htf_temp_in)
