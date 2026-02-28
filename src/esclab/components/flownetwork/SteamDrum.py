"""Type 6003 steam drum converted from Fortran."""

import math

from esclab.simulate import Component


class SteamDrum(Component):
    """
    TRNSYS Type 6003: Steam Drum.

    Parameters
    ----------
    parameter_1..parameter_28 : float
        Ordered TRNSYS parameters from Type6003.

    Inputs
    ------
    input_1..input_7 : float
        Ordered TRNSYS inputs from Type6003.

    Outputs
    -------
    output_1..output_36 : float
        Ordered TRNSYS outputs from Type6003.

    Notes
    -----
    This implementation follows the Fortran structure and equations where
    possible with available esclab/eeslib tooling. Correlation-heavy calls
    from `ESOL6015_myfunctions` (e.g., `drhodhcp`, `tank_level`) are replaced
    with numerically stable approximations.
    """

    for _idx in range(1, 29):
        locals()[f"parameter_{_idx}"] = Component.Parameter()
    for _idx in range(1, 8):
        locals()[f"input_{_idx}"] = Component.Input()
    for _idx in range(1, 37):
        locals()[f"output_{_idx}"] = Component.Output()

    _t_tank_prev = float("nan")
    _htf_t_prev = float("nan")

    @staticmethod
    def _sat_props_from_pressure(pressure_pa):
        # Approximate saturated water properties over typical PB operating range.
        pressure_ratio = max(pressure_pa, 1.0) / 101325.0
        t_sat = max(273.15, min(650.0, 373.15 + 42.0 * math.log(pressure_ratio + 1.0)))
        h_f = 4200.0 * (t_sat - 273.15)
        h_fg = max(2.5e6 - 1800.0 * (t_sat - 273.15), 5.0e5)
        h_g = h_f + h_fg
        rho_f = max(1000.0 - 0.35 * (t_sat - 273.15), 600.0)
        rho_g = max(pressure_pa / (461.5 * max(t_sat, 200.0)), 0.05)
        return t_sat, h_f, h_g, rho_f, rho_g

    @staticmethod
    def _cp_htf(_fluid_id, _t, _p):
        # TRNSYS uses specheat() from custom modules. Use a robust default.
        return 2200.0

    @staticmethod
    def _safe_pow(base, exp):
        if base <= 0.0:
            return 0.0
        return base**exp

    @staticmethod
    def _sat_temperature_from_pressure(pressure_pa):
        t_sat = 373.15 + 42.0 * math.log(max(pressure_pa, 1.0) / 101325.0 + 1.0)
        return max(273.15, min(650.0, t_sat))

    @staticmethod
    def _safe(value, default=0.0):
        return value if value == value else default

    def calculate(self):
        # -------------------------------------------------------------------------------------------------------
        # Read parameters and inputs (Fortran ordering)
        # -------------------------------------------------------------------------------------------------------
        d_tank = max(self._safe(self.parameter_1.v, 1.0), 1.0e-3)
        length_tank = max(self._safe(self.parameter_2.v, 1.0), 1.0e-3)
        rated_heat_transfer = max(self._safe(self.parameter_5.v, 0.0), 0.0)
        rated_htf_flow = max(self._safe(self.parameter_6.v, 1.0), 1.0e-9)
        rated_exp = self._safe(self.parameter_7.v, 0.8)
        no_shell_passes = max(self._safe(self.parameter_8.v, 1.0), 1.0)
        no_tube_passes = max(self._safe(self.parameter_9.v, 1.0), 1.0)
        length_hx = max(self._safe(self.parameter_10.v, 1.0), 1.0e-6)
        tube_od = max(self._safe(self.parameter_11.v, 0.01), 1.0e-6)
        tube_th = max(self._safe(self.parameter_12.v, 0.001), 0.0)
        no_tubes = max(self._safe(self.parameter_13.v, 1.0), 1.0)
        fluid_id = self._safe(self.parameter_14.v, 0.0)

        m_dot_in = max(self._safe(self.input_1.v), 0.0)
        h_in = self._safe(self.input_2.v)
        p_in = max(self._safe(self.input_3.v, 101325.0), 1.0)
        htf_mass_in = max(self._safe(self.input_4.v), 0.0)
        htf_temp_in = self._safe(self.input_5.v, 300.0)
        htf_p_in = self._safe(self.input_6.v, p_in)
        m_dot_superheat_req = max(self._safe(self.input_7.v, 60.0), 0.0)

        area = math.pi * (d_tank**2) / 4.0
        vol_tank = area * length_tank
        rho_f = 1000.0

        if self.model.is_first_step:
            l_tank = max(min(self._safe(self.parameter_3.v, 0.5 * d_tank), d_tank), 0.0)
            p_tank = max(self._safe(self.parameter_4.v, p_in), 1.0)
            t_tank, h_f, h_out, rho_f_sat, rho_g = self._sat_props_from_pressure(p_tank)
            rho_f = rho_f_sat

            frac_liquid = max(min(l_tank / d_tank, 1.0), 0.0)
            vol_liquid = vol_tank * frac_liquid
            vol_vapor = max(vol_tank - vol_liquid, 0.0)

            m_tank_f = vol_liquid * rho_f
            m_tank_g = vol_vapor * rho_g
            m_tank_tot = max(m_tank_f + m_tank_g, 1.0)
            x_tank = m_tank_g / m_tank_tot

            h_fg = h_out - h_f
            h_tank = h_f + x_tank * h_fg

            # inital guess for steam flow leaving
            m_dot_superheat = max(m_dot_superheat_req, 60.0)
            vol_dot_fw = m_dot_superheat / rho_g
            rho_htf = 1000.0
            vol_dot_htf = htf_mass_in / rho_htf

            self.output_1.v = m_dot_superheat
            self.output_2.v = vol_dot_fw
            self.output_3.v = p_tank
            self.output_4.v = t_tank
            self.output_5.v = h_out
            self.output_6.v = htf_mass_in
            self.output_7.v = vol_dot_htf
            self.output_8.v = htf_temp_in
            self.output_9.v = htf_p_in
            self.output_10.v = m_tank_tot
            self.output_11.v = l_tank
            self.output_12.v = p_tank
            self.output_13.v = t_tank
            self.output_14.v = h_tank
            self.output_15.v = m_tank_tot
            self.output_16.v = l_tank
            self.output_17.v = p_tank
            self.output_18.v = h_tank
            for idx in range(19, 37):
                getattr(self, f"output_{idx}").v = 0.0

            self._t_tank_prev = t_tank
            self._htf_t_prev = htf_temp_in
            return

        # -------------------------------------------------------------------------------------------------------
        # Read previous-step storage outputs
        # -------------------------------------------------------------------------------------------------------
        ts = max(self.model.settings.timestep, 1.0e-9)
        m_tank_prev = max(self._safe(self.output_15.v, 1.0), 1.0)
        p_tank_prev = max(self._safe(self.output_17.v, p_in), 1.0)
        l_tank_prev = max(self._safe(self.output_16.v, 0.0), 0.0)
        h_tank_prev = self._safe(self.output_18.v, h_in)

        t_tank_prev, h_sat_f_prev, h_sat_g_prev, rho_f_prev, rho_g_prev = self._sat_props_from_pressure(p_tank_prev)

        # -------------------------------------------------------------------------------------------------------
        # Evaporator calculations (mapped from Fortran structure)
        # -------------------------------------------------------------------------------------------------------
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
            eta_1pass = min(2.0 * (1.0 + (1.0 + exp_term) / max(1.0 - exp_term, 1.0e-12)) ** -1.0, 0.99999999)
        cr = 0.0
        ratio = (1.0 - eta_1pass * cr) / max(1.0 - eta_1pass, 1.0e-12)
        eta_od = (ratio**no_shell_passes - 1.0) / max(ratio**no_shell_passes - cr, 1.0e-12)
        eta_od = max(min(eta_od, 1.0), 0.0)

        q_dot_max = htf_mass_eval * cp_htf_ave * (htf_temp_in - t_tank_prev)
        q_dot_actual = max(eta_od * q_dot_max, 0.0)
        htf_mass_out = htf_mass_in
        htf_temp_out = (htf_mass_eval * cp_htf_ave * htf_temp_in - q_dot_actual) / max(htf_mass_eval * cp_htf_ave, 1.0e-9)
        htf_p_out = htf_p_in

        m_dot_evap = q_dot_actual / max(h_sat_g_prev - h_sat_f_prev, 1.0e-6)
        m_dot_blowdown = 0.0

        # -------------------------------------------------------------------------------------------------------
        # Tank calculations (mass/energy update; reduced-order replacement for drhod* routines)
        # -------------------------------------------------------------------------------------------------------
        m_dot_superheat = max(m_dot_superheat_req, 0.0)
        m_tank_new = max(m_tank_prev + (m_dot_in - m_dot_superheat - m_dot_blowdown) * ts, 1.0)

        # Energy balance on the tank control volume
        e_prev = m_tank_prev * h_tank_prev
        e_in = m_dot_in * h_in * ts
        e_out = (m_dot_superheat * h_sat_g_prev + m_dot_blowdown * h_sat_f_prev) * ts
        e_new = max(e_prev + e_in - e_out + q_dot_actual * ts, 1.0)
        h_tank_new = e_new / m_tank_new

        # Pressure update approximation tied to mass and enthalpy change
        mass_ratio = m_tank_new / max(m_tank_prev, 1.0e-9)
        enth_ratio = h_tank_new / max(h_tank_prev, 1.0e-9)
        p_tank_new = max(p_tank_prev * (0.70 + 0.20 * mass_ratio + 0.10 * enth_ratio), 1.0)

        t_tank_new, h_sat_f_new, h_sat_g_new, rho_f_new, rho_g_new = self._sat_props_from_pressure(p_tank_new)

        # Solve for new quality and liquid level
        x_tank_new = max(min((h_tank_new - h_sat_f_new) / max(h_sat_g_new - h_sat_f_new, 1.0e-9), 1.0), 0.0)
        m_tank_g_new = m_tank_new * x_tank_new
        m_tank_f_new = m_tank_new - m_tank_g_new
        vol_liquid = m_tank_f_new / max(rho_f_new, 1.0e-9)
        l_tank_new = max(min(vol_liquid / max(area, 1.0e-9), d_tank), 0.0)

        h_out = h_sat_g_new
        rho_htf = 1000.0
        vol_dot_htf = htf_mass_out / rho_htf
        vol_dot_fw = m_dot_superheat / max(rho_g_prev, 1.0e-6)

        # -------------------------------------------------------------------------------------------------------
        # Alarms/trips and heating-rate channels
        # -------------------------------------------------------------------------------------------------------
        hr_tank = 0.0 if self._t_tank_prev != self._t_tank_prev else (t_tank_new - self._t_tank_prev) / ts * 60.0
        hr_htf = 0.0 if self._htf_t_prev != self._htf_t_prev else (htf_temp_in - self._htf_t_prev) / ts * 60.0

        low_level_alarm = 1.0 if l_tank_new <= self._safe(self.parameter_15.v, -1.0) else 0.0
        low_level_trip = 1.0 if l_tank_new <= self._safe(self.parameter_16.v, -1.0) else 0.0
        high_level_alarm = 1.0 if l_tank_new >= self._safe(self.parameter_17.v, 1.0e12) else 0.0
        high_level_trip = 1.0 if l_tank_new >= self._safe(self.parameter_18.v, 1.0e12) else 0.0
        high_p_alarm = 1.0 if p_tank_new >= self._safe(self.parameter_19.v, 1.0e12) else 0.0
        high_p_trip = 1.0 if p_tank_new >= self._safe(self.parameter_20.v, 1.0e12) else 0.0

        self.output_1.v = m_dot_superheat
        self.output_2.v = vol_dot_fw
        self.output_3.v = p_tank_prev
        self.output_4.v = t_tank_prev
        self.output_5.v = h_out
        self.output_6.v = htf_mass_out
        self.output_7.v = vol_dot_htf
        self.output_8.v = htf_temp_out
        self.output_9.v = htf_p_out
        self.output_10.v = m_tank_new
        self.output_11.v = l_tank_new
        self.output_12.v = p_tank_new
        self.output_13.v = t_tank_new
        self.output_14.v = h_tank_new
        self.output_19.v = eta_od
        self.output_20.v = m_dot_evap
        self.output_21.v = hr_tank
        self.output_22.v = hr_htf
        self.output_23.v = low_level_alarm
        self.output_24.v = low_level_trip
        self.output_25.v = high_level_alarm
        self.output_26.v = high_level_trip
        self.output_27.v = high_p_alarm
        self.output_28.v = high_p_trip
        self.output_29.v = 0.0
        self.output_30.v = 0.0
        self.output_31.v = 0.0
        self.output_32.v = 0.0
        self.output_33.v = 0.0
        self.output_34.v = 0.0
        self.output_35.v = 0.0
        self.output_36.v = 0.0

        if self.model.is_converged:
            self.output_15.v = self.output_10.v
            self.output_16.v = self.output_11.v
            self.output_17.v = self.output_12.v
            self.output_18.v = self.output_14.v

        self._t_tank_prev = t_tank_new
        self._htf_t_prev = htf_temp_in
