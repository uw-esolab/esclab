"""Type 6011 deaerator-pump converted from Fortran."""

import math

from esclab.simulate import Component


class DeaeratorPump(Component):
    """
    TRNSYS Type 6011: Deaerator-Pump.

    Uses the Fortran step ordering:
    1) solve each pump operating point from pump/system curves,
    2) combine pump outlet conditions,
    3) solve LP extraction request,
    4) update deaerator tank states and level.

    Notes
    -----
    The original Fortran uses `ESOL6015_myfunctions` derivatives (`drhodhcp`,
    `dudpch`, etc.) for an RK4 pressure-enthalpy state solve. This conversion
    preserves the same I/O and control structure, but replaces that block with a
    stable reduced-order mass/energy update suitable for iterative classroom use.
    """

    p1_head_a = Component.Parameter()
    p1_head_b = Component.Parameter()
    p1_head_c = Component.Parameter()
    p1_eta_a = Component.Parameter()
    p1_eta_b = Component.Parameter()
    p1_eta_c = Component.Parameter()
    p1_eta_d = Component.Parameter()
    p1_npsh_a = Component.Parameter()
    p1_npsh_b = Component.Parameter()
    p1_npsh_c = Component.Parameter()
    p1_npsh_d = Component.Parameter()

    p2_head_a = Component.Parameter()
    p2_head_b = Component.Parameter()
    p2_head_c = Component.Parameter()
    p2_eta_a = Component.Parameter()
    p2_eta_b = Component.Parameter()
    p2_eta_c = Component.Parameter()
    p2_eta_d = Component.Parameter()
    p2_npsh_a = Component.Parameter()
    p2_npsh_b = Component.Parameter()
    p2_npsh_c = Component.Parameter()
    p2_npsh_d = Component.Parameter()

    p3_head_a = Component.Parameter()
    p3_head_b = Component.Parameter()
    p3_head_c = Component.Parameter()
    p3_eta_a = Component.Parameter()
    p3_eta_b = Component.Parameter()
    p3_eta_c = Component.Parameter()
    p3_eta_d = Component.Parameter()
    p3_npsh_a = Component.Parameter()
    p3_npsh_b = Component.Parameter()
    p3_npsh_c = Component.Parameter()
    p3_npsh_d = Component.Parameter()

    tank_diameter = Component.Parameter()
    tank_length = Component.Parameter()
    tank_to_pump_length = Component.Parameter()
    tank_pressure_initial = Component.Parameter()
    tank_level_initial = Component.Parameter()
    vent_setting = Component.Parameter()
    lpb1_flow_max = Component.Parameter()
    lpb1_ramp_limit = Component.Parameter()
    extraction_tolerance = Component.Parameter()
    ll_alarm_threshold = Component.Parameter()
    ll_trip_threshold = Component.Parameter()
    hl_alarm_threshold = Component.Parameter()
    hl_trip_threshold = Component.Parameter()

    in_turbine_on = Component.Input()
    in_pump1_power = Component.Input()
    in_pump2_power = Component.Input()
    in_pump3_power = Component.Input()
    in_pump_speed = Component.Input()
    in_fw_m_dot = Component.Input()
    in_fw_p = Component.Input()
    in_fw_h = Component.Input()
    in_sd_p = Component.Input()
    in_sys_piping_p = Component.Input()
    in_tb_m_dot = Component.Input()
    in_unused_12 = Component.Input()
    in_tb_h = Component.Input()
    in_hpfwh_m_dot = Component.Input()
    in_unused_15 = Component.Input()
    in_hpfwh_h = Component.Input()
    in_unused_17 = Component.Input()
    in_lpb1_h = Component.Input()

    out_pump_m_dot = Component.Output()
    out_pump_vol_dot = Component.Output()
    out_pump_p = Component.Output()
    out_pump_h = Component.Output()
    out_pump_t = Component.Output()
    out_lpb1_m_dot = Component.Output()
    out_vent_m_dot = Component.Output()
    out_vent_h = Component.Output()
    out_state_m_start = Component.Output()
    out_state_p_start = Component.Output()
    out_state_l_start = Component.Output()
    out_state_h_start = Component.Output()
    out_tank_m = Component.Output()
    out_tank_t = Component.Output()
    out_tank_l = Component.Output()
    out_tank_p = Component.Output()
    out_tank_h = Component.Output()
    out_p1_m_dot = Component.Output()
    out_p2_m_dot = Component.Output()
    out_p3_m_dot = Component.Output()
    out_pump_w_total = Component.Output()
    out_p1_w_dot = Component.Output()
    out_p1_eta = Component.Output()
    out_p2_w_dot = Component.Output()
    out_p2_eta = Component.Output()
    out_p3_w_dot = Component.Output()
    out_p3_eta = Component.Output()
    out_p1_point_1x = Component.Output()
    out_p1_point_1y = Component.Output()
    out_p1_point_2x = Component.Output()
    out_p2_point_1x = Component.Output()
    out_p2_point_1y = Component.Output()
    out_p2_point_2x = Component.Output()
    out_p3_point_1x = Component.Output()
    out_p3_point_1y = Component.Output()
    out_p3_point_2x = Component.Output()
    out_p1_cav = Component.Output()
    out_p2_cav = Component.Output()
    out_p3_cav = Component.Output()
    out_ll_alarm = Component.Output()
    out_ll_trip = Component.Output()
    out_hl_alarm = Component.Output()
    out_hl_trip = Component.Output()

    p1_curve_c = p1_head_c

    @staticmethod
    def _safe(value, default=0.0):
        return value if value == value else default

    @staticmethod
    def _sat_temperature_from_pressure(pressure_pa):
        return max(273.15, min(650.0, 373.15 + 42.0 * math.log(max(pressure_pa, 1.0) / 101325.0 + 1.0)))

    @staticmethod
    def _sat_props_from_pressure(pressure_pa):
        t_sat = DeaeratorPump._sat_temperature_from_pressure(pressure_pa)
        h_f = 4200.0 * (t_sat - 273.15)
        h_fg = max(2.5e6 - 1800.0 * (t_sat - 273.15), 5.0e5)
        h_g = h_f + h_fg
        rho_f = max(1000.0 - 0.35 * (t_sat - 273.15), 600.0)
        rho_g = max(pressure_pa / (461.5 * max(t_sat, 200.0)), 0.05)
        return t_sat, h_f, h_g, rho_f, rho_g

    @staticmethod
    def _water_temp_from_h(h):
        return 273.15 + h / 4200.0

    @staticmethod
    def _cylinder_segment_area(radius, level):
        level = min(max(level, 0.0), 2.0 * radius)
        if level <= 0.0:
            return 0.0
        if level >= 2.0 * radius:
            return math.pi * radius**2
        return (
            math.acos((radius - level) / radius) * radius**2
            - (radius - level) * math.sqrt(max(2.0 * radius * level - level**2, 0.0))
        )

    @staticmethod
    def _npshr(q_dot, speed, npsh_a, npsh_b, npsh_c, npsh_d):
        return (
            npsh_a * speed**3 * q_dot**3
            + npsh_b * speed**2 * q_dot**2
            + npsh_c * speed * q_dot
            + npsh_d
        )

    def _pump_step(self, p_idx, power_on, speed_factor, rho_f, p_tank, l_tank, l_tank2pump, p_piping_sys, p_sd, p_pump_prev):
        pump_blocks = {
            1: {
                "head_a": self.p1_head_a,
                "head_b": self.p1_head_b,
                "head_c": self.p1_head_c,
                "eta_a": self.p1_eta_a,
                "eta_b": self.p1_eta_b,
                "eta_c": self.p1_eta_c,
                "eta_d": self.p1_eta_d,
                "npsh_a": self.p1_npsh_a,
                "npsh_b": self.p1_npsh_b,
                "npsh_c": self.p1_npsh_c,
                "npsh_d": self.p1_npsh_d,
                "m_dot_out": self.out_p1_m_dot,
                "point_1x_out": self.out_p1_point_1x,
                "point_1y_out": self.out_p1_point_1y,
                "point_2x_out": self.out_p1_point_2x,
            },
            2: {
                "head_a": self.p2_head_a,
                "head_b": self.p2_head_b,
                "head_c": self.p2_head_c,
                "eta_a": self.p2_eta_a,
                "eta_b": self.p2_eta_b,
                "eta_c": self.p2_eta_c,
                "eta_d": self.p2_eta_d,
                "npsh_a": self.p2_npsh_a,
                "npsh_b": self.p2_npsh_b,
                "npsh_c": self.p2_npsh_c,
                "npsh_d": self.p2_npsh_d,
                "m_dot_out": self.out_p2_m_dot,
                "point_1x_out": self.out_p2_point_1x,
                "point_1y_out": self.out_p2_point_1y,
                "point_2x_out": self.out_p2_point_2x,
            },
            3: {
                "head_a": self.p3_head_a,
                "head_b": self.p3_head_b,
                "head_c": self.p3_head_c,
                "eta_a": self.p3_eta_a,
                "eta_b": self.p3_eta_b,
                "eta_c": self.p3_eta_c,
                "eta_d": self.p3_eta_d,
                "npsh_a": self.p3_npsh_a,
                "npsh_b": self.p3_npsh_b,
                "npsh_c": self.p3_npsh_c,
                "npsh_d": self.p3_npsh_d,
                "m_dot_out": self.out_p3_m_dot,
                "point_1x_out": self.out_p3_point_1x,
                "point_1y_out": self.out_p3_point_1y,
                "point_2x_out": self.out_p3_point_2x,
            },
        }
        pump = pump_blocks[p_idx]

        coef_a = self._safe(pump["head_a"].v, -100.0)
        coef_b = self._safe(pump["head_b"].v, 100.0)
        coef_c = self._safe(pump["head_c"].v, 10.0)
        eta_a = self._safe(pump["eta_a"].v, 0.0)
        eta_b = self._safe(pump["eta_b"].v, 0.0)
        eta_c = self._safe(pump["eta_c"].v, 0.0)
        eta_d = self._safe(pump["eta_d"].v, 0.7)

        npsh_a = self._safe(pump["npsh_a"].v, 0.0)
        npsh_b = self._safe(pump["npsh_b"].v, 0.0)
        npsh_c = self._safe(pump["npsh_c"].v, 0.0)
        npsh_d = self._safe(pump["npsh_d"].v, 0.0)

        if power_on != 1.0:
            return {
                "m_dot": 0.0,
                "p_out": p_tank,
                "h_out": self._sat_props_from_pressure(p_tank)[1],
                "w_dot": 0.0,
                "eta": 0.0,
                "q": 0.0,
                "point_1x": 0.0,
                "point_1y": 0.0,
                "point_2x": 0.0,
                "cav": 0.0,
            }

        speed = max(min(speed_factor, 1.0), 0.01)
        if p_idx == 1:
            b = coef_b * speed
            c = coef_c * speed**2
            eta_a_eff = eta_a / speed**4
            eta_b_eff = eta_b / speed**3
            eta_c_eff = eta_c / speed**2
            eta_d_eff = eta_d / speed
        else:
            b = coef_b
            c = coef_c
            eta_a_eff = eta_a / 100.0
            eta_b_eff = eta_b / 100.0
            eta_c_eff = eta_c / 100.0
            eta_d_eff = eta_d / 100.0

        a = coef_a
        q_prev = max(self._safe(pump["m_dot_out"].v, 1.0e-6) / max(rho_f, 1.0e-9), 1.0e-9)

        point_1x = self._safe(pump["point_1x_out"].v, q_prev)
        point_1y = self._safe(pump["point_1y_out"].v, 0.0)
        point_2x = self._safe(pump["point_2x_out"].v, q_prev)

        discr = max(b * b - 4.0 * a * c, 0.0)
        if abs(a) < 1.0e-12:
            q_max = max(-c / max(b, 1.0e-12), 1.0e-6)
        else:
            q_max = max(
                (-b + math.sqrt(discr)) / (2.0 * a),
                (-b - math.sqrt(discr)) / (2.0 * a),
                1.0e-6,
            )

        p_pump_head = p_pump_prev / rho_f / 9.81 - p_tank / rho_f / 9.81 - l_tank - l_tank2pump
        system_headloss = (p_pump_prev - p_piping_sys) / rho_f / 9.81

        if self.model.timestep_iteration == 0:
            err_prev = p_pump_head - p_sd / rho_f / 9.81 - system_headloss + p_tank / rho_f / 9.81 + l_tank + l_tank2pump
            if abs(err_prev) <= 1.0:
                q_new = q_prev
            elif err_prev > 0.0:
                q_new = min(q_prev + 0.001, q_max)
            else:
                q_new = max(q_prev - 0.001, 1.0e-6)
            point_1x = q_prev
            point_1y = err_prev
            point_2x = q_new
            err_new = err_prev
        else:
            err_new = p_pump_head - p_sd / rho_f / 9.81 - system_headloss + p_tank / rho_f / 9.81 + l_tank + l_tank2pump
            if abs(point_2x - point_1x) <= 1.0e-12:
                slope = 0.0
            else:
                slope = (err_new - point_1y) / (point_2x - point_1x)
            intercept = point_1y - slope * point_1x

            if abs(err_new) <= 1.0:
                q_new = point_2x
            elif abs(slope) > 1.0e-12:
                q_star = -intercept / slope
                if err_new >= 0.0:
                    q_star = max(q_star, point_2x)
                    q_star = min(q_star, q_max)
                    dq = abs(q_star - point_2x)
                    q_new = point_2x + 0.2 * dq
                else:
                    q_star = min(q_star, point_2x)
                    q_star = max(q_star, 1.0e-6)
                    dq = abs(q_star - point_2x)
                    q_new = point_2x - 0.2 * dq
            else:
                if err_new >= 0.0:
                    q_new = min(point_2x + 0.0001, q_max)
                else:
                    q_new = max(point_2x - 0.0001, 1.0e-6)

            point_1x = point_2x
            point_1y = err_new
            point_2x = q_new

        p_pump_in = p_tank + rho_f * 9.81 * (l_tank + l_tank2pump)
        m_dot = max(q_new * rho_f, 1.0e-6)
        p_out = (a * q_new**2 + b * q_new + c) * rho_f * 9.81 + p_pump_in

        eta = max(
            eta_a_eff * q_new**4 + eta_b_eff * q_new**3 + eta_c_eff * q_new**2 + eta_d_eff * q_new,
            0.2 if p_idx == 1 else 0.01,
        )

        _, h_f, _, _, _ = self._sat_props_from_pressure(p_tank)
        h_in = h_f
        w_dot = max(p_out - p_pump_in, 0.0) * q_new / max(eta, 1.0e-3)
        h_out = h_in + w_dot / max(m_dot, 1.0e-9)

        d_pump_inlet = 0.1524
        vel = q_new / (math.pi / 4.0 * d_pump_inlet**2)
        npsha = p_pump_in / 1000.0 / 9.81 + vel**2 / (2.0 * 9.81)
        npshr = self._npshr(q_new, speed, npsh_a, npsh_b, npsh_c, npsh_d)
        cav = 1.0 if npsha <= npshr else 0.0

        return {
            "m_dot": m_dot,
            "p_out": p_out,
            "h_out": h_out,
            "w_dot": w_dot,
            "eta": eta,
            "q": q_new,
            "point_1x": point_1x,
            "point_1y": point_1y,
            "point_2x": point_2x,
            "cav": cav,
        }

    def calculate(self):
        # Parameters
        d_tank = max(self._safe(self.tank_diameter.v, 1.0), 1.0e-3)
        length_tank = max(self._safe(self.tank_length.v, 1.0), 1.0e-3)
        length_tank2pump = max(self._safe(self.tank_to_pump_length.v, 0.0), 0.0)
        p_tank_ini = max(self._safe(self.tank_pressure_initial.v, 101325.0), 1.0)
        l_tank_ini = max(self._safe(self.tank_level_initial.v, 0.0), 0.0)

        vent_param = max(self._safe(self.vent_setting.v, 0.0), 0.0)
        m_dot_lpb1_max = max(self._safe(self.lpb1_flow_max.v, 0.0), 0.0)
        da_ss_lpb1 = max(self._safe(self.lpb1_ramp_limit.v, 0.0), 0.0)
        extraction_tol = max(self._safe(self.extraction_tolerance.v, 0.0), 0.0)
        ll_alarm_th = self._safe(self.ll_alarm_threshold.v, -1.0)
        ll_trip_th = self._safe(self.ll_trip_threshold.v, -1.0)
        hl_alarm_th = self._safe(self.hl_alarm_threshold.v, 1.0e9)
        hl_trip_th = self._safe(self.hl_trip_threshold.v, 1.0e9)

        # Inputs
        turbine_on = self._safe(self.in_turbine_on.v, 0.0)
        power_p1 = self._safe(self.in_pump1_power.v, 0.0)
        power_p2 = self._safe(self.in_pump2_power.v, 0.0)
        power_p3 = self._safe(self.in_pump3_power.v, 0.0)
        pump_speed = max(self._safe(self.in_pump_speed.v, 0.0), 0.0)

        m_dot_fw_in = max(self._safe(self.in_fw_m_dot.v, 0.0), 0.0)
        p_fw_in = max(self._safe(self.in_fw_p.v, 101325.0), 1.0)
        h_fw_in = self._safe(self.in_fw_h.v, 1.0e6)
        p_sd = max(self._safe(self.in_sd_p.v, p_fw_in), 1.0)
        p_piping_sys = max(self._safe(self.in_sys_piping_p.v, p_fw_in), 1.0)

        m_dot_tb = max(self._safe(self.in_tb_m_dot.v, 0.0), 0.0)
        h_tb = self._safe(self.in_tb_h.v, h_fw_in)
        m_dot_hpfwh = max(self._safe(self.in_hpfwh_m_dot.v, 0.0), 0.0)
        h_hpfwh = self._safe(self.in_hpfwh_h.v, h_fw_in)
        h_lpb1 = self._safe(self.in_lpb1_h.v, h_fw_in)

        ts = max(self.model.settings.timestep * 3600.0, 1.0e-9)
        radius = d_tank / 2.0
        area_liq_ini = self._cylinder_segment_area(radius, min(l_tank_ini, d_tank))
        vol_tank = math.pi / 4.0 * d_tank**2 * max(length_tank - d_tank, 0.0) + 4.0 / 3.0 * math.pi * radius**3

        if self.model.is_first_step:
            t_sat, h_f, h_g, rho_f, rho_g = self._sat_props_from_pressure(p_tank_ini)
            m_tank_f = area_liq_ini * length_tank * rho_f
            m_tank_g = max(vol_tank - area_liq_ini * length_tank, 0.0) * rho_g
            m_tank = max(m_tank_f + m_tank_g, 1.0)
            x_tank = min(max(m_tank_g / max(m_tank, 1.0e-9), 0.0), 1.0)
            h_tank = h_f + x_tank * (h_g - h_f)

            m_dot_lpb1 = m_dot_lpb1_max if turbine_on == 1.0 else 0.0
            p_pump_out = self._safe(self.p1_curve_c.v, 1.0) * max(pump_speed, 0.01) ** 2 * rho_f * 9.81 + p_tank_ini + (l_tank_ini + length_tank2pump) * rho_f * 9.81
            m_dot_pump = 1.0e-4
            q_dot_pump = m_dot_pump / rho_f
            eta_p1 = 0.01
            w_dot_total = max(p_pump_out - (p_tank_ini + (l_tank_ini + length_tank2pump) * rho_f * 9.81), 0.0) * q_dot_pump / eta_p1
            h_pump_out = h_f + w_dot_total / max(m_dot_pump, 1.0e-9)

            self.out_pump_m_dot.v = m_dot_pump
            self.out_pump_vol_dot.v = q_dot_pump
            self.out_pump_p.v = p_pump_out
            self.out_pump_h.v = h_pump_out
            self.out_pump_t.v = self._water_temp_from_h(h_pump_out)
            self.out_lpb1_m_dot.v = m_dot_lpb1
            self.out_vent_m_dot.v = 0.0
            self.out_vent_h.v = 0.0
            self.out_state_m_start.v = m_tank
            self.out_state_p_start.v = p_tank_ini
            self.out_state_l_start.v = l_tank_ini
            self.out_state_h_start.v = h_tank
            self.out_tank_m.v = m_tank
            self.out_tank_t.v = t_sat
            self.out_tank_l.v = l_tank_ini
            self.out_tank_p.v = p_tank_ini
            self.out_tank_h.v = h_tank
            self.out_p1_m_dot.v = 1.0e-4
            self.out_p2_m_dot.v = 0.0
            self.out_p3_m_dot.v = 0.0
            self.out_pump_w_total.v = w_dot_total
            self.out_p1_w_dot.v = w_dot_total
            self.out_p1_eta.v = eta_p1
            self.out_p2_w_dot.v = 0.0
            self.out_p2_eta.v = 0.0
            self.out_p3_w_dot.v = 0.0
            self.out_p3_eta.v = 0.0
            self.out_p1_point_1x.v = 0.0
            self.out_p1_point_1y.v = 0.0
            self.out_p1_point_2x.v = 0.0
            self.out_p2_point_1x.v = 0.0
            self.out_p2_point_1y.v = 0.0
            self.out_p2_point_2x.v = 0.0
            self.out_p3_point_1x.v = 0.0
            self.out_p3_point_1y.v = 0.0
            self.out_p3_point_2x.v = 0.0
            self.out_p1_cav.v = 0.0
            self.out_p2_cav.v = 0.0
            self.out_p3_cav.v = 0.0
            self.out_ll_alarm.v = 1.0 if l_tank_ini <= ll_alarm_th else 0.0
            self.out_ll_trip.v = 1.0 if l_tank_ini <= ll_trip_th else 0.0
            self.out_hl_alarm.v = 1.0 if l_tank_ini >= hl_alarm_th else 0.0
            self.out_hl_trip.v = 1.0 if l_tank_ini >= hl_trip_th else 0.0
            return

        # Current tank states at beginning of timestep/iteration
        m_tank = max(self._safe(self.out_state_m_start.v, 1.0), 1.0)
        p_tank = max(self._safe(self.out_state_p_start.v, p_tank_ini), 1.0)
        l_tank = max(self._safe(self.out_state_l_start.v, l_tank_ini), 0.0)
        h_tank = self._safe(self.out_state_h_start.v, h_fw_in)
        _, h_sat_f, h_sat_g, rho_tank_f, _ = self._sat_props_from_pressure(p_tank)

        p_pump_prev = max(self._safe(self.out_pump_p.v, p_tank), 1.0)

        # Step 1-3: Pump curve solves
        p1 = self._pump_step(1, power_p1, pump_speed, rho_tank_f, p_tank, l_tank, length_tank2pump, p_piping_sys, p_sd, p_pump_prev)
        p2 = self._pump_step(2, power_p2, 1.0, rho_tank_f, p_tank, l_tank, length_tank2pump, p_piping_sys, p_sd, p_pump_prev)
        p3 = self._pump_step(3, power_p3, 1.0, rho_tank_f, p_tank, l_tank, length_tank2pump, p_piping_sys, p_sd, p_pump_prev)

        # Step 4: combine pump flows
        m_dot_pump = p1["m_dot"] + p2["m_dot"] + p3["m_dot"]
        if m_dot_pump > 1.0e-9:
            vol_dot_pump = m_dot_pump / max(rho_tank_f, 1.0e-9)
            h_pump_out = (p1["m_dot"] * p1["h_out"] + p2["m_dot"] * p2["h_out"] + p3["m_dot"] * p3["h_out"]) / m_dot_pump
            p_pump_out = (p1["m_dot"] * p1["p_out"] + p2["m_dot"] * p2["p_out"] + p3["m_dot"] * p3["p_out"]) / m_dot_pump
        else:
            m_dot_pump = 1.0e-8
            vol_dot_pump = 0.0
            h_pump_out = self._safe(self.out_pump_h.v, h_fw_in)
            p_pump_out = self._safe(self.out_pump_p.v, p_tank)
        t_pump_out = self._water_temp_from_h(h_pump_out)
        w_dot_total = p1["w_dot"] + p2["w_dot"] + p3["w_dot"]

        # Step 5: LP bleed request
        if turbine_on == 1.0:
            m_dot_lpb1_prev = max(self._safe(self.out_lpb1_m_dot.v, 0.0), 0.0)
            h_mix_target = h_sat_f + 0.05 * (h_sat_g - h_sat_f)
            denom = max(h_mix_target - h_lpb1, -1.0e-6)
            if abs(denom) < 1.0e-6:
                m_dot_lpb1_req = m_dot_lpb1_prev
            else:
                m_dot_lpb1_req = (
                    m_dot_fw_in * (h_fw_in - h_mix_target)
                    + m_dot_tb * (h_tb - h_mix_target)
                    + m_dot_hpfwh * (h_hpfwh - h_mix_target)
                ) / denom
            m_dot_lpb1_req = min(max(m_dot_lpb1_req, 0.0), m_dot_lpb1_max)

            delta = m_dot_lpb1_req - m_dot_lpb1_prev
            if abs(delta) > extraction_tol:
                if delta > 0.0:
                    m_dot_lpb1 = min(m_dot_lpb1_prev + abs(delta) * 0.1, m_dot_lpb1_prev + da_ss_lpb1, m_dot_lpb1_max)
                else:
                    m_dot_lpb1 = max(m_dot_lpb1_prev - abs(delta) * 0.1, m_dot_lpb1_prev - da_ss_lpb1, 0.0)
            else:
                m_dot_lpb1 = m_dot_lpb1_prev
        else:
            m_dot_lpb1 = 0.0

        # Step 6: venting and mixture into deaerator
        m_dot_in = m_dot_fw_in + m_dot_tb + m_dot_hpfwh + m_dot_lpb1
        if m_dot_in > 1.0e-9:
            h_in = (
                m_dot_fw_in * h_fw_in
                + m_dot_tb * h_tb
                + m_dot_hpfwh * h_hpfwh
                + m_dot_lpb1 * h_lpb1
            ) / m_dot_in
        else:
            h_in = h_tank

        if vent_param <= 1.0:
            m_dot_vent = vent_param * m_dot_in
        else:
            m_dot_vent = vent_param
        m_dot_vent = max(min(m_dot_vent, m_dot_in), 0.0)
        h_vent = h_sat_g

        # Step 7-8: reduced-order tank state update and new level
        m_tank_new = max(m_tank + (m_dot_in - m_dot_pump - m_dot_vent) * ts, 1.0)
        e_tank_prev = m_tank * h_tank
        e_in = m_dot_in * h_in * ts
        e_out = (m_dot_pump * h_sat_f + m_dot_vent * h_vent) * ts
        h_tank_new = max((e_tank_prev + e_in - e_out) / m_tank_new, 1.0)

        mass_ratio = m_tank_new / max(m_tank, 1.0e-9)
        enthalpy_ratio = h_tank_new / max(h_tank, 1.0e-9)
        p_tank_new = max(p_tank * (0.7 + 0.2 * mass_ratio + 0.1 * enthalpy_ratio), 1.0)

        t_tank_new, _, _, rho_f_new, _ = self._sat_props_from_pressure(p_tank_new)
        # Invert level with bisection using cylinder-segment area model
        l_low, l_high = 0.0, d_tank
        m_liq = m_tank_new * max(1.0 - min(max((h_tank_new - h_sat_f) / max(h_sat_g - h_sat_f, 1.0), 0.0), 1.0), 0.0)
        vol_liq_target = m_liq / max(rho_f_new, 1.0e-9)
        for _ in range(40):
            l_mid = 0.5 * (l_low + l_high)
            a_mid = self._cylinder_segment_area(radius, l_mid)
            v_mid = a_mid * length_tank
            if v_mid > vol_liq_target:
                l_high = l_mid
            else:
                l_low = l_mid
        l_tank_new = 0.5 * (l_low + l_high)

        ll_alarm = 1.0 if l_tank_new <= ll_alarm_th else 0.0
        ll_trip = 1.0 if l_tank_new <= ll_trip_th else 0.0
        hl_alarm = 1.0 if l_tank_new >= hl_alarm_th else 0.0
        hl_trip = 1.0 if l_tank_new >= hl_trip_th else 0.0

        # Outputs
        self.out_pump_m_dot.v = m_dot_pump
        self.out_pump_vol_dot.v = vol_dot_pump
        self.out_pump_p.v = p_pump_out
        self.out_pump_h.v = h_pump_out
        self.out_pump_t.v = t_pump_out
        self.out_lpb1_m_dot.v = m_dot_lpb1
        self.out_vent_m_dot.v = m_dot_vent
        self.out_vent_h.v = h_vent

        self.out_tank_m.v = m_tank_new
        self.out_tank_t.v = t_tank_new
        self.out_tank_l.v = l_tank_new
        self.out_tank_p.v = p_tank_new
        self.out_tank_h.v = h_tank_new

        self.out_p1_m_dot.v = p1["m_dot"]
        self.out_p2_m_dot.v = p2["m_dot"]
        self.out_p3_m_dot.v = p3["m_dot"]
        self.out_pump_w_total.v = w_dot_total
        self.out_p1_w_dot.v = p1["w_dot"]
        self.out_p1_eta.v = p1["eta"]
        self.out_p2_w_dot.v = p2["w_dot"]
        self.out_p2_eta.v = p2["eta"]
        self.out_p3_w_dot.v = p3["w_dot"]
        self.out_p3_eta.v = p3["eta"]

        self.out_p1_point_1x.v = p1["point_1x"]
        self.out_p1_point_1y.v = p1["point_1y"]
        self.out_p1_point_2x.v = p1["point_2x"]
        self.out_p2_point_1x.v = p2["point_1x"]
        self.out_p2_point_1y.v = p2["point_1y"]
        self.out_p2_point_2x.v = p2["point_2x"]
        self.out_p3_point_1x.v = p3["point_1x"]
        self.out_p3_point_1y.v = p3["point_1y"]
        self.out_p3_point_2x.v = p3["point_2x"]

        # Fortran does cavitation checks in end-of-timestep block; map to converged state.
        self.out_p1_cav.v = p1["cav"]
        self.out_p2_cav.v = p2["cav"]
        self.out_p3_cav.v = p3["cav"]
        self.out_ll_alarm.v = ll_alarm
        self.out_ll_trip.v = ll_trip
        self.out_hl_alarm.v = hl_alarm
        self.out_hl_trip.v = hl_trip

        # Persist beginning-of-timestep states after convergence.
        if self.model.is_converged:
            self.out_state_m_start.v = self.out_tank_m.v
            self.out_state_p_start.v = self.out_tank_p.v
            self.out_state_l_start.v = self.out_tank_l.v
            self.out_state_h_start.v = self.out_tank_h.v
