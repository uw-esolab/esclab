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

    for _idx in range(1, 47):
        locals()[f"parameter_{_idx}"] = Component.Parameter()
    for _idx in range(1, 19):
        locals()[f"input_{_idx}"] = Component.Input()
    for _idx in range(1, 44):
        locals()[f"output_{_idx}"] = Component.Output()

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
        # Pump coefficient blocks: (1-11), (12-22), (23-33)
        p_base = 1 + (p_idx - 1) * 11
        coef_a = self._safe(getattr(self, f"parameter_{p_base}").v, -100.0)
        coef_b = self._safe(getattr(self, f"parameter_{p_base + 1}").v, 100.0)
        coef_c = self._safe(getattr(self, f"parameter_{p_base + 2}").v, 10.0)
        eta_a = self._safe(getattr(self, f"parameter_{p_base + 3}").v, 0.0)
        eta_b = self._safe(getattr(self, f"parameter_{p_base + 4}").v, 0.0)
        eta_c = self._safe(getattr(self, f"parameter_{p_base + 5}").v, 0.0)
        eta_d = self._safe(getattr(self, f"parameter_{p_base + 6}").v, 0.7)

        # NPSH curve for cavitation signal (outputs 37-39)
        npsh_a = self._safe(getattr(self, f"parameter_{p_base + 7}").v, 0.0)
        npsh_b = self._safe(getattr(self, f"parameter_{p_base + 8}").v, 0.0)
        npsh_c = self._safe(getattr(self, f"parameter_{p_base + 9}").v, 0.0)
        npsh_d = self._safe(getattr(self, f"parameter_{p_base + 10}").v, 0.0)

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
        q_prev = max(self._safe(getattr(self, f"output_{17 + p_idx}").v, 1.0e-6) / max(rho_f, 1.0e-9), 1.0e-9)

        out_base = {1: 28, 2: 31, 3: 34}[p_idx]
        point_1x = self._safe(getattr(self, f"output_{out_base}").v, q_prev)
        point_1y = self._safe(getattr(self, f"output_{out_base + 1}").v, 0.0)
        point_2x = self._safe(getattr(self, f"output_{out_base + 2}").v, q_prev)

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
        d_tank = max(self._safe(self.parameter_34.v, 1.0), 1.0e-3)
        length_tank = max(self._safe(self.parameter_35.v, 1.0), 1.0e-3)
        length_tank2pump = max(self._safe(self.parameter_36.v, 0.0), 0.0)
        p_tank_ini = max(self._safe(self.parameter_37.v, 101325.0), 1.0)
        l_tank_ini = max(self._safe(self.parameter_38.v, 0.0), 0.0)

        vent_param = max(self._safe(self.parameter_39.v, 0.0), 0.0)
        m_dot_lpb1_max = max(self._safe(self.parameter_40.v, 0.0), 0.0)
        da_ss_lpb1 = max(self._safe(self.parameter_41.v, 0.0), 0.0)
        extraction_tol = max(self._safe(self.parameter_42.v, 0.0), 0.0)
        ll_alarm_th = self._safe(self.parameter_43.v, -1.0)
        ll_trip_th = self._safe(self.parameter_44.v, -1.0)
        hl_alarm_th = self._safe(self.parameter_45.v, 1.0e9)
        hl_trip_th = self._safe(self.parameter_46.v, 1.0e9)

        # Inputs
        turbine_on = self._safe(self.input_1.v, 0.0)
        power_p1 = self._safe(self.input_2.v, 0.0)
        power_p2 = self._safe(self.input_3.v, 0.0)
        power_p3 = self._safe(self.input_4.v, 0.0)
        pump_speed = max(self._safe(self.input_5.v, 0.0), 0.0)

        m_dot_fw_in = max(self._safe(self.input_6.v, 0.0), 0.0)
        p_fw_in = max(self._safe(self.input_7.v, 101325.0), 1.0)
        h_fw_in = self._safe(self.input_8.v, 1.0e6)
        p_sd = max(self._safe(self.input_9.v, p_fw_in), 1.0)
        p_piping_sys = max(self._safe(self.input_10.v, p_fw_in), 1.0)

        m_dot_tb = max(self._safe(self.input_11.v, 0.0), 0.0)
        h_tb = self._safe(self.input_13.v, h_fw_in)
        m_dot_hpfwh = max(self._safe(self.input_14.v, 0.0), 0.0)
        h_hpfwh = self._safe(self.input_16.v, h_fw_in)
        h_lpb1 = self._safe(self.input_18.v, h_fw_in)

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
            p_pump_out = self._safe(self.parameter_3.v, 1.0) * max(pump_speed, 0.01) ** 2 * rho_f * 9.81 + p_tank_ini + (l_tank_ini + length_tank2pump) * rho_f * 9.81
            m_dot_pump = 1.0e-4
            q_dot_pump = m_dot_pump / rho_f
            eta_p1 = 0.01
            w_dot_total = max(p_pump_out - (p_tank_ini + (l_tank_ini + length_tank2pump) * rho_f * 9.81), 0.0) * q_dot_pump / eta_p1
            h_pump_out = h_f + w_dot_total / max(m_dot_pump, 1.0e-9)

            self.output_1.v = m_dot_pump
            self.output_2.v = q_dot_pump
            self.output_3.v = p_pump_out
            self.output_4.v = h_pump_out
            self.output_5.v = self._water_temp_from_h(h_pump_out)
            self.output_6.v = m_dot_lpb1
            self.output_7.v = 0.0
            self.output_8.v = 0.0
            self.output_9.v = m_tank
            self.output_10.v = p_tank_ini
            self.output_11.v = l_tank_ini
            self.output_12.v = h_tank
            self.output_13.v = m_tank
            self.output_14.v = t_sat
            self.output_15.v = l_tank_ini
            self.output_16.v = p_tank_ini
            self.output_17.v = h_tank
            self.output_18.v = 1.0e-4
            self.output_19.v = 0.0
            self.output_20.v = 0.0
            self.output_21.v = w_dot_total
            self.output_22.v = w_dot_total
            self.output_23.v = eta_p1
            self.output_24.v = 0.0
            self.output_25.v = 0.0
            self.output_26.v = 0.0
            self.output_27.v = 0.0
            for idx in range(28, 37):
                getattr(self, f"output_{idx}").v = 0.0
            self.output_37.v = 0.0
            self.output_38.v = 0.0
            self.output_39.v = 0.0
            self.output_40.v = 1.0 if l_tank_ini <= ll_alarm_th else 0.0
            self.output_41.v = 1.0 if l_tank_ini <= ll_trip_th else 0.0
            self.output_42.v = 1.0 if l_tank_ini >= hl_alarm_th else 0.0
            self.output_43.v = 1.0 if l_tank_ini >= hl_trip_th else 0.0
            return

        # Current tank states at beginning of timestep/iteration
        m_tank = max(self._safe(self.output_9.v, 1.0), 1.0)
        p_tank = max(self._safe(self.output_10.v, p_tank_ini), 1.0)
        l_tank = max(self._safe(self.output_11.v, l_tank_ini), 0.0)
        h_tank = self._safe(self.output_12.v, h_fw_in)
        _, h_sat_f, h_sat_g, rho_tank_f, _ = self._sat_props_from_pressure(p_tank)

        p_pump_prev = max(self._safe(self.output_3.v, p_tank), 1.0)

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
            h_pump_out = self._safe(self.output_4.v, h_fw_in)
            p_pump_out = self._safe(self.output_3.v, p_tank)
        t_pump_out = self._water_temp_from_h(h_pump_out)
        w_dot_total = p1["w_dot"] + p2["w_dot"] + p3["w_dot"]

        # Step 5: LP bleed request
        if turbine_on == 1.0:
            m_dot_lpb1_prev = max(self._safe(self.output_6.v, 0.0), 0.0)
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
        self.output_1.v = m_dot_pump
        self.output_2.v = vol_dot_pump
        self.output_3.v = p_pump_out
        self.output_4.v = h_pump_out
        self.output_5.v = t_pump_out
        self.output_6.v = m_dot_lpb1
        self.output_7.v = m_dot_vent
        self.output_8.v = h_vent

        self.output_13.v = m_tank_new
        self.output_14.v = t_tank_new
        self.output_15.v = l_tank_new
        self.output_16.v = p_tank_new
        self.output_17.v = h_tank_new

        self.output_18.v = p1["m_dot"]
        self.output_19.v = p2["m_dot"]
        self.output_20.v = p3["m_dot"]
        self.output_21.v = w_dot_total
        self.output_22.v = p1["w_dot"]
        self.output_23.v = p1["eta"]
        self.output_24.v = p2["w_dot"]
        self.output_25.v = p2["eta"]
        self.output_26.v = p3["w_dot"]
        self.output_27.v = p3["eta"]

        self.output_28.v = p1["point_1x"]
        self.output_29.v = p1["point_1y"]
        self.output_30.v = p1["point_2x"]
        self.output_31.v = p2["point_1x"]
        self.output_32.v = p2["point_1y"]
        self.output_33.v = p2["point_2x"]
        self.output_34.v = p3["point_1x"]
        self.output_35.v = p3["point_1y"]
        self.output_36.v = p3["point_2x"]

        # Fortran does cavitation checks in end-of-timestep block; map to converged state.
        self.output_37.v = p1["cav"]
        self.output_38.v = p2["cav"]
        self.output_39.v = p3["cav"]
        self.output_40.v = ll_alarm
        self.output_41.v = ll_trip
        self.output_42.v = hl_alarm
        self.output_43.v = hl_trip

        # Persist beginning-of-timestep states after convergence.
        if self.model.is_converged:
            self.output_9.v = self.output_13.v
            self.output_10.v = self.output_16.v
            self.output_11.v = self.output_15.v
            self.output_12.v = self.output_17.v
