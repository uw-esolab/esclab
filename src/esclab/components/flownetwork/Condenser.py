"""Type 6007 condenser converted from Fortran."""

import math

from eeslib import fluid_properties as fp
from esclab.components.flownetwork.SimplePipe import FricFactor_IC
from esclab.simulate import Component


class Condenser(Component):
    """
    TRNSYS Type 6007: Condenser.

    Preserves Fortran input/output ordering and maps the 3 major calculation
    blocks:
    1) pump flow/pressure solve,
    2) condenser HX heat transfer,
    3) hotwell tank update with reservoir-flow correction.
    """

    tank_height = Component.Parameter()
    tank_area = Component.Parameter()
    tank_level_initial = Component.Parameter()
    tank_pressure_initial = Component.Parameter()
    pump_curve_a = Component.Parameter()
    pump_curve_b = Component.Parameter()
    pump_curve_c = Component.Parameter()
    reserved_parameter_8 = Component.Parameter()
    pump_eta_a = Component.Parameter()
    pump_eta_b = Component.Parameter()
    pump_eta_c = Component.Parameter()
    pump_eta_d = Component.Parameter()
    hx_ttd = Component.Parameter()
    hx_tube_count = Component.Parameter()
    hx_tube_length = Component.Parameter()
    hx_tube_id = Component.Parameter()
    hx_tube_thickness = Component.Parameter()

    steam_m_dot_in = Component.Input()
    steam_h_in = Component.Input()
    steam_p_in = Component.Input()
    cooling_m_dot_in = Component.Input()
    cooling_h_in = Component.Input()
    cooling_p_in = Component.Input()
    deaerator_pressure = Component.Input()
    system_losses = Component.Input()
    pump_speed_signal = Component.Input()

    out_pump_m_dot = Component.Output()
    out_pump_vol_dot = Component.Output()
    out_pump_t = Component.Output()
    out_pump_h = Component.Output()
    out_pump_p = Component.Output()
    out_cooling_m_dot = Component.Output()
    out_cooling_vol_dot = Component.Output()
    out_cooling_t = Component.Output()
    out_cooling_p = Component.Output()
    out_cooling_h = Component.Output()
    out_tank_level = Component.Output()
    out_tank_p = Component.Output()
    out_tank_h = Component.Output()
    out_tank_t = Component.Output()
    out_tank_m = Component.Output()
    out_state_level_prev = Component.Output()
    out_state_p_prev = Component.Output()
    out_state_h_prev = Component.Output()
    out_state_m_prev = Component.Output()
    out_pump_solver_err = Component.Output()
    out_pump_solver_q_prev = Component.Output()
    out_pump_solver_q_anchor = Component.Output()
    out_pump_eta = Component.Output()
    out_pump_w_dot = Component.Output()
    out_reservoir_m_dot = Component.Output()
    out_aux_h = Component.Output()
    out_cond_q_dot = Component.Output()
    out_ff = Component.Output()

    @staticmethod
    def _safe(value, default=0.0):
        return value if value == value else default

    def _solve_pump_flow(
        self, p_tank, l_tank, p_da, sys_losses, pump_speed, q_prev, point_1x, point_1y
    ):
        rho_f = 1000.0
        coef_a = self._safe(self.pump_curve_a.v, 0.0)
        coef_b = self._safe(self.pump_curve_b.v, 0.0)
        coef_c = self._safe(self.pump_curve_c.v, 1.0)

        a = coef_a
        b = coef_b * max(pump_speed, 1.0e-8)
        c = coef_c * max(pump_speed, 1.0e-8) ** 2
        discr = max(b * b - 4.0 * a * c, 0.0)
        if abs(a) < 1.0e-12:
            q_max = max(1.0e-6, -c / max(b, 1.0e-12))
        else:
            q_max = max(
                (-b + math.sqrt(discr)) / (2.0 * a),
                (-b - math.sqrt(discr)) / (2.0 * a),
                1.0e-6,
            )

        p_pump_in = p_tank + rho_f * 9.81 * l_tank
        p_pump_prev = (a * q_prev * q_prev + b * q_prev + c) * rho_f * 9.81 + p_pump_in
        pump_head = p_pump_prev / rho_f / 9.81 - p_tank / rho_f / 9.81 - l_tank
        system_headloss = (p_pump_prev - sys_losses) / rho_f / 9.81
        point_2y = (
            pump_head
            - p_da / rho_f / 9.81
            - system_headloss
            + p_tank / rho_f / 9.81
            + l_tank
        )

        if point_1x != point_1x:
            point_1x = q_prev
            point_1y = point_2y

        if abs(point_2y) <= 1.0:
            q_new = q_prev
        else:
            if abs(q_prev - point_1x) > 1.0e-12:
                slope = (point_2y - point_1y) / (q_prev - point_1x)
                if abs(slope) > 1.0e-12:
                    q_new = q_prev - point_2y / slope
                else:
                    q_new = q_prev + (0.001 if point_2y >= 0.0 else -0.001)
            else:
                q_new = q_prev + (0.001 if point_2y >= 0.0 else -0.001)
            q_new = min(max(q_new, 1.0e-9), q_max)

        p_pump_out = (a * q_new * q_new + b * q_new + c) * rho_f * 9.81 + p_pump_in
        eta_a = self._safe(self.pump_eta_a.v, 0.0)
        eta_b = self._safe(self.pump_eta_b.v, 0.0)
        eta_c = self._safe(self.pump_eta_c.v, 0.0)
        eta_d = self._safe(self.pump_eta_d.v, 0.5)
        spd = max(pump_speed, 1.0e-8)
        flow_ratio = q_new / spd
        eta_pump = max(
            eta_a * flow_ratio**4 + eta_b * flow_ratio**3 + eta_c * flow_ratio**2 + eta_d * flow_ratio,
            0.01,
        )
        w_dot_pump = max(p_pump_out - p_pump_in, 0.0) * q_new / eta_pump

        point_1x_new = q_prev
        point_1y_new = point_2y
        return q_new, p_pump_in, p_pump_out, eta_pump, w_dot_pump, point_1x_new, point_1y_new

    def _condenser_hx(self, m_dot_cool_in, h_cool_in, p_cool_in, p_tank_prev, ff_guess):
        ttd = self._safe(self.hx_ttd.v, 5.0)
        no_tubes = max(self._safe(self.hx_tube_count.v, 1.0), 1.0)
        length_tubes = max(self._safe(self.hx_tube_length.v, 1.0), 1.0e-9)
        id_tube = max(self._safe(self.hx_tube_id.v, 0.02), 1.0e-6)
        th_tube = max(self._safe(self.hx_tube_thickness.v, 0.001), 0.0)
        od_tube = id_tube + 2.0 * th_tube

        p_sat = min(max(p_tank_prev, 1.0), 2.2e7)
        t_sat = float(fp.temperature("water", P=p_sat, Q=0.0))
        h_sat_f = float(fp.enthalpy("water", P=p_sat, Q=0.0))
        h_sat_g = float(fp.enthalpy("water", P=p_sat, Q=1.0))
        rho_sat_f = max(float(fp.density("water", P=p_sat, Q=0.0)), 1.0e-6)
        rho_sat_g = max(float(fp.density("water", P=p_sat, Q=1.0)), 1.0e-6)
        t_cool_in = float(fp.temperature("water", P=max(p_cool_in, 1.0), h=max(h_cool_in, 1.0)))

        if m_dot_cool_in <= 1.0:
            return (
                m_dot_cool_in,
                m_dot_cool_in / 1000.0,
                t_cool_in,
                p_cool_in,
                h_cool_in,
                0.0,
                max(ff_guess, 1.0e-6),
            )

        t_cool_out_ttd = min(max(t_sat - ttd, t_cool_in), 342.0)
        h_cool_out_ttd = float(fp.enthalpy("water", P=max(p_cool_in, 1.0), T=max(t_cool_out_ttd, 273.15)))

        # Inside convection (cooling water side)
        a_s_inner = math.pi * id_tube * length_tubes * no_tubes
        a_s_outer = math.pi * od_tube * length_tubes * no_tubes
        vel = m_dot_cool_in / no_tubes / (math.pi / 4.0 * id_tube**2) / 1000.0
        mu_water = 0.001
        k_water = 0.6
        re_cw = 1000.0 * vel * id_tube / max(mu_water, 1.0e-12)
        if re_cw > 2300.0:
            ff = FricFactor_IC(0.0, re_cw, max(ff_guess, 0.05))
            pr = mu_water * 4200.0 / k_water
            nu = ((ff / 8.0) * (re_cw - 1000.0) * pr) / max(
                1.0 + 12.7 * (ff / 8.0) ** 0.5 * (pr ** (2.0 / 3.0) - 1.0),
                1.0e-9,
            )
            h_bar_cool = nu * k_water / id_tube
        else:
            ff = max(ff_guess, 0.05)
            nu = 3.67
            h_bar_cool = nu * k_water / id_tube
        r_cool = 1.0 / max(h_bar_cool * a_s_inner, 1.0e-12)

        # Outside condensation coefficient (reduced-order from Incropera form)
        t_cw = 0.5 * (t_cool_in + t_cool_out_ttd)
        mu_l = 3.0e-4
        k_l = 0.65
        cp_l = 4200.0
        h_fg_mod = (h_sat_g - h_sat_f) + 0.68 * cp_l * max(t_sat - t_cw, 1.0)
        h_bar_cond = 0.729 * (
            (
                9.81
                * rho_sat_f
                * max(rho_sat_f - rho_sat_g, 1.0e-6)
                * k_l**3
                * h_fg_mod
            )
            / max(mu_l * max(t_sat - t_cw, 1.0) * od_tube, 1.0e-12)
        ) ** 0.25
        r_cond = 1.0 / max(h_bar_cond * a_s_outer, 1.0e-12)

        r_tot = r_cool + r_cond
        q_dot_r = (t_sat - t_cw) / max(r_tot, 1.0e-12)
        q_dot_cw = m_dot_cool_in * (h_cool_out_ttd - h_cool_in)
        q_dot = max(min(q_dot_r, q_dot_cw), 0.0)

        h_cool_out = h_cool_in + q_dot / max(m_dot_cool_in, 1.0e-12)
        t_cool_out = float(fp.temperature("water", P=max(p_cool_in, 1.0), h=max(h_cool_out, 1.0)))
        p_cool_out = p_cool_in
        vol_cool_out = m_dot_cool_in / 1000.0
        return m_dot_cool_in, vol_cool_out, t_cool_out, p_cool_out, h_cool_out, q_dot, max(ff, 1.0e-6)

    def _tank_step(
        self,
        m_tank_prev,
        p_tank_prev,
        h_tank_prev,
        m_dot_in,
        h_in,
        m_dot_pump,
        h_pump,
        m_dot_res,
        q_dot,
        ts,
        area_tank,
    ):
        p_sat_prev = min(max(p_tank_prev, 1.0), 2.2e7)
        h_sat_f_prev = float(fp.enthalpy("water", P=p_sat_prev, Q=0.0))
        e_prev = m_tank_prev * h_tank_prev
        e_in = (m_dot_in * h_in + m_dot_res * h_sat_f_prev) * ts
        e_out = (m_dot_pump * h_pump + q_dot) * ts

        m_tank_new = max(m_tank_prev + (m_dot_in - m_dot_pump + m_dot_res) * ts, 1.0)
        e_new = max(e_prev + e_in - e_out, 1.0)
        h_tank_new = e_new / m_tank_new

        mass_ratio = m_tank_new / max(m_tank_prev, 1.0e-9)
        enth_ratio = h_tank_new / max(h_tank_prev, 1.0e-9)
        p_tank_new = max(p_tank_prev * (0.72 + 0.20 * mass_ratio + 0.08 * enth_ratio), 1.0)
        t_tank_new = float(fp.temperature("water", P=max(p_tank_new, 1.0), h=max(h_tank_new, 1.0)))
        l_tank_new = max(m_tank_new / (area_tank * 1000.0), 0.0)
        return m_tank_new, p_tank_new, h_tank_new, t_tank_new, l_tank_new

    def _solve_reservoir_flow(
        self,
        l_target,
        m_tank_prev,
        p_tank_prev,
        h_tank_prev,
        m_dot_in,
        h_in,
        m_dot_pump,
        h_pump,
        q_dot,
        ts,
        area_tank,
        vol_tank,
    ):
        rho_f = 1000.0
        res_min = -(m_tank_prev - m_dot_pump * ts) + 0.01
        res_max = vol_tank * rho_f - m_tank_prev
        res_tol = 1.0e-3

        m_res_a = max(min(0.0, res_max), res_min)
        state_a = self._tank_step(
            m_tank_prev,
            p_tank_prev,
            h_tank_prev,
            m_dot_in,
            h_in,
            m_dot_pump,
            h_pump,
            m_res_a,
            q_dot,
            ts,
            area_tank,
        )
        err_a = state_a[4] - l_target

        if abs(err_a) <= res_tol:
            return m_res_a, state_a

        m_res_b = max(min(m_res_a + (10.0 if err_a < 0.0 else -10.0), res_max), res_min)
        best = (m_res_a, state_a, abs(err_a))

        for _ in range(120):
            state_b = self._tank_step(
                m_tank_prev,
                p_tank_prev,
                h_tank_prev,
                m_dot_in,
                h_in,
                m_dot_pump,
                h_pump,
                m_res_b,
                q_dot,
                ts,
                area_tank,
            )
            err_b = state_b[4] - l_target
            if abs(err_b) < best[2]:
                best = (m_res_b, state_b, abs(err_b))
            if abs(err_b) <= res_tol:
                return m_res_b, state_b

            if abs(m_res_b - m_res_a) <= 1.0e-6:
                m_res_c = m_res_b + (-10.0 if err_b > 0.0 else 10.0)
            else:
                slope = (err_b - err_a) / (m_res_b - m_res_a)
                if abs(slope) <= 1.0e-12:
                    m_res_c = m_res_b + (-10.0 if err_b > 0.0 else 10.0)
                else:
                    m_res_c = m_res_b - err_b / slope
            m_res_c = max(min(m_res_c, res_max), res_min)
            m_res_a, err_a = m_res_b, err_b
            m_res_b = m_res_b + 0.2 * (m_res_c - m_res_b)

        return best[0], best[1]

    def calculate(self):
        # Parameters
        height_tank = max(self._safe(self.tank_height.v, 1.0), 1.0e-6)
        area_tank = max(self._safe(self.tank_area.v, 1.0), 1.0e-6)
        l_tank_ini = self._safe(self.tank_level_initial.v, 0.0)
        p_tank_ini = max(self._safe(self.tank_pressure_initial.v, 101325.0), 1.0)

        # Inputs
        m_dot_in = max(self._safe(self.steam_m_dot_in.v, 0.0), 0.0)
        h_in = self._safe(self.steam_h_in.v, 1.0e6)
        p_in = max(self._safe(self.steam_p_in.v, p_tank_ini), 1.0)
        m_dot_cool_in = max(self._safe(self.cooling_m_dot_in.v, 0.0), 0.0)
        h_cool_in = self._safe(self.cooling_h_in.v, 1.0e6)
        p_cool_in = self._safe(self.cooling_p_in.v, p_in)
        p_da = self._safe(self.deaerator_pressure.v, p_in)
        sys_losses = self._safe(self.system_losses.v, p_in)
        pump_speed = max(self._safe(self.pump_speed_signal.v, 0.0), 0.0)

        rho_f = 1000.0
        vol_tank = height_tank * area_tank

        if self.model.is_first_step:
            p_sat_ini = min(max(p_tank_ini, 1.0), 2.2e7)
            t_tank = float(fp.temperature("water", P=p_sat_ini, Q=0.0))
            h_tank = float(fp.enthalpy("water", P=p_sat_ini, Q=0.0))
            m_tank = max(area_tank * max(l_tank_ini, 0.0) * rho_f, 1.0)

            p_pump_out = (
                self._safe(self.pump_curve_c.v, 1.0)
                * max(pump_speed, 1.0e-6) ** 2
                * rho_f
                * 9.81
                + p_tank_ini
                + l_tank_ini * rho_f * 9.81
            )
            self.out_pump_m_dot.v = 1.0e-6
            self.out_pump_vol_dot.v = 0.0
            self.out_pump_t.v = t_tank
            self.out_pump_h.v = h_tank
            self.out_pump_p.v = p_pump_out

            self.out_cooling_m_dot.v = m_dot_cool_in
            self.out_cooling_vol_dot.v = m_dot_cool_in / rho_f
            self.out_cooling_t.v = 300.0
            self.out_cooling_p.v = p_cool_in
            self.out_cooling_h.v = h_cool_in

            self.out_tank_level.v = l_tank_ini
            self.out_tank_p.v = p_tank_ini
            self.out_tank_h.v = h_tank
            self.out_tank_t.v = t_tank
            self.out_tank_m.v = m_tank
            self.out_state_level_prev.v = l_tank_ini
            self.out_state_p_prev.v = p_tank_ini
            self.out_state_h_prev.v = h_tank
            self.out_state_m_prev.v = m_tank
            self.out_pump_solver_err.v = 0.0
            self.out_pump_solver_q_prev.v = 0.0
            self.out_pump_solver_q_anchor.v = 0.0
            self.out_pump_eta.v = 0.0
            self.out_pump_w_dot.v = 0.0
            self.out_reservoir_m_dot.v = 0.0
            self.out_aux_h.v = 0.0
            self.out_cond_q_dot.v = 0.0
            self.out_ff.v = 0.1
            return

        # STEP 1: Pump-flow solve
        p_tank_prev = max(self._safe(self.out_state_p_prev.v, p_tank_ini), 1.0)
        l_tank_prev = max(self._safe(self.out_state_level_prev.v, l_tank_ini), 0.0)
        h_tank_prev = self._safe(
            self.out_state_h_prev.v,
            float(fp.enthalpy("water", P=min(max(p_tank_prev, 1.0), 2.2e7), Q=0.0)),
        )
        m_tank_prev = max(
            self._safe(self.out_state_m_prev.v, area_tank * l_tank_prev * rho_f),
            1.0,
        )
        q_prev = max(
            self._safe(self.out_pump_solver_q_prev.v, self._safe(self.out_pump_m_dot.v, 1.0e-6) / rho_f),
            1.0e-9,
        )
        point_1x = self._safe(self.out_pump_solver_q_anchor.v, q_prev)
        point_1y = self._safe(self.out_pump_solver_err.v, 0.0)

        q_new, _, p_pump_out, eta_pump, w_dot_pump, p1x_new, p1y_new = self._solve_pump_flow(
            p_tank_prev,
            l_tank_prev,
            p_da,
            sys_losses,
            pump_speed,
            q_prev,
            point_1x,
            point_1y,
        )
        m_dot_pump = q_new * rho_f
        h_pump_in = float(fp.enthalpy("water", P=min(max(p_tank_prev, 1.0), 2.2e7), Q=0.0))
        h_pump_out = h_pump_in + w_dot_pump / max(m_dot_pump, 1.0e-9)
        t_pump_out = float(fp.temperature("water", P=max(p_pump_out, 1.0), h=max(h_pump_out, 1.0)))

        # STEP 2: Condenser HX calculations
        ff_guess = self._safe(self.out_ff.v, 0.1)
        m_dot_cool_out, vol_cool_out, t_cool_out, p_cool_out, h_cool_out, q_dot, ff_new = self._condenser_hx(
            m_dot_cool_in,
            h_cool_in,
            p_cool_in,
            p_tank_prev,
            ff_guess,
        )

        # STEP 3: Hotwell tank + reservoir flow solve
        ts = max(self.model.settings.timestep, 1.0e-9)
        t_crit = 0.1

        if ts <= t_crit:
            m_dot_res, tank_state = self._solve_reservoir_flow(
                l_tank_ini,
                m_tank_prev,
                p_tank_prev,
                h_tank_prev,
                m_dot_in,
                h_in,
                m_dot_pump,
                h_pump_in,
                q_dot,
                ts,
                area_tank,
                vol_tank,
            )
            m_tank_new, p_tank_new, h_tank_new, t_tank_new, l_tank_new = tank_state
        else:
            ts_sub_n = max(int(math.ceil(ts / t_crit)), 1)
            ts_sub = ts / ts_sub_n
            m_res_accum = 0.0
            m_curr, p_curr, h_curr = m_tank_prev, p_tank_prev, h_tank_prev
            l_curr = l_tank_prev
            t_curr = float(fp.temperature("water", P=max(p_curr, 1.0), h=max(h_curr, 1.0)))
            for _ in range(ts_sub_n):
                m_dot_res_sub, tank_state = self._solve_reservoir_flow(
                    l_tank_ini,
                    m_curr,
                    p_curr,
                    h_curr,
                    m_dot_in,
                    h_in,
                    m_dot_pump,
                    h_pump_in,
                    q_dot,
                    ts_sub,
                    area_tank,
                    vol_tank,
                )
                m_curr, p_curr, h_curr, t_curr, l_curr = tank_state
                m_res_accum += m_dot_res_sub
            m_dot_res = m_res_accum / ts_sub_n
            m_tank_new, p_tank_new, h_tank_new, t_tank_new, l_tank_new = (
                m_curr,
                p_curr,
                h_curr,
                t_curr,
                l_curr,
            )

        # Outputs mapping
        self.out_pump_m_dot.v = m_dot_pump
        self.out_pump_vol_dot.v = q_new
        self.out_pump_t.v = t_pump_out
        self.out_pump_h.v = h_pump_out
        self.out_pump_p.v = p_pump_out
        self.out_cooling_m_dot.v = m_dot_cool_out
        self.out_cooling_vol_dot.v = vol_cool_out
        self.out_cooling_t.v = t_cool_out
        self.out_cooling_p.v = p_cool_out
        self.out_cooling_h.v = h_cool_out

        self.out_tank_level.v = l_tank_new
        self.out_tank_p.v = p_tank_new
        self.out_tank_h.v = h_tank_new
        self.out_tank_t.v = t_tank_new
        self.out_tank_m.v = m_tank_new

        self.out_pump_solver_err.v = p1y_new
        self.out_pump_solver_q_prev.v = q_new
        self.out_pump_solver_q_anchor.v = p1x_new
        self.out_pump_eta.v = eta_pump
        self.out_pump_w_dot.v = w_dot_pump
        self.out_reservoir_m_dot.v = m_dot_res
        self.out_aux_h.v = self._safe(self.out_aux_h.v, h_cool_out)
        self.out_cond_q_dot.v = q_dot
        self.out_ff.v = ff_new

        if self.model.is_converged:
            self.out_state_level_prev.v = self.out_tank_level.v
            self.out_state_p_prev.v = self.out_tank_p.v
            self.out_state_h_prev.v = self.out_tank_h.v
            self.out_state_m_prev.v = self.out_tank_m.v
