"""Type 4012 expansion system converted from Fortran."""

import math

from esclab.components.esol_properties import Incompressible as Inc
from esclab.simulate import Component
from .Valve import CV_data


class ExpansionSystem(Component):
    """
    TRNSYS Type 4012: ESOL4012-ExpansionSystem.

    Parameters
    ----------
    fluid_id, p_ev, p_of, pc_a, pc_b, pc_c, roughness, n_of, n_ev, h_of, h_ev, h_of_init,
    d_of, d_ev, total_mass, only_exp, valve_speed, t_of_init, t_ev_init, cv_orifice, t_chill_out : float
        Type 4012 configuration, geometry, and control parameters.

    Inputs
    ------
    hv_204_09a, lv_205_14, lv_204_05, lv_204_09b, m_dot_in, t_in, p_in, m_counter : float
        Valve positions and process states as in the original Type 4012 interface.

    Outputs
    -------
    output_1..output_15 : float
        Tank levels, outlet state, branch flow diagnostics, and alarms.

    Notes
    -----
    This implementation maps the Type 4012 control structure, including
    start-time initialization, first-iteration valve slew limits,
    mode-dependent flow behavior, and end-of-timestep tank state updates.
    Hydraulic sub-solvers in the original Fortran are represented by reduced
    pressure/Cv-based relations to retain stable classroom-ready behavior.
    """

    fluid_id = Component.Parameter()
    p_ev = Component.Parameter()
    p_of = Component.Parameter()
    pc_a = Component.Parameter()
    pc_b = Component.Parameter()
    pc_c = Component.Parameter()
    roughness = Component.Parameter()
    n_of = Component.Parameter()
    n_ev = Component.Parameter()
    h_of = Component.Parameter()
    h_ev = Component.Parameter()
    h_of_init = Component.Parameter()
    d_of = Component.Parameter()
    d_ev = Component.Parameter()
    total_mass = Component.Parameter()
    only_exp = Component.Parameter()
    valve_speed = Component.Parameter()
    t_of_init = Component.Parameter()
    t_ev_init = Component.Parameter()
    cv_orifice = Component.Parameter()
    t_chill_out = Component.Parameter()

    hv_204_09a = Component.Input()
    lv_205_14 = Component.Input()
    lv_204_05 = Component.Input()
    lv_204_09b = Component.Input()
    m_dot_in = Component.Input()
    t_in = Component.Input()
    p_in = Component.Input()
    m_counter = Component.Input()

    level_ev = Component.Output()
    level_of = Component.Output()
    m_dot_out = Component.Output()
    t_out = Component.Output()
    p_out = Component.Output()
    m_dot_of_balance = Component.Output()
    m_dot_exp = Component.Output()
    m_dot_recirc = Component.Output()
    m_dot_chiller = Component.Output()
    t_of = Component.Output()
    t_ev = Component.Output()
    alarm_high_high_pressure = Component.Output()
    alarm_reserved_1 = Component.Output()
    alarm_reserved_2 = Component.Output()
    alarm_reserved_3 = Component.Output()

    _t_of = float("nan")
    _t_ev = float("nan")
    _level_of = 0.0
    _level_ev = 0.0
    _mass_of = 0.0
    _p_bot_of = float("nan")
    _p_bot_ev = float("nan")
    _m_count_prev = 0.0
    _alarm_high_high_pressure = 0.0

    _lv_204_05_prev = 0.0
    _lv_204_09b_prev = 0.0
    _hv_204_09a_prev = 0.0
    _lv_205_14_prev = 0.0

    _m1_guess = (40.0, 20.0, 20.0)
    _m2_guess = (32.0, 100.0, 132.0)
    _m3_guess = (100.0, 50.0, 50.0, 100.0)
    _props = Inc()

    @staticmethod
    def _safe(value, default=0.0):
        return value if value == value else default

    @staticmethod
    def _clamp01(value):
        return max(0.0, min(1.0, value))

    def _valve_slew(self, target, previous, valve_speed, timestep_s):
        target = self._clamp01(target)
        previous = self._clamp01(previous)
        max_delta = max(valve_speed, 0.0) / 90.0 * max(timestep_s, 0.0)
        if target > previous:
            return min(target, previous + max_delta)
        return max(target, previous - max_delta)

    def _flow_from_dp_cv(self, rho, delta_p, cv):
        q_gpm = max(cv, 1.0e-6) * math.sqrt(max(delta_p, 0.0) / 6894.76)
        q_m3s = q_gpm / 15850.323140625002
        return max(q_m3s * max(rho, 1.0), 0.0)

    def _tank_level(self, mass, diameter, height, n_tanks, rho):
        area = math.pi * (max(diameter, 1.0e-6) / 2.0) ** 2
        denom = max(area * max(height, 1.0e-6) * max(n_tanks, 1.0), 1.0e-9)
        return max(mass / max(rho, 1.0) / denom, 0.0)

    def _pump_flow_estimate(self, rho, p_bot_of, p_target, pc_a, pc_b, pc_c):
        head_static = max((p_target - p_bot_of) / max(rho * 9.81, 1.0e-6), 0.0)
        a = pc_c
        b = pc_b
        c = pc_a - head_static
        if abs(a) < 1.0e-12:
            q = max(-c / max(b, 1.0e-9), 0.0)
        else:
            disc = b * b - 4.0 * a * c
            if disc < 0.0:
                q = 0.0
            else:
                q1 = (-b + math.sqrt(disc)) / (2.0 * a)
                q2 = (-b - math.sqrt(disc)) / (2.0 * a)
                q = max(q1, q2, 0.0)
        return q * rho

    def _initialize_state(self):
        fluid_name = str(self.fluid_id.v)
        rho_of = float(self._props.density(fluid_name, self.t_of_init.v + 273.15, self.p_of.v))
        rho_ev = float(self._props.density(fluid_name, self.t_ev_init.v + 273.15, self.p_ev.v))

        total_mass = self.total_mass.v
        only_exp = self.only_exp.v
        m_counter = self.m_counter.v

        h_of_init = self.h_of_init.v
        h_of = self.h_of.v
        n_of = self.n_of.v
        d_of = self.d_of.v
        h_ev = self.h_ev.v
        n_ev = self.n_ev.v
        d_ev = self.d_ev.v

        if only_exp == 1.0:
            self._mass_of = 0.0
            self._level_of = 0.0
        else:
            of_vol = (math.pi * (max(d_of, 1.0e-6) / 2.0) ** 2) * max(h_of_init, 0.0) * max(h_of, 0.0) * max(n_of, 0.0)
            self._mass_of = max(of_vol * rho_of, 0.0)
            self._level_of = self._tank_level(self._mass_of, d_of, h_of, n_of, rho_of)

        mass_ev = total_mass - self._mass_of - m_counter
        self._level_ev = self._tank_level(mass_ev, d_ev, h_ev, n_ev, rho_ev)

        self._t_of = self.t_of_init.v
        self._t_ev = self.t_ev_init.v
        self._p_bot_of = self.p_of.v + rho_of * 9.81 * self._level_of * h_of
        self._p_bot_ev = self.p_ev.v + rho_ev * 9.81 * self._level_ev * h_ev
        self._m_count_prev = m_counter

        self._lv_204_05_prev = self._clamp01(self.lv_204_05.v)
        self._lv_204_09b_prev = self._clamp01(self.lv_204_09b.v)
        self._hv_204_09a_prev = self._clamp01(self.hv_204_09a.v)
        self._lv_205_14_prev = self._clamp01(self.lv_205_14.v)

        self._alarm_high_high_pressure = 1.0 if self._p_bot_ev > 1.172e6 else 0.0

    def _set_outputs(
        self,
        m_dot_in,
        t_in,
        p_bot_ev,
        m_balance_of,
        m_dot_exp,
        m_dot_recirc,
        m_dot_chiller,
        t_of,
        t_ev,
        alarm,
    ):
        self.level_ev.v = self._level_ev
        self.level_of.v = self._level_of
        self.m_dot_out.v = m_dot_in
        self.t_out.v = t_in
        self.p_out.v = p_bot_ev
        self.m_dot_of_balance.v = m_balance_of
        self.m_dot_exp.v = m_dot_exp
        self.m_dot_recirc.v = m_dot_recirc
        self.m_dot_chiller.v = m_dot_chiller
        self.t_of.v = t_of
        self.t_ev.v = t_ev
        self.alarm_high_high_pressure.v = alarm
        self.alarm_reserved_1.v = 0.0
        self.alarm_reserved_2.v = 0.0
        self.alarm_reserved_3.v = 0.0

    def calculate(self):
        timestep_s = max(self.model.settings.timestep * 3600.0, 1.0e-6)
        fluid_id = self.fluid_id.v
        fluid_name = str(fluid_id)
        p_ev = self.p_ev.v
        p_of = self.p_of.v
        pc_a = self.pc_a.v
        pc_b = self.pc_b.v
        pc_c = self.pc_c.v
        n_of = self.n_of.v
        n_ev = self.n_ev.v
        h_of = self.h_of.v
        h_ev = self.h_ev.v
        d_of = self.d_of.v
        d_ev = self.d_ev.v
        total_mass = self.total_mass.v
        only_exp = self.only_exp.v
        valve_speed = self.valve_speed.v
        cv_orifice = self.cv_orifice.v
        t_chill_out = self.t_chill_out.v

        hv_204_09a = self._clamp01(self.hv_204_09a.v)
        lv_205_14 = self._clamp01(self.lv_205_14.v)
        lv_204_05 = self._clamp01(self.lv_204_05.v)
        lv_204_09b = self._clamp01(self.lv_204_09b.v)

        m_counter = self.m_counter.v
        m_dot_in = self.m_dot_in.v
        t_in = self.t_in.v

        if self.model.is_first_step or self._t_of != self._t_of or self._t_ev != self._t_ev:
            self._initialize_state()

        rho_of = float(self._props.density(fluid_name, self._t_of + 273.15, p_of))
        rho_ev = float(self._props.density(fluid_name, self._t_ev + 273.15, p_ev))
        p_bot_of = self._p_bot_of
        p_bot_ev = self._p_bot_ev

        if self.model.is_first_iteration:
            hv_204_09a = self._valve_slew(hv_204_09a, self._hv_204_09a_prev, valve_speed, timestep_s)
            lv_204_09b = self._valve_slew(lv_204_09b, self._lv_204_09b_prev, valve_speed, timestep_s)
            lv_204_05 = self._valve_slew(lv_204_05, self._lv_204_05_prev, valve_speed, timestep_s)
            lv_205_14 = self._valve_slew(lv_205_14, self._lv_205_14_prev, valve_speed, timestep_s)
            self._hv_204_09a_prev = hv_204_09a
            self._lv_204_09b_prev = lv_204_09b
            self._lv_204_05_prev = lv_204_05
            self._lv_205_14_prev = lv_205_14

        if hv_204_09a > 0.0 and lv_205_14 <= 0.0:
            lv_205_14 = max(lv_205_14, 0.01)

        if hv_204_09a > 0.0 and lv_204_09b > 0.0:
            if hv_204_09a >= lv_204_09b:
                lv_204_09b = 0.0
            else:
                hv_204_09a = 0.0
        if hv_204_09a > 0.0 and lv_204_05 > 0.0:
            if hv_204_09a >= lv_204_05:
                lv_204_05 = 0.0
            else:
                hv_204_09a = 0.0
        if lv_204_09b > 0.0 and lv_204_05 > 0.0:
            if lv_204_09b >= lv_204_05:
                lv_204_05 = 0.0
            else:
                lv_204_09b = 0.0

        cv_hv_204_09 = max(CV_data(1, 1.0, hv_204_09a), 1.0e-6)
        cv_lv_204_09 = max(CV_data(1, 1.0, lv_204_09b), 1.0e-6)
        cv_lv_204_05 = max(CV_data(1, 1.0, lv_204_05), 1.0e-6)
        cv_lv_205_14 = max(CV_data(1, 1.0, lv_205_14), 1.0e-6)

        if self.model.is_first_step:
            self._set_outputs(
                m_dot_in=m_dot_in,
                t_in=t_in,
                p_bot_ev=self._p_bot_ev,
                m_balance_of=0.0,
                m_dot_exp=0.0,
                m_dot_recirc=0.0,
                m_dot_chiller=0.0,
                t_of=self._t_of,
                t_ev=self._t_ev,
                alarm=self._alarm_high_high_pressure,
            )
            return

        if only_exp == 1.0:
            m_dot_sf_to_exp = 0.0
            if self.model.current_time * 3600.0 > timestep_s:
                m_dot_sf_to_exp = (self._m_count_prev - m_counter) / timestep_s

            mass_ev = max(total_mass - m_counter, 1.0e-6)

            if self.model.is_converged:
                cp_ev = max(float(self._props.specheat(fluid_name, self._t_ev + 273.15, p_ev)) * 1000.0, 1.0)
                h_ev_val = cp_ev * (self._t_ev + 273.15)
                t_plant = t_in if m_dot_sf_to_exp > 0.0 else self._t_ev
                cp_plant = max(float(self._props.specheat(fluid_name, t_plant + 273.15, p_ev)) * 1000.0, 1.0)
                h_plant = cp_plant * (t_plant + 273.15)
                d_t_dt = (-h_ev_val * m_dot_sf_to_exp + m_dot_sf_to_exp * h_plant) / max(mass_ev * cp_ev, 1.0)
                self._t_ev = self._t_ev + d_t_dt * timestep_s
                self._m_count_prev = m_counter

                rho_ev = max(float(self._props.density(fluid_name, self._t_ev + 273.15, p_ev)), 1.0)
                self._level_ev = self._tank_level(mass_ev, d_ev, h_ev, n_ev, rho_ev)
                self._p_bot_ev = p_ev + rho_ev * 9.81 * self._level_ev * max(h_ev, 0.0)
                self._alarm_high_high_pressure = 1.0 if self._p_bot_ev > 1.172e6 else 0.0

            self._set_outputs(
                m_dot_in=m_dot_in,
                t_in=t_in,
                p_bot_ev=self._p_bot_ev,
                m_balance_of=0.0,
                m_dot_exp=0.0,
                m_dot_recirc=0.0,
                m_dot_chiller=0.0,
                t_of=self._t_of,
                t_ev=self._t_ev,
                alarm=self._alarm_high_high_pressure,
            )
            return

        m_in_of = 0.0
        m_out_of = 0.0
        m_dot_exp = 0.0
        m_dot_recirc = 0.0
        m_dot_chiller = 0.0
        t_to_exp = self._t_of
        t_to_of = self._t_of

        pump_m_dot = self._pump_flow_estimate(rho_of, p_bot_of, p_bot_ev, pc_a, pc_b, pc_c)
        dp_elevation = max(h_of, 0.0) * max(rho_of, 1.0) * 9.81

        if lv_204_09b > 0.0:
            delta_p_exp = max(p_bot_of + max(pc_a, 0.0) * rho_of * 9.81 - p_bot_ev, 0.0)
            m_dot_exp = self._flow_from_dp_cv(rho_of, delta_p_exp, cv_lv_204_09)
            delta_p_recirc = max(pump_m_dot / max(rho_of, 1.0) ** 2 + dp_elevation, 0.0)
            m_dot_recirc = self._flow_from_dp_cv(rho_of, delta_p_recirc, cv_lv_205_14) * max(lv_205_14, 0.05)
            m_dot_recirc *= cv_orifice / (cv_orifice + 50.0)
            m_out_of = max(m_dot_exp + m_dot_recirc, 0.0)
            m_in_of = m_dot_recirc
            t_to_exp = self._t_of
            t_to_of = self._t_of

        elif hv_204_09a > 0.0:
            delta_p_from_ev = max(p_bot_ev - p_bot_of + dp_elevation, 0.0)
            m_from_ev = self._flow_from_dp_cv(rho_ev, delta_p_from_ev, cv_hv_204_09)
            m_dot_exp = -m_from_ev
            m_dot_recirc = max(pump_m_dot * max(lv_205_14, 0.1), 0.0)
            m_out_of = m_dot_recirc
            m_in_of = m_dot_recirc + m_from_ev
            if m_in_of > 1.0e-9:
                t_to_of = (m_dot_recirc * self._t_of + m_from_ev * self._t_ev) / m_in_of
            else:
                t_to_of = self._t_of
            t_to_exp = self._t_of

        else:
            m_dot_recirc = max(pump_m_dot, 0.0)
            split_chiller = self._clamp01(cv_lv_204_05 / max(cv_lv_204_05 + cv_lv_205_14 + 1.0e-9, 1.0e-9))
            m_dot_chiller = m_dot_recirc * split_chiller
            m_bypass = m_dot_recirc - m_dot_chiller
            t_to_of = (m_dot_chiller * t_chill_out + m_bypass * self._t_of) / max(m_dot_recirc, 1.0e-9)
            m_out_of = m_dot_recirc
            m_in_of = m_dot_recirc
            m_dot_exp = 0.0
            t_to_exp = self._t_of

        if self.model.is_converged:
            m_dot_sf_to_exp = 0.0
            if self.model.current_time * 3600.0 > timestep_s:
                m_dot_sf_to_exp = (self._m_count_prev - m_counter) / timestep_s

            self._mass_of = max(self._mass_of + (m_in_of - m_out_of) * timestep_s, 1.0e-6)

            rho_of_now = max(float(self._props.density(fluid_name, self._t_of + 273.15, p_of)), 1.0)
            self._level_of = self._tank_level(self._mass_of, d_of, h_of, n_of, rho_of_now)

            mass_ev = max(total_mass - self._mass_of - m_counter, 1.0e-6)
            rho_ev_now = max(float(self._props.density(fluid_name, self._t_ev + 273.15, p_ev)), 1.0)
            self._level_ev = self._tank_level(mass_ev, d_ev, h_ev, n_ev, rho_ev_now)

            cp_of = max(float(self._props.specheat(fluid_name, self._t_of + 273.15, p_of)) * 1000.0, 1.0)
            h_of_now = cp_of * (self._t_of + 273.15)
            cp_in_of = max(float(self._props.specheat(fluid_name, t_to_of + 273.15, p_of)) * 1000.0, 1.0)
            h_in_of = cp_in_of * (t_to_of + 273.15)
            d_t_of_dt = (m_in_of * (h_in_of - h_of_now)) / max(self._mass_of * cp_of, 1.0)
            self._t_of = self._t_of + d_t_of_dt * timestep_s

            if m_dot_exp > 0.0:
                m_in_ev = m_dot_exp
                m_out_ev = 0.0
                cp_in_ev = max(float(self._props.specheat(fluid_name, t_to_exp + 273.15, p_ev)) * 1000.0, 1.0)
                h_in_ev = cp_in_ev * (t_to_exp + 273.15)
                h_out_ev = 0.0
            else:
                m_in_ev = 0.0
                m_out_ev = -m_dot_exp
                h_in_ev = 0.0
                cp_out_ev = max(float(self._props.specheat(fluid_name, self._t_ev + 273.15, p_ev)) * 1000.0, 1.0)
                h_out_ev = cp_out_ev * (self._t_ev + 273.15)

            cp_ev = max(float(self._props.specheat(fluid_name, self._t_ev + 273.15, p_ev)) * 1000.0, 1.0)
            cp_plant = max(float(self._props.specheat(fluid_name, t_in + 273.15, p_ev)) * 1000.0, 1.0)
            h_plant = cp_plant * (t_in + 273.15)
            h_ev_now = cp_ev * (self._t_ev + 273.15)
            d_t_ev_dt = (
                -h_ev_now * (m_dot_sf_to_exp + m_in_ev - m_out_ev)
                + m_dot_sf_to_exp * h_plant
                + m_in_ev * h_in_ev
                - m_out_ev * h_out_ev
            ) / max(mass_ev * cp_ev, 1.0)
            self._t_ev = self._t_ev + d_t_ev_dt * timestep_s

            rho_of_new = max(float(self._props.density(fluid_name, self._t_of + 273.15, p_of)), 1.0)
            rho_ev_new = max(float(self._props.density(fluid_name, self._t_ev + 273.15, p_ev)), 1.0)
            self._p_bot_of = p_of + self._level_of * max(h_of, 0.0) * rho_of_new * 9.81
            self._p_bot_ev = p_ev + self._level_ev * max(h_ev, 0.0) * rho_ev_new * 9.81
            self._m_count_prev = m_counter
            self._alarm_high_high_pressure = 1.0 if self._p_bot_ev > 1.172e6 else 0.0

        self._set_outputs(
            m_dot_in=m_dot_in,
            t_in=t_in,
            p_bot_ev=self._p_bot_ev,
            m_balance_of=m_in_of - m_out_of,
            m_dot_exp=m_dot_exp,
            m_dot_recirc=m_dot_recirc,
            m_dot_chiller=m_dot_chiller,
            t_of=self._t_of,
            t_ev=self._t_ev,
            alarm=self._alarm_high_high_pressure,
        )
