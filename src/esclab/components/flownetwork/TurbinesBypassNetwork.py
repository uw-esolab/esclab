"""Type 6028 turbines and bypass network converted from Fortran."""

import math

from esclab.simulate import Component


class TurbinesBypassNetwork(Component):
    """
    TRNSYS Type 6028: Turbines & Bypass Network.

    Maintains the Fortran I/O topology and staged structure while using a
    stable reduced-order turbine + reheater + piping approximation. This pass
    adds the end-of-timestep alarm/trip channels and heating-rate calculations
    (outputs 96-123) from the original implementation.
    """

    for _idx in range(1, 117):
        locals()[f"parameter_{_idx}"] = Component.Parameter()
    for _idx in range(1, 37):
        locals()[f"input_{_idx}"] = Component.Input()
    for _idx in range(1, 124):
        locals()[f"output_{_idx}"] = Component.Output()

    @staticmethod
    def _safe(value, default=0.0):
        return value if value == value else default

    @staticmethod
    def _sat_temperature_from_pressure(pressure_pa):
        return max(273.15, min(650.0, 373.15 + 42.0 * math.log(max(pressure_pa, 1.0) / 101325.0 + 1.0)))

    @staticmethod
    def _alarm_trip(value, alarm_limit, trip_limit, high=True):
        if high:
            alarm = 1.0 if value >= alarm_limit else 0.0
            trip = 1.0 if value >= trip_limit else 0.0
        else:
            alarm = 1.0 if value <= alarm_limit else 0.0
            trip = 1.0 if value <= trip_limit else 0.0
        return alarm, trip

    @staticmethod
    def _alarm_timed_trip(value, alarm_limit, timed_limit, trip_limit):
        if value < alarm_limit:
            return 0.0, 0.0, 0.0
        if value < timed_limit:
            return 1.0, 0.0, 0.0
        if value < trip_limit:
            return 1.0, 1.0, 0.0
        return 1.0, 1.0, 1.0

    @staticmethod
    def _set_all_outputs_zero(component):
        for idx in range(1, 124):
            getattr(component, f"output_{idx}").v = 0.0

    def _ensure_rate_buffers(self, n_int, htf_t_in, steam_t):
        if not hasattr(self, "_htf_rate_hist") or len(self._htf_rate_hist) != n_int:
            self._htf_rate_hist = [htf_t_in] * n_int
        if not hasattr(self, "_steam_rate_hist") or len(self._steam_rate_hist) != n_int:
            self._steam_rate_hist = [steam_t] * n_int

    @staticmethod
    def _compute_heating_rate(hist, current_t, ts_sec):
        n = len(hist)
        if n <= 0:
            return 0.0
        rate_sum = 0.0
        for idx in range(1, n):
            rate_sum += (hist[idx] - hist[idx - 1]) / ts_sec
        rate_sum += (current_t - hist[-1]) / ts_sec
        return rate_sum / n * 60.0

    @staticmethod
    def _advance_history(hist, new_val):
        if len(hist) > 1:
            hist[:-1] = hist[1:]
        hist[-1] = new_val

    def calculate(self):
        ts_sec = max(self.model.settings.timestep * 3600.0, 1.0e-9)
        n_int = max(int(math.ceil(60.0 / ts_sec)), 1)

        turbine_on = 1.0 if self._safe(self.input_1.v, 0.0) == 1.0 else 0.0
        m_dot_sgt_in = max(self._safe(self.input_2.v, 0.0), 0.0)
        h_sgt_in = self._safe(self.input_3.v, 2.8e6)
        p_sgt_in = max(self._safe(self.input_4.v, 8.0e6), 1.0)

        # Manual valve position inputs
        hpt_cv_vpi = max(min(self._safe(self.input_5.v, 1.0), 1.0), 0.0)
        hp_bypass_vpi = max(min(self._safe(self.input_6.v, 0.0), 1.0), 0.0)
        hp_aux_vpi = max(min(self._safe(self.input_7.v, 0.0), 1.0), 0.0)
        hp_warmup_vpi = max(min(self._safe(self.input_8.v, 0.0), 1.0), 0.0)
        hp_drain_vpi = max(min(self._safe(self.input_9.v, 0.0), 1.0), 0.0)
        aux_da_vpi = max(min(self._safe(self.input_10.v, 0.0), 1.0), 0.0)
        lp_bypass_vpi = max(min(self._safe(self.input_11.v, 0.0), 1.0), 0.0)
        lp_aux_vpi = max(min(self._safe(self.input_12.v, 0.0), 1.0), 0.0)
        lp_warmup_vpi = max(min(self._safe(self.input_13.v, 0.0), 1.0), 0.0)
        lp_drain_vpi = max(min(self._safe(self.input_14.v, 0.0), 1.0), 0.0)

        # Reheater and condenser side inputs
        m_dot_htf = max(self._safe(self.input_23.v, 0.0), 0.0)
        t_htf_in = self._safe(self.input_24.v, 560.0)
        p_htf_in = max(self._safe(self.input_25.v, 2.0e6), 1.0)
        p_cond = max(self._safe(self.input_26.v, 1.0e5), 1.0)

        self._set_all_outputs_zero(self)

        # Outputs 1-11 are control/status passthroughs
        self.output_1.v = turbine_on
        self.output_2.v = hpt_cv_vpi
        self.output_3.v = hp_bypass_vpi
        self.output_4.v = hp_aux_vpi
        self.output_5.v = hp_warmup_vpi
        self.output_6.v = hp_drain_vpi
        self.output_7.v = aux_da_vpi
        self.output_8.v = lp_bypass_vpi
        self.output_9.v = lp_aux_vpi
        self.output_10.v = lp_warmup_vpi
        self.output_11.v = lp_drain_vpi

        # Basic reduced-order turbine split model
        if turbine_on == 1.0:
            m_dot_hpt_in = m_dot_sgt_in * hpt_cv_vpi
            m_dot_hpt_s1 = 0.52 * m_dot_hpt_in
            m_dot_hpt_s2 = m_dot_hpt_in - m_dot_hpt_s1
            m_dot_hpt_exh = m_dot_hpt_s2

            p_hpt1 = max(0.65 * p_sgt_in, p_cond)
            p_hpt2 = max(0.30 * p_sgt_in, p_cond)
            p_hpt_exh = p_hpt2
            h_hpt1 = h_sgt_in - 1.3e5
            h_hpt2 = h_hpt1 - 1.1e5
            h_hpt_exh = h_hpt2
            t_hpt2 = 500.0

            m_dot_lpt_s1 = m_dot_hpt_exh
            m_dot_lpt_s2 = 0.75 * m_dot_lpt_s1
            m_dot_lpt_s3 = 0.82 * m_dot_lpt_s2
            m_dot_lpt_s4 = 0.90 * m_dot_lpt_s3
            m_dot_lpt_exh = m_dot_lpt_s4

            p_lpt1 = max(0.65 * p_hpt_exh, p_cond)
            p_lpt2 = max(0.45 * p_hpt_exh, p_cond)
            p_lpt3 = max(0.28 * p_hpt_exh, p_cond)
            p_lpt4 = max(0.15 * p_hpt_exh, p_cond)
            h_lpt1 = h_hpt_exh - 7.0e4
            h_lpt2 = h_lpt1 - 6.0e4
            h_lpt3 = h_lpt2 - 5.0e4
            h_lpt4 = h_lpt3 - 4.0e4
            h_lpt_exh = h_lpt4
            t_lpt_exh = 350.0

            w_hpt = max((h_sgt_in - h_hpt_exh) * m_dot_hpt_in, 0.0)
            w_lpt = max((h_hpt_exh - h_lpt_exh) * m_dot_lpt_exh, 0.0)
            w_dot_total = w_hpt + w_lpt
        else:
            m_dot_hpt_in = 0.0
            m_dot_hpt_s1 = 0.0
            m_dot_hpt_s2 = 0.0
            m_dot_hpt_exh = 0.0
            m_dot_lpt_s1 = 0.0
            m_dot_lpt_s2 = 0.0
            m_dot_lpt_s3 = 0.0
            m_dot_lpt_s4 = 0.0
            m_dot_lpt_exh = 0.0

            p_hpt1 = p_sgt_in
            p_hpt2 = p_sgt_in
            p_hpt_exh = p_sgt_in
            p_lpt1 = p_cond
            p_lpt2 = p_cond
            p_lpt3 = p_cond
            p_lpt4 = p_cond
            h_hpt1 = h_sgt_in
            h_hpt2 = h_sgt_in
            h_hpt_exh = h_sgt_in
            h_lpt1 = h_sgt_in
            h_lpt2 = h_sgt_in
            h_lpt3 = h_sgt_in
            h_lpt4 = h_sgt_in
            h_lpt_exh = h_sgt_in
            t_hpt2 = 0.0
            t_lpt_exh = 0.0
            w_dot_total = 0.0

        # Bypass + drain flows (outputs 83-91)
        m_dot_hp_bypass = hp_bypass_vpi * 0.2 * m_dot_sgt_in
        m_dot_hp_aux = hp_aux_vpi * 0.1 * m_dot_sgt_in
        m_dot_hp_warmup = hp_warmup_vpi * 0.1 * m_dot_sgt_in
        m_dot_hp_drain = hp_drain_vpi * 0.08 * m_dot_sgt_in
        m_dot_aux_da = aux_da_vpi * 0.08 * m_dot_sgt_in

        m_dot_lp_bypass = lp_bypass_vpi * 0.2 * m_dot_lpt_exh
        m_dot_lp_aux = lp_aux_vpi * 0.1 * m_dot_lpt_exh
        m_dot_lp_warmup = lp_warmup_vpi * 0.1 * m_dot_lpt_exh
        m_dot_lp_drain = lp_drain_vpi * 0.08 * m_dot_lpt_exh

        h_lp_bypass = h_lpt_exh
        h_lp_warmup = h_lpt_exh
        h_hp_warmup = h_hpt_exh

        # Steam separator + reheater surrogate
        m_dot_ss_in = m_dot_hpt_exh + m_dot_hp_bypass
        if m_dot_ss_in > 0.0:
            h_ss_in = (m_dot_hpt_exh * h_hpt_exh + m_dot_hp_bypass * h_sgt_in) / m_dot_ss_in
            if h_ss_in > 2.5e6:
                m_dot_ss_drain = 0.0
                m_dot_ss_steam = m_dot_ss_in
                h_ss_drain = 0.0
                h_ss_steam = h_ss_in
                t_ss = 450.0
            else:
                x = min(max((h_ss_in - 4.0e5) / max(2.5e6 - 4.0e5, 1.0), 0.0), 1.0)
                m_dot_ss_drain = m_dot_ss_in * (1.0 - x)
                m_dot_ss_steam = m_dot_ss_in - m_dot_ss_drain
                h_ss_drain = 4.0e5
                h_ss_steam = 2.5e6
                t_ss = 373.15
        else:
            m_dot_ss_drain = 0.0
            m_dot_ss_steam = 0.0
            h_ss_drain = 0.0
            h_ss_steam = 0.0
            t_ss = 0.0

        if m_dot_ss_steam > 0.0 and m_dot_htf > 0.0:
            cp_htf = 2200.0
            cp_steam = 4200.0
            q_cap_htf = m_dot_htf * cp_htf * max(t_htf_in - t_ss, 0.0)
            q_cap_steam = m_dot_ss_steam * cp_steam * max(t_htf_in - t_ss, 0.0)
            q_dot_hx = 0.7 * min(q_cap_htf, q_cap_steam)
            h_hx_out = h_ss_steam + q_dot_hx / max(m_dot_ss_steam, 1.0e-9)
            t_hx_out = t_ss + q_dot_hx / max(m_dot_ss_steam * cp_steam, 1.0)
            t_htf_out = t_htf_in - q_dot_hx / max(m_dot_htf * cp_htf, 1.0)
            eta_od = min(max(q_dot_hx / max(q_cap_steam, 1.0), 0.0), 1.0)
        else:
            q_dot_hx = 0.0
            h_hx_out = h_ss_steam
            t_hx_out = t_ss
            t_htf_out = t_htf_in
            eta_od = 0.0

        # Pipe states (surrogate from outputs if available)
        p_hpmain = self._safe(self.output_71.v, p_sgt_in)
        t_hpmain = self._safe(self.output_72.v, 520.0 if turbine_on == 1.0 else 0.0)
        p_lpmain = self._safe(self.output_75.v, p_cond)
        t_lpmain = self._safe(self.output_76.v, 360.0 if turbine_on == 1.0 else 0.0)
        p_aux = self._safe(self.output_79.v, p_cond)
        t_aux = self._safe(self.output_80.v, 340.0 if turbine_on == 1.0 else 0.0)

        # Turbine-seal trip logic (output 95)
        p_ts_req = max(self._safe(self.parameter_44.v, 1.0e5), 1.0)
        p_ts_min = 0.7 * p_ts_req
        m_dot_turbine_seals = max(self._safe(self.parameter_43.v, 0.0), 0.0)
        if turbine_on != 1.0:
            if p_aux > p_ts_req:
                trip_ts = 0.0
            else:
                trip_ts = 1.0
        else:
            trip_ts = 0.0

        # Smoothed steam-drum demand (output 12)
        m_dot_sgt_req_raw = m_dot_hp_bypass + m_dot_hp_aux + m_dot_hpt_in
        m_dot_sgt_prev = self._safe(self.output_12.v, m_dot_sgt_req_raw)
        if abs(m_dot_sgt_prev - m_dot_sgt_req_raw) < 0.1:
            m_dot_sgt_req = m_dot_sgt_prev
        elif m_dot_sgt_prev > m_dot_sgt_req_raw:
            m_dot_sgt_req = m_dot_sgt_prev - abs(m_dot_sgt_prev - m_dot_sgt_req_raw) * 0.6
        else:
            m_dot_sgt_req = m_dot_sgt_prev + abs(m_dot_sgt_prev - m_dot_sgt_req_raw) * 0.6

        # Condenser + DA combined outputs
        m_dot_cond = m_dot_lpt_exh + m_dot_lp_bypass + m_dot_lp_warmup + m_dot_hp_warmup
        if m_dot_cond > 0.0:
            h_cond = (
                m_dot_lpt_exh * h_lpt_exh
                + m_dot_lp_bypass * h_lp_bypass
                + m_dot_lp_warmup * h_lp_warmup
                + m_dot_hp_warmup * h_hp_warmup
            ) / m_dot_cond
        else:
            h_cond = 0.0

        m_dot_da = m_dot_aux_da + m_dot_ss_drain
        if m_dot_da > 0.0:
            h_da = (m_dot_aux_da * h_sgt_in + m_dot_ss_drain * h_ss_drain) / m_dot_da
        else:
            h_da = 0.0

        # Main outputs 12-95 mapping
        self.output_12.v = m_dot_sgt_req
        self.output_13.v = m_dot_cond
        self.output_14.v = h_cond
        self.output_15.v = m_dot_da
        self.output_16.v = h_da
        self.output_17.v = w_dot_total

        self.output_18.v = m_dot_htf
        self.output_19.v = m_dot_htf / 1000.0
        self.output_20.v = p_htf_in
        self.output_21.v = t_htf_out

        self.output_22.v = m_dot_hpt_in
        self.output_23.v = p_hpt_exh
        self.output_24.v = h_hpt_exh
        self.output_25.v = m_dot_hpt_in
        self.output_26.v = p_hpt1
        self.output_27.v = h_hpt1
        self.output_28.v = m_dot_hpt_s2
        self.output_29.v = p_hpt2
        self.output_30.v = h_hpt2
        self.output_31.v = m_dot_hpt_exh
        self.output_32.v = p_hpt_exh
        self.output_33.v = h_hpt_exh
        self.output_34.v = t_hpt2

        self.output_35.v = m_dot_ss_drain
        self.output_36.v = m_dot_ss_drain / 1000.0
        self.output_37.v = p_lpmain
        self.output_38.v = h_ss_drain
        self.output_39.v = t_ss

        self.output_40.v = m_dot_ss_steam
        self.output_41.v = m_dot_ss_steam / 1000.0
        self.output_42.v = p_lpmain
        self.output_43.v = h_ss_steam
        self.output_44.v = t_ss

        self.output_45.v = m_dot_ss_steam
        self.output_46.v = m_dot_ss_steam / 1000.0
        self.output_47.v = t_hx_out
        self.output_48.v = p_lpmain
        self.output_49.v = q_dot_hx
        self.output_50.v = eta_od

        self.output_51.v = m_dot_lpt_s1
        self.output_52.v = p_lpt1
        self.output_53.v = 430.0 if turbine_on == 1.0 else 0.0
        self.output_54.v = h_lpt1
        self.output_55.v = m_dot_lpt_s2
        self.output_56.v = p_lpt2
        self.output_57.v = 400.0 if turbine_on == 1.0 else 0.0
        self.output_58.v = h_lpt2
        self.output_59.v = m_dot_lpt_s3
        self.output_60.v = p_lpt3
        self.output_61.v = 380.0 if turbine_on == 1.0 else 0.0
        self.output_62.v = h_lpt3
        self.output_63.v = m_dot_lpt_s4
        self.output_64.v = 360.0 if turbine_on == 1.0 else 0.0
        self.output_65.v = p_lpt4
        self.output_66.v = h_lpt4
        self.output_67.v = m_dot_lpt_exh
        self.output_68.v = m_dot_lpt_exh / 1000.0
        self.output_69.v = t_lpt_exh
        self.output_70.v = h_lpt_exh

        self.output_71.v = p_hpmain
        self.output_72.v = t_hpmain
        self.output_73.v = 100.0 if turbine_on == 1.0 else 0.0
        self.output_74.v = self._safe(self.output_74.v, t_hpmain)

        self.output_75.v = p_lpmain
        self.output_76.v = t_lpmain
        self.output_77.v = 100.0 if turbine_on == 1.0 else 0.0
        self.output_78.v = self._safe(self.output_78.v, t_lpmain)

        self.output_79.v = p_aux
        self.output_80.v = t_aux
        self.output_81.v = 100.0 if turbine_on == 1.0 else 0.0
        self.output_82.v = self._safe(self.output_82.v, t_aux)

        self.output_83.v = m_dot_hp_bypass
        self.output_84.v = m_dot_hp_aux
        self.output_85.v = m_dot_hp_drain
        self.output_86.v = m_dot_hp_warmup
        self.output_87.v = m_dot_lp_aux
        self.output_88.v = m_dot_lp_bypass
        self.output_89.v = m_dot_lp_drain
        self.output_90.v = m_dot_lp_warmup
        self.output_91.v = m_dot_aux_da

        self.output_92.v = self._safe(self.output_92.v, 0.1)
        self.output_93.v = self._safe(self.output_93.v, 0.1)
        self.output_94.v = self._safe(self.output_94.v, 0.1)
        self.output_95.v = trip_ts

        # End-of-timestep alarm/trip channels (96-123)
        if turbine_on != 1.0:
            for idx in range(96, 108):
                getattr(self, f"output_{idx}").v = 0.0
        else:
            hpt_sh_alarm_fl = self._safe(self.parameter_85.v, 50.0)
            hpt_sh_trip_fl = self._safe(self.parameter_86.v, 20.0)
            hpt_sh_alarm_pl = self._safe(self.parameter_87.v, 50.0)
            hpt_sh_trip_pl = self._safe(self.parameter_88.v, 20.0)
            partial_load = self._safe(self.parameter_89.v, 1.0e8)

            hpt_high_t_alarm = self._safe(self.parameter_90.v, 850.0)
            hpt_high_t_timed = self._safe(self.parameter_91.v, 900.0)
            hpt_high_t_trip = self._safe(self.parameter_92.v, 950.0)
            hpt_exh_p_alarm = self._safe(self.parameter_93.v, 3.0e6)
            hpt_exh_p_trip = self._safe(self.parameter_94.v, 4.0e6)

            lpt_sh_alarm = self._safe(self.parameter_95.v, 50.0)
            lpt_sh_trip = self._safe(self.parameter_96.v, 20.0)
            lpt_high_t_alarm = self._safe(self.parameter_97.v, 750.0)
            lpt_high_t_timed = self._safe(self.parameter_98.v, 800.0)
            lpt_high_t_trip = self._safe(self.parameter_99.v, 850.0)

            t_sat_hp = self._sat_temperature_from_pressure(max(self.output_71.v, 1.0))
            superheat_hp = self.output_72.v - t_sat_hp
            if self.output_17.v > partial_load:
                alarm, trip = self._alarm_trip(superheat_hp, hpt_sh_alarm_fl, hpt_sh_trip_fl, high=False)
            else:
                alarm, trip = self._alarm_trip(superheat_hp, hpt_sh_alarm_pl, hpt_sh_trip_pl, high=False)
            self.output_96.v = alarm
            self.output_97.v = trip

            alarm, timed_trip, trip = self._alarm_timed_trip(self.output_72.v, hpt_high_t_alarm, hpt_high_t_timed, hpt_high_t_trip)
            self.output_98.v = alarm
            self.output_99.v = timed_trip
            self.output_100.v = trip

            alarm, trip = self._alarm_trip(self.output_75.v, hpt_exh_p_alarm, hpt_exh_p_trip, high=True)
            self.output_101.v = alarm
            self.output_102.v = trip

            t_sat_lp = self._sat_temperature_from_pressure(max(self.output_75.v, 1.0))
            superheat_lp = self.output_76.v - t_sat_lp
            alarm, trip = self._alarm_trip(superheat_lp, lpt_sh_alarm, lpt_sh_trip, high=False)
            self.output_103.v = alarm
            self.output_104.v = trip

            alarm, timed_trip, trip = self._alarm_timed_trip(self.output_76.v, lpt_high_t_alarm, lpt_high_t_timed, lpt_high_t_trip)
            self.output_105.v = alarm
            self.output_106.v = timed_trip
            self.output_107.v = trip

        # Reheater alarms/trips
        hi_htf_in_alarm = self._safe(self.parameter_103.v, 900.0)
        hi_htf_in_trip = self._safe(self.parameter_104.v, 950.0)
        lo_htf_out_alarm = self._safe(self.parameter_105.v, 450.0)
        lo_htf_out_trip = self._safe(self.parameter_106.v, 400.0)
        hi_htf_flow_alarm = self._safe(self.parameter_107.v, 1.0e9)
        hi_htf_flow_trip = self._safe(self.parameter_108.v, 1.0e9)
        hi_htf_p_alarm = self._safe(self.parameter_109.v, 1.0e9)
        hi_htf_p_trip = self._safe(self.parameter_110.v, 1.0e9)
        hi_dt_alarm = self._safe(self.parameter_111.v, 200.0)
        hi_dt_trip = self._safe(self.parameter_112.v, 250.0)

        alarm, trip = self._alarm_trip(t_htf_in, hi_htf_in_alarm, hi_htf_in_trip, high=True)
        self.output_110.v = alarm
        self.output_111.v = trip

        alarm, trip = self._alarm_trip(t_htf_out, lo_htf_out_alarm, lo_htf_out_trip, high=False)
        self.output_112.v = alarm
        self.output_113.v = trip

        alarm, trip = self._alarm_trip(m_dot_htf, hi_htf_flow_alarm, hi_htf_flow_trip, high=True)
        self.output_114.v = alarm
        self.output_115.v = trip

        alarm, trip = self._alarm_trip(p_htf_in, hi_htf_p_alarm, hi_htf_p_trip, high=True)
        self.output_116.v = alarm
        self.output_117.v = trip

        delta_t = t_htf_in - self.output_44.v
        alarm, trip = self._alarm_trip(delta_t, hi_dt_alarm, hi_dt_trip, high=True)
        self.output_118.v = alarm
        self.output_119.v = trip

        # Heating-rate channels (108,109,120-123)
        self._ensure_rate_buffers(n_int, t_htf_in, self.output_44.v)
        hr_htf = self._compute_heating_rate(self._htf_rate_hist, t_htf_in, ts_sec)
        hr_steam = self._compute_heating_rate(self._steam_rate_hist, self.output_44.v, ts_sec)

        hr_htf_alarm = self._safe(self.parameter_113.v, 999.0)
        hr_htf_trip = self._safe(self.parameter_114.v, 9999.0)
        hr_steam_alarm = self._safe(self.parameter_115.v, 999.0)
        hr_steam_trip = self._safe(self.parameter_116.v, 9999.0)

        self.output_108.v = hr_htf
        self.output_109.v = hr_steam
        self.output_120.v = 1.0 if abs(hr_htf) >= hr_htf_alarm else 0.0
        self.output_121.v = 1.0 if abs(hr_htf) >= hr_htf_trip else 0.0
        self.output_122.v = 1.0 if abs(hr_steam) >= hr_steam_alarm else 0.0
        self.output_123.v = 1.0 if abs(hr_steam) >= hr_steam_trip else 0.0

        self._advance_history(self._htf_rate_hist, t_htf_in)
        self._advance_history(self._steam_rate_hist, self.output_44.v)
