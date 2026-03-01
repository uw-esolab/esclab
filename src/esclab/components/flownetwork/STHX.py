"""Type 6017 steam-to-HTF heat exchanger converted from Fortran."""

import math

from eeslib import fluid_properties as fp
from esclab.components.esol_properties import Incompressible as Inc

from esclab.simulate import Component


class STHX(Component):
    """
    TRNSYS Type 6017: ESOL6017-STHX.

    Object: ESOL6017-STHX
    Simulation Studio Model: ESOL6017-STHX

    Parameters
    ----------
    heat_transfer_rated, m_dot_fw_rated, m_dot_htf_rated, rated_exp,
    no_shell_passes, no_tube_passes, length_hx, tube_od, tube_th, no_tubes,
    high_htf_temp_in_alarm, high_htf_temp_in_trip,
    low_htf_temp_in_alarm, low_htf_temp_in_trip,
    low_htf_temp_out_alarm, low_htf_temp_out_trip,
    high_htf_flow_in_alarm, high_htf_flow_in_trip,
    high_htf_pressure_in_alarm, high_htf_pressure_in_trip,
    low_fw_temp_in_alarm, low_fw_temp_in_trip,
    high_hr_htf_in_alarm, high_hr_htf_in_trip,
    high_hr_fw_out_alarm, high_hr_fw_out_trip : float
        Ordered TRNSYS Type 6017 parameter map.
    fluid_id : str
        HTF fluid name passed to :class:`esclab.components.esol_properties.Incompressible`.

    Inputs
    ------
    m_dot_fw, p_fw, h_fw, m_dot_htf, p_htf, t_htf : float
        Feedwater and HTF inlet states.

    Outputs
    -------
    m_dot_fw_out, vol_dot_fw_out, p_fw_out, h_fw_out, t_fw_out,
    m_dot_htf_out, vol_dot_htf_out, p_htf_out, t_htf_out,
    q_dot_hx, eta_od,
    hr_htf_in, hr_fw_out,
    alarm_high_htf_temp_in, trip_high_htf_temp_in,
    alarm_low_htf_temp_in, trip_low_htf_temp_in,
    alarm_low_htf_temp_out, trip_low_htf_temp_out,
    alarm_high_htf_flow_in, trip_high_htf_flow_in,
    alarm_high_htf_pressure_in, trip_high_htf_pressure_in,
    alarm_low_fw_temp_in, trip_low_fw_temp_in,
    alarm_high_hr_htf_in, trip_high_hr_htf_in,
    alarm_high_hr_fw_out, trip_high_hr_fw_out : float
        HX state plus alarm/trip channels.
    """

    heat_transfer_rated = Component.Parameter()
    m_dot_fw_rated = Component.Parameter()
    m_dot_htf_rated = Component.Parameter()
    rated_exp = Component.Parameter()
    no_shell_passes = Component.Parameter()
    no_tube_passes = Component.Parameter()
    length_hx = Component.Parameter()
    tube_od = Component.Parameter()
    tube_th = Component.Parameter()
    no_tubes = Component.Parameter()
    fluid_id = Component.Parameter()
    high_htf_temp_in_alarm = Component.Parameter()
    high_htf_temp_in_trip = Component.Parameter()
    low_htf_temp_in_alarm = Component.Parameter()
    low_htf_temp_in_trip = Component.Parameter()
    low_htf_temp_out_alarm = Component.Parameter()
    low_htf_temp_out_trip = Component.Parameter()
    high_htf_flow_in_alarm = Component.Parameter()
    high_htf_flow_in_trip = Component.Parameter()
    high_htf_pressure_in_alarm = Component.Parameter()
    high_htf_pressure_in_trip = Component.Parameter()
    low_fw_temp_in_alarm = Component.Parameter()
    low_fw_temp_in_trip = Component.Parameter()
    high_hr_htf_in_alarm = Component.Parameter()
    high_hr_htf_in_trip = Component.Parameter()
    high_hr_fw_out_alarm = Component.Parameter()
    high_hr_fw_out_trip = Component.Parameter()

    m_dot_fw = Component.Input()
    p_fw = Component.Input()
    h_fw = Component.Input()
    m_dot_htf = Component.Input()
    p_htf = Component.Input()
    t_htf = Component.Input()

    m_dot_fw_out = Component.Output()
    vol_dot_fw_out = Component.Output()
    p_fw_out = Component.Output()
    h_fw_out = Component.Output()
    t_fw_out = Component.Output()
    m_dot_htf_out = Component.Output()
    vol_dot_htf_out = Component.Output()
    p_htf_out = Component.Output()
    t_htf_out = Component.Output()
    q_dot_hx = Component.Output()
    eta_od = Component.Output()
    hr_htf_in = Component.Output()
    hr_fw_out = Component.Output()
    alarm_high_htf_temp_in = Component.Output()
    trip_high_htf_temp_in = Component.Output()
    alarm_low_htf_temp_in = Component.Output()
    trip_low_htf_temp_in = Component.Output()
    alarm_low_htf_temp_out = Component.Output()
    trip_low_htf_temp_out = Component.Output()
    alarm_high_htf_flow_in = Component.Output()
    trip_high_htf_flow_in = Component.Output()
    alarm_high_htf_pressure_in = Component.Output()
    trip_high_htf_pressure_in = Component.Output()
    alarm_low_fw_temp_in = Component.Output()
    trip_low_fw_temp_in = Component.Output()
    alarm_high_hr_htf_in = Component.Output()
    trip_high_hr_htf_in = Component.Output()
    alarm_high_hr_fw_out = Component.Output()
    trip_high_hr_fw_out = Component.Output()

    _htf_in_hist = []
    _fw_out_hist = []
    _props = Inc()

    @staticmethod
    def _safe(value, default=0.0):
        return value if value == value else default

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
        if len(hist) > 0:
            hist[-1] = new_val

    def _set_outputs_1_to_11(self, m_dot_fw, vol_dot_fw, p_fw, h_fw_out, t_fw_out, m_dot_htf, vol_dot_htf, p_htf, t_htf_out, q_dot_hx, eta_od):
        self.m_dot_fw_out.v = m_dot_fw
        self.vol_dot_fw_out.v = vol_dot_fw
        self.p_fw_out.v = p_fw
        self.h_fw_out.v = h_fw_out
        self.t_fw_out.v = t_fw_out
        self.m_dot_htf_out.v = m_dot_htf
        self.vol_dot_htf_out.v = vol_dot_htf
        self.p_htf_out.v = p_htf
        self.t_htf_out.v = t_htf_out
        self.q_dot_hx.v = q_dot_hx
        self.eta_od.v = eta_od

    def _set_alarm_outputs_zero(self):
        self.hr_htf_in.v = 0.0
        self.hr_fw_out.v = 0.0
        self.alarm_high_htf_temp_in.v = 0.0
        self.trip_high_htf_temp_in.v = 0.0
        self.alarm_low_htf_temp_in.v = 0.0
        self.trip_low_htf_temp_in.v = 0.0
        self.alarm_low_htf_temp_out.v = 0.0
        self.trip_low_htf_temp_out.v = 0.0
        self.alarm_high_htf_flow_in.v = 0.0
        self.trip_high_htf_flow_in.v = 0.0
        self.alarm_high_htf_pressure_in.v = 0.0
        self.trip_high_htf_pressure_in.v = 0.0
        self.alarm_low_fw_temp_in.v = 0.0
        self.trip_low_fw_temp_in.v = 0.0
        self.alarm_high_hr_htf_in.v = 0.0
        self.trip_high_hr_htf_in.v = 0.0
        self.alarm_high_hr_fw_out.v = 0.0
        self.trip_high_hr_fw_out.v = 0.0

    def _run_end_of_timestep_checks(self, m_dot_htf, p_htf, t_htf, p_fw, h_fw):
        # Perform Any "After Convergence" manipulations that may be required at end of timestep.
        ts_sec = max(self.model.settings.timestep, 1.0e-9)
        n_int = max(int(math.ceil(60.0 / ts_sec)), 1)

        t_htf_out = self._safe(self.t_htf_out.v, t_htf)
        t_fw_eval = float(fp.temperature("water", P=max(p_fw, 1.0), h=max(h_fw, 1.0)))
        t_fw_out = self._safe(self.t_fw_out.v, t_fw_eval)
        t_fw = t_fw_eval

        if len(self._htf_in_hist) != n_int:
            # At beginning of simulation set initial temperature arrays for Heating Rates.
            # Heating Rate for HTF Inlet is stored over the last minute.
            self._htf_in_hist = [t_htf] * n_int
        if len(self._fw_out_hist) != n_int:
            # Heating Rate for FW Outlet is stored over the last minute.
            self._fw_out_hist = [t_fw_out] * n_int

        # !!!!!High or Low Tempreature Alarms!!!!
        # High HTF Temp In Check
        alarm, trip = self._alarm_trip(
            t_htf,
            self._safe(self.high_htf_temp_in_alarm.v),
            self._safe(self.high_htf_temp_in_trip.v),
            high=True,
        )
        self.alarm_high_htf_temp_in.v = alarm
        self.trip_high_htf_temp_in.v = trip

        # Low HTF Temp In Check
        alarm, trip = self._alarm_trip(
            t_htf,
            self._safe(self.low_htf_temp_in_alarm.v),
            self._safe(self.low_htf_temp_in_trip.v),
            high=False,
        )
        self.alarm_low_htf_temp_in.v = alarm
        self.trip_low_htf_temp_in.v = trip

        # Low HTF Temp Out Check
        alarm, trip = self._alarm_trip(
            t_htf_out,
            self._safe(self.low_htf_temp_out_alarm.v),
            self._safe(self.low_htf_temp_out_trip.v),
            high=False,
        )
        self.alarm_low_htf_temp_out.v = alarm
        self.trip_low_htf_temp_out.v = trip

        # High HTF Flow entering HX Check
        alarm, trip = self._alarm_trip(
            m_dot_htf,
            self._safe(self.high_htf_flow_in_alarm.v),
            self._safe(self.high_htf_flow_in_trip.v),
            high=True,
        )
        self.alarm_high_htf_flow_in.v = alarm
        self.trip_high_htf_flow_in.v = trip

        # High HTF Pressure entering HX Check
        alarm, trip = self._alarm_trip(
            p_htf,
            self._safe(self.high_htf_pressure_in_alarm.v),
            self._safe(self.high_htf_pressure_in_trip.v),
            high=True,
        )
        self.alarm_high_htf_pressure_in.v = alarm
        self.trip_high_htf_pressure_in.v = trip

        # Low FW Temperature entering HX Check
        alarm, trip = self._alarm_trip(
            t_fw,
            self._safe(self.low_fw_temp_in_alarm.v),
            self._safe(self.low_fw_temp_in_trip.v),
            high=False,
        )
        self.alarm_low_fw_temp_in.v = alarm
        self.trip_low_fw_temp_in.v = trip

        # !!!!!!High Heating Rates Alarms and Trips!!!!!
        # HTF Inlet Heating Rate over 1 minute time
        hr_htf = self._compute_heating_rate(self._htf_in_hist, t_htf, ts_sec)
        self.hr_htf_in.v = hr_htf
        alarm, trip = self._alarm_trip(
            abs(hr_htf),
            self._safe(self.high_hr_htf_in_alarm.v),
            self._safe(self.high_hr_htf_in_trip.v),
            high=True,
        )
        self.alarm_high_hr_htf_in.v = alarm
        self.trip_high_hr_htf_in.v = trip

        # Outlet FW Heating Rate over 1 minute time
        hr_fw = self._compute_heating_rate(self._fw_out_hist, t_fw_out, ts_sec)
        self.hr_fw_out.v = hr_fw
        alarm, trip = self._alarm_trip(
            abs(hr_fw),
            self._safe(self.high_hr_fw_out_alarm.v),
            self._safe(self.high_hr_fw_out_trip.v),
            high=True,
        )
        self.alarm_high_hr_fw_out.v = alarm
        self.trip_high_hr_fw_out.v = trip

        self._advance_history(self._htf_in_hist, t_htf)
        self._advance_history(self._fw_out_hist, t_fw_out)

    def calculate(self):
        # Read parameters and inputs (Fortran ordering).
        heat_transfer_rated = max(self._safe(self.heat_transfer_rated.v), 0.0)
        m_dot_htf_rated = max(self._safe(self.m_dot_htf_rated.v, 1.0), 1.0e-9)
        rated_exp = self._safe(self.rated_exp.v, 0.8)
        no_shell_passes = max(self._safe(self.no_shell_passes.v, 1.0), 1.0)
        no_tube_passes = max(self._safe(self.no_tube_passes.v, 1.0), 1.0)
        length_hx = max(self._safe(self.length_hx.v, 1.0), 1.0e-9)
        tube_od = max(self._safe(self.tube_od.v, 0.02), 1.0e-9)
        tube_th = max(self._safe(self.tube_th.v, 0.001), 0.0)
        no_tubes = max(self._safe(self.no_tubes.v, 1.0), 1.0)
        fluid_id = self._safe(self.fluid_id.v, "Nitrate Salt")
        fluid_name = str(fluid_id)

        m_dot_fw = max(self._safe(self.m_dot_fw.v), 0.0)
        p_fw = self._safe(self.p_fw.v)
        h_fw = self._safe(self.h_fw.v)
        m_dot_htf = max(self._safe(self.m_dot_htf.v), 0.0)
        p_htf = self._safe(self.p_htf.v)
        t_htf = self._safe(self.t_htf.v, 500.0)
        if t_htf == 0.0:
            # default value until actual temperature enters
            t_htf = 500.0

        is_start_time = (
            abs(self._safe(getattr(self.model, "time", 0.0)) - self._safe(getattr(self.model.settings, "start_time", 0.0))) < 1.0e-12
        )
        if is_start_time and not self.model.is_converged:
            # Do all of the first timestep manipulations here.
            # Set the initial values of outputs 1..11.
            self._set_outputs_1_to_11(m_dot_fw, 0.0, p_fw, h_fw, 0.0, m_dot_htf, 0.0, p_htf, t_htf, 0.0, 0.0)
            self._set_alarm_outputs_zero()
            return

        if self.model.is_converged:
            # Do any end-of-timestep manipulations here.
            self._run_end_of_timestep_checks(m_dot_htf, p_htf, t_htf, p_fw, h_fw)
            return

        if m_dot_htf > 0.01:
            if m_dot_fw > 0.01:
                if p_fw > 0.0:
                    # Hydraulic calculations are reasonable.
                    t_sat = float(fp.temperature("water", P=max(p_fw, 1.0), Q=0.0))
                    h_sat_f = float(fp.enthalpy("water", P=max(p_fw, 1.0), Q=0.0))
                    h_sat_g = float(fp.enthalpy("water", P=max(p_fw, 1.0), Q=1.0))
                    t_fw = float(fp.temperature("water", P=max(p_fw, 1.0), h=max(h_fw, 1.0)))
                    cp_htf = float(self._props.specheat(fluid_name, t_htf, p_htf))

                    # solve for specific heat values of fw and htf
                    if abs(t_fw - t_sat) > 1.0:
                        cp_fw = float(fp.specheat("water", T=max(t_fw, 273.15), P=max(p_fw, 1.0)))
                    else:
                        # T_fw is too close to saturation temperature to give accurate cp value
                        cp_fw = float(fp.specheat("water", T=max(0.5 * (t_htf + t_fw), 273.15), P=max(p_fw, 1.0)))

                    # Find effectiveness of heat exchanger based on inlet conditions of HTF and FW.
                    # Surface area of the feedwater side of the heat exchanger.
                    a_s = math.pi * max(tube_od - 2.0 * tube_th, 1.0e-9) * length_hx * no_tube_passes * no_tubes
                    ua_rated = heat_transfer_rated * a_s
                    ua_od = ua_rated * (m_dot_htf / m_dot_htf_rated) ** rated_exp

                    cap_fw = m_dot_fw * cp_fw
                    cap_htf = m_dot_htf * cp_htf
                    cap_min = max(min(cap_fw, cap_htf), 1.0e-9)
                    cap_max = max(cap_fw, cap_htf, 1.0e-9)
                    cr = max(cap_min / cap_max, 0.001)
                    ntu_od = ua_od / cap_min

                    sqrt_term = math.sqrt(1.0 + cr**2)
                    exp_term = math.exp(-ntu_od * sqrt_term)
                    ratio_term = (1.0 + exp_term) / max(1.0 - exp_term, 1.0e-9)
                    eta_1pass = 2.0 / max(1.0 + cr + sqrt_term * ratio_term, 1.0e-9)

                    shell_ratio = (1.0 - eta_1pass * cr) / max(1.0 - eta_1pass, 1.0e-9)
                    shell_pow = shell_ratio**no_shell_passes
                    eta_od = (shell_pow - 1.0) / max(shell_pow - cr, 1.0e-9)
                    eta_od = max(min(eta_od, 1.0), 0.0)

                    h_fw_out_s = float(fp.enthalpy("water", P=max(p_fw, 1.0), T=max(t_htf, 273.15)))

                    # Check that HTF temp is higher than FW temp.
                    if t_fw < t_htf:
                        # Heat Transfer is going the correct way.
                        q_dot_hx = min(
                            eta_od * m_dot_fw * max(h_fw_out_s - h_fw, 0.0),
                            eta_od * m_dot_htf * cp_htf * max(t_htf - t_fw, 0.0),
                        )
                        h_fw_out = max((m_dot_fw * h_fw + q_dot_hx) / max(m_dot_fw, 1.0e-9), h_fw)
                        # Flow is entering subcooled: make sure it is not passing T_sat.
                        if h_fw < h_sat_f and h_fw_out > h_sat_f:
                            h_fw_out = h_sat_f
                            q_dot_hx = m_dot_fw * (h_fw_out - h_fw)
                        t_fw_out = float(fp.temperature("water", P=max(p_fw, 1.0), h=max(h_fw_out, 1.0)))
                        # energy balance on HTF side of heat exchanger
                        t_htf_out = (m_dot_htf * cp_htf * t_htf - q_dot_hx) / max(m_dot_htf * cp_htf, 1.0e-9)
                    else:
                        # htf temperature is lower than feedwater temperature,
                        # heat transfer is going the wrong way
                        q_dot_hx = min(
                            eta_od * m_dot_fw * max(h_fw - h_fw_out_s, 0.0),
                            eta_od * m_dot_htf * cp_htf * max(t_fw - t_htf, 0.0),
                        )
                        if m_dot_fw < 1.0:
                            q_dot_hx = 0.0
                            # enthalpy out = enthalpy in
                            h_fw_out = h_fw
                            t_fw_out = t_fw
                        else:
                            h_fw_out = (m_dot_fw * h_fw - q_dot_hx) / max(m_dot_fw, 1.0e-9)
                            # if this type is the superheater, do not allow it to go below saturation
                            if h_fw >= h_sat_g:
                                h_fw_out = max(h_fw_out, h_sat_g)
                            t_fw_out = float(fp.temperature("water", P=max(p_fw, 1.0), h=max(h_fw_out, 1.0)))
                        # energy balance on HTF side of heat exchanger
                        t_htf_out = (m_dot_htf * cp_htf * t_htf + q_dot_hx) / max(m_dot_htf * cp_htf, 1.0e-9)

                    # calculating the volumetric flow rates
                    rho_fw = float(fp.density("water", P=max(p_fw, 1.0), h=max(h_fw_out, 1.0)))
                    vol_dot_fw = m_dot_fw / rho_fw if rho_fw > 0.0 else 0.0

                    rho_htf = float(self._props.density(fluid_name, t_htf_out, p_htf))
                    vol_dot_htf = m_dot_htf / max(rho_htf, 1.0e-9)

                    self._set_outputs_1_to_11(
                        m_dot_fw,
                        vol_dot_fw,
                        p_fw,
                        h_fw_out,
                        t_fw_out,
                        m_dot_htf,
                        vol_dot_htf,
                        p_htf,
                        t_htf_out,
                        q_dot_hx,
                        eta_od,
                    )
                else:
                    # Pressure is not possible, need to wait for next iteration to compute temperatures.
                    # Keep values the same.
                    t_fw = float(fp.temperature("water", P=max(p_fw, 1.0), h=max(h_fw, 1.0)))
                    self._set_outputs_1_to_11(
                        m_dot_fw,
                        0.0,
                        p_fw,
                        h_fw,
                        t_fw,
                        m_dot_htf,
                        0.0,
                        p_htf,
                        t_htf,
                        0.0,
                        0.0,
                    )
            else:
                # no FW Flow entering the system: set HTF outputs as-through
                t_htf_out = t_htf
                rho_htf = float(self._props.density(fluid_name, t_htf, p_htf))
                vol_dot_htf = m_dot_htf / max(rho_htf, 1.0e-9)

                t_fw_out = float(fp.temperature("water", P=max(p_fw, 1.0), h=max(h_fw, 1.0)))
                vol_dot_fw = m_dot_fw
                h_fw_out = h_fw

                self._set_outputs_1_to_11(
                    m_dot_fw,
                    vol_dot_fw,
                    p_fw,
                    h_fw_out,
                    t_fw_out,
                    m_dot_htf,
                    vol_dot_htf,
                    p_htf,
                    t_htf_out,
                    0.0,
                    0.0,
                )
        else:
            # no HTF Flow entering the system, set feedwater outlet the same as inlet
            t_htf_out = t_htf
            vol_dot_htf = 0.0

            rho_fw = float(fp.density("water", P=max(p_fw, 1.0), h=max(h_fw, 1.0)))
            vol_dot_fw = m_dot_fw / rho_fw if rho_fw > 0.0 else 0.0
            t_fw_out = float(fp.temperature("water", P=max(p_fw, 1.0), h=max(h_fw, 1.0)))

            self._set_outputs_1_to_11(
                m_dot_fw,
                vol_dot_fw,
                p_fw,
                h_fw,
                t_fw_out,
                m_dot_htf,
                vol_dot_htf,
                p_htf,
                t_htf_out,
                0.0,
                0.0,
            )
