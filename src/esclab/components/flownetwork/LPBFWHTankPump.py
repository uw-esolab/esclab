"""Type 6014 LPBFWH tank-pump converted from Fortran."""

import math

from eeslib import fluid_properties as fp

from esclab.simulate import Component
from .TurbinesBypassNetwork import pb_cv_data


class LPBFWHTankPump(Component):
    """
    TRNSYS Type 6014: LPBFWH Tank-Pump.

    Parameters
    ----------
    h_tank : float
        Height of the vertical cylindrical receiver tank [m].
    d_tank : float
        Diameter of the vertical cylindrical receiver tank [m].
    perc_tank_ini : float
        Initial tank fill fraction [0-1].
    t_tank_ini : float
        Initial tank temperature [K].
    valve_speed : float
        Maximum valve slew rate [deg/s].
    d_valve : float
        Valve diameter [m].
    valve_type : float
        Valve type identifier (1=concentric butterfly, 2=triple offset butterfly).
    pc_coef_a, pc_coef_b, pc_coef_c : float
        Pump head polynomial coefficients where
        ``pump_head = a*flow^2 + b*flow + c`` and flow is volumetric [m^3/s].
    eta_coef_a, eta_coef_b, eta_coef_c, eta_coef_d : float
        Pump efficiency polynomial coefficients where
        ``eta = a*flow^4 + b*flow^3 + c*flow^2 + d*flow``.

    Inputs
    ------
    auto_control : float
        Manual/automatic mode flag (1=automatic).
    pump_power : float
        Pump power state (1=on, otherwise off).
    vp_input : float
        Valve position request [0-1].
    turbine_on : float
        Turbine on signal (retained for direct Fortran mapping).
    m_dot_fw_in : float
        Feedwater mass flow to be mixed with pump discharge [kg/s].
    h_fw_in : float
        Feedwater specific enthalpy [J/kg].
    p_fw_in : float
        Feedwater pressure [Pa].
    m_dot_bfwh : float
        BFWH inflow mass flow to receiver tank [kg/s].
    h_bfwh : float
        BFWH inflow specific enthalpy [J/kg].
    p_bfwh : float
        BFWH inflow pressure [Pa].
    pid_signal : float
        PID valve position signal [0-1] used in automatic mode.

    Outputs
    -------
    vp_output : float
        Valve position output [0-1].
    m_dot_fw_out : float
        Feedwater outlet mass flow [kg/s].
    h_fw_out : float
        Feedwater outlet specific enthalpy [J/kg].
    p_fw_out : float
        Feedwater outlet pressure [Pa].
    m_dot_pump : float
        Pump discharge mass flow from tank [kg/s].
    p_pump_out : float
        Pump outlet pressure [Pa].
    w_dot_pump : float
        Pump power consumption [W].
    l_tank_start : float
        Tank liquid level at start of timestep [m].
    l_tank_end : float
        Tank liquid level at end of timestep [m].
    t_tank_start : float
        Tank temperature at start of timestep [K].
    t_tank_end : float
        Tank temperature at end of timestep [K].
    m_tank_start : float
        Tank mass at start of timestep [kg].
    m_tank_end : float
        Tank mass at end of timestep [kg].
    h_tank_start : float
        Tank specific enthalpy at start of timestep [J/kg].
    h_tank_end : float
        Tank specific enthalpy at end of timestep [J/kg].
    """

    # Model Parameters (Fortran indices 1..14)
    h_tank = Component.Parameter()
    d_tank = Component.Parameter()
    perc_tank_ini = Component.Parameter()
    t_tank_ini = Component.Parameter()
    valve_speed = Component.Parameter()
    d_valve = Component.Parameter()
    valve_type = Component.Parameter()
    pc_coef_a = Component.Parameter()
    pc_coef_b = Component.Parameter()
    pc_coef_c = Component.Parameter()
    eta_coef_a = Component.Parameter()
    eta_coef_b = Component.Parameter()
    eta_coef_c = Component.Parameter()
    eta_coef_d = Component.Parameter()

    # Model Inputs (Fortran indices 1..11)
    auto_control = Component.Input()
    pump_power = Component.Input()
    vp_input = Component.Input()
    turbine_on = Component.Input()
    m_dot_fw_in = Component.Input()
    h_fw_in = Component.Input()
    p_fw_in = Component.Input()
    m_dot_bfwh = Component.Input()
    h_bfwh = Component.Input()
    p_bfwh = Component.Input()
    pid_signal = Component.Input()

    # Model Outputs (Fortran indices 1..15)
    vp_output = Component.Output()
    m_dot_fw_out = Component.Output()
    h_fw_out = Component.Output()
    p_fw_out = Component.Output()
    m_dot_pump = Component.Output()
    p_pump_out = Component.Output()
    w_dot_pump = Component.Output()
    l_tank_start = Component.Output()
    l_tank_end = Component.Output()
    t_tank_start = Component.Output()
    t_tank_end = Component.Output()
    m_tank_start = Component.Output()
    m_tank_end = Component.Output()
    h_tank_start = Component.Output()
    h_tank_end = Component.Output()

    @staticmethod
    def _safe(value, default=0.0):
        try:
            value_f = float(value)
        except Exception:
            return float(default)
        if not math.isfinite(value_f):
            return float(default)
        return value_f

    def calculate(self):
        # Do Any "After Convergence" Manipulations That May Be Required at the End of Each Timestep
        if self.model.is_converged:
            # !!!!UPDATE TANK CONDITIONS!!!!
            self.l_tank_start.v = self._safe(self.l_tank_end.v, 0.0)
            self.t_tank_start.v = self._safe(self.t_tank_end.v, self._safe(self.t_tank_ini.v, 300.0))
            self.m_tank_start.v = self._safe(self.m_tank_end.v, 0.0)
            self.h_tank_start.v = self._safe(self.h_tank_end.v, self._safe(self.h_bfwh.v, 0.0))
            return

        # Do All of the First Timestep Manipulations Here - There Are No Iterations at the Initial Time
        if self.model.is_first_step:
            h_tank = self._safe(self.h_tank.v, 1.0)
            d_tank = max(self._safe(self.d_tank.v, 1.0), 1.0e-6)
            perc_tank_ini = self._safe(self.perc_tank_ini.v, 0.0)
            t_tank_ini = min(max(self._safe(self.t_tank_ini.v, 300.0), 273.15), 373.15)
            vp_input = self._safe(self.vp_input.v, 0.0)

            # determine initial tank enthalpy
            p_tank = 101325.0  # Pressure of the tank is set at atmospheric pressure
            h_tank_start = float(fp.enthalpy("water", P=max(p_tank, 1.0), T=max(t_tank_ini, 273.15)))

            # determine the amount of mass in the tank at the start of the timestep
            rho_water = 1000.0
            cross_section = math.pi * (d_tank / 2.0) ** 2
            l_tank_start = max(perc_tank_ini, 0.0) * h_tank
            l_tank_start = min(max(l_tank_start, 0.0), h_tank)
            m_tank_start = l_tank_start * cross_section * rho_water

            self.vp_output.v = vp_input
            self.l_tank_start.v = l_tank_start
            self.l_tank_end.v = l_tank_start
            self.t_tank_start.v = t_tank_ini
            self.t_tank_end.v = t_tank_ini
            self.m_tank_start.v = m_tank_start
            self.m_tank_end.v = m_tank_start
            self.h_tank_start.v = h_tank_start
            self.h_tank_end.v = h_tank_start
            return

        # Read the Inputs
        auto_control = self._safe(self.auto_control.v, 0.0)
        pump_power = self._safe(self.pump_power.v, 0.0)
        vp_input = self._safe(self.vp_input.v, 0.0)
        m_dot_fw_in = self._safe(self.m_dot_fw_in.v, 0.0)
        h_fw_in = self._safe(self.h_fw_in.v, 0.0)
        p_fw_in = self._safe(self.p_fw_in.v, 101325.0)
        m_dot_bfwh = self._safe(self.m_dot_bfwh.v, 0.0)
        h_bfwh = self._safe(self.h_bfwh.v, 0.0)
        pid_signal = self._safe(self.pid_signal.v, vp_input)

        # Read parameters used in normal iterations
        h_tank = max(self._safe(self.h_tank.v, 1.0), 1.0e-6)
        d_tank = max(self._safe(self.d_tank.v, 1.0), 1.0e-6)
        valve_speed = max(self._safe(self.valve_speed.v, 0.0), 0.0)
        d_valve = max(self._safe(self.d_valve.v, 1.0e-3), 1.0e-6)
        valve_type = int(round(self._safe(self.valve_type.v, 1.0)))
        pc_coef_a = self._safe(self.pc_coef_a.v, 0.0)
        pc_coef_b = self._safe(self.pc_coef_b.v, 0.0)
        pc_coef_c = self._safe(self.pc_coef_c.v, 0.0)
        eta_coef_a = self._safe(self.eta_coef_a.v, 0.0)
        eta_coef_b = self._safe(self.eta_coef_b.v, 0.0)
        eta_coef_c = self._safe(self.eta_coef_c.v, 0.0)
        eta_coef_d = self._safe(self.eta_coef_d.v, 0.7)

        p_tank = 101325.0
        rho_water = 1000.0

        # retrieve initial tank values
        vp_output = self._safe(self.vp_output.v, vp_input)
        l_tank_start = self._safe(self.l_tank_start.v, 0.0)
        t_tank_start = self._safe(self.t_tank_start.v, self._safe(self.t_tank_ini.v, 300.0))
        m_tank_start = max(self._safe(self.m_tank_start.v, 0.0), 0.0)
        h_tank_start = self._safe(self.h_tank_start.v, h_bfwh)
        ts = max(self._safe(self.model.settings.timestep, 0.0) * 3600.0, 0.0)

        # max pump head and max pump flow based on pump curve
        tol = 1000.0  # tolerance for while loops
        learning_rate = 0.4
        coef_c_adj = max(pc_coef_c - (p_fw_in - p_tank) / rho_water / 9.81, 1.0)
        disc = pc_coef_b**2 - 4.0 * pc_coef_a * coef_c_adj
        if abs(pc_coef_a) > 1.0e-12 and disc >= 0.0:
            m_dot_pump_max = ((-pc_coef_b - math.sqrt(disc)) / (2.0 * pc_coef_a)) * rho_water
            if m_dot_pump_max <= 0.0:
                m_dot_pump_max = ((-pc_coef_b + math.sqrt(disc)) / (2.0 * pc_coef_a)) * rho_water
        elif abs(pc_coef_b) > 1.0e-12:
            m_dot_pump_max = max(-coef_c_adj / pc_coef_b, 1.0e-6) * rho_water
        else:
            m_dot_pump_max = 1.0
        m_dot_pump_max = max(m_dot_pump_max, 1.0e-6)
        m_dot_pump_min = 1.0e-8

        # !!!!!Update valve position based on PID controller or manual control if First Iteration in Timestep!!!!!
        if auto_control == 1.0:
            # Automatic Mode
            pump_power = 1.0
            if self.model.is_first_iteration:
                requested_vp = pid_signal
                if self.model.is_first_step:
                    vp_output = requested_vp
                else:
                    vp_input_d = requested_vp * 90.0
                    vp_output_d = vp_output * 90.0
                    if vp_input_d > vp_output_d:
                        vp_output = min(vp_output_d + valve_speed * ts, vp_input_d) / 90.0
                    elif vp_input_d < vp_output_d:
                        vp_output = max(vp_output_d - valve_speed * ts, vp_input_d) / 90.0
                    else:
                        vp_output = requested_vp
                if vp_output <= 1.0e-5:
                    vp_output = 1.0e-5
        else:
            # Manual Control
            if self.model.is_first_iteration:
                requested_vp = vp_input
                if self.model.is_first_step:
                    vp_output = requested_vp
                else:
                    vp_input_d = requested_vp * 90.0
                    vp_output_d = vp_output * 90.0
                    if vp_input_d > vp_output_d:
                        vp_output = min(vp_output_d + valve_speed * ts, vp_input_d) / 90.0
                    elif vp_input_d < vp_output_d:
                        vp_output = max(vp_output_d - valve_speed * ts, vp_input_d) / 90.0
                    else:
                        vp_output = requested_vp
                if vp_output <= 1.0e-5:
                    vp_output = 1.0e-5

        if pump_power != 1.0:
            # !!!Pump is OFF!!!
            m_dot_pump = 0.0
            p_pump_out = p_fw_in
            h_pump_out = h_fw_in
            w_dot_pump = 0.0
        else:
            # !!!Pump is ON!!!
            # Find the amount of flow leaving the pump based on the valve position
            cv = max(pb_cv_data(valve_type, d_valve, vp_output), 1.0e-9)
            m_dot_pump = max(self._safe(self.m_dot_pump.v, m_dot_pump_min), m_dot_pump_min)
            error = tol + 1.0
            while_iterations = 0
            m_dot_pump_prev = m_dot_pump
            error_prev = error
            while abs(error) > tol and while_iterations < 200:
                while_iterations += 1
                vol_in = m_dot_pump / rho_water  # Volumetric flow entering the valve [m^3/s]
                vol_in_gpm = max(vol_in * 15850.3, 1.0e-9)  # Volumetric flow entering valve [GPM]
                delta_p_psi = (vol_in_gpm / cv) ** 2  # pressure drop in psi across valve
                delta_p = delta_p_psi * 6894.76  # pressure drop in Pa across valve
                p_pump_guess = p_fw_in + delta_p - p_tank - l_tank_start * 1000.0 * 9.81
                p_pump_out = (pc_coef_a * vol_in**2 + pc_coef_b * vol_in + pc_coef_c) * 9.81 * rho_water + p_tank
                error = p_pump_guess - p_pump_out

                if m_dot_pump < 0.001 and error > 0.0:
                    m_dot_pump = m_dot_pump_min
                    error = tol / 2.0
                if abs(m_dot_pump - m_dot_pump_max) < 0.001 and error < 0.0:
                    m_dot_pump = m_dot_pump_max
                    error = tol / 2.0

                if while_iterations == 1:
                    if abs(error) > tol:
                        m_dot_pump_prev = m_dot_pump
                        error_prev = error
                        if error < 0.0:
                            m_dot_pump = min(m_dot_pump + 1.0, m_dot_pump_max)
                        else:
                            m_dot_pump = max(m_dot_pump - 1.0, m_dot_pump_min)
                else:
                    if abs(error) > tol:
                        if m_dot_pump_prev != m_dot_pump:
                            slope = (error_prev - error) / (m_dot_pump_prev - m_dot_pump)
                            y_int = error_prev - slope * m_dot_pump_prev
                            if slope != 0.0:
                                m_dot_new = -y_int / slope
                                m_dot_new = min(max(m_dot_new, m_dot_pump_min), m_dot_pump_max)
                                m_dot_pump_prev = m_dot_pump
                                error_prev = error
                                m_dot_pump = m_dot_pump + (m_dot_new - m_dot_pump) * learning_rate
                            else:
                                m_dot_pump_prev = m_dot_pump
                                error_prev = error
                                if error < 0.0:
                                    m_dot_pump = min(m_dot_pump + 1.0, m_dot_pump_max)
                                else:
                                    m_dot_pump = max(m_dot_pump - 1.0, m_dot_pump_min)
                        else:
                            m_dot_pump_prev = m_dot_pump
                            error_prev = error
                            if error < 0.0:
                                m_dot_pump = min(m_dot_pump + 1.0, m_dot_pump_max)
                            else:
                                m_dot_pump = max(m_dot_pump - 1.0, m_dot_pump_min)

            if while_iterations >= 200 and abs(error) > tol:
                # TODO: revisit root-finding robustness if this loop frequently reaches iteration cap.
                m_dot_pump = min(max(m_dot_pump, m_dot_pump_min), m_dot_pump_max)

        # !!!!COMPLETE PUMP CALCS!!!!
        if m_dot_pump != m_dot_pump_min and m_dot_pump > 0.0:
            h_pump_in = h_tank_start
            p_pump_in = p_tank
            flow = m_dot_pump / rho_water
            eta_pump = max(eta_coef_a * flow**4 + eta_coef_b * flow**3 + eta_coef_c * flow**2 + eta_coef_d * flow, 0.2)
            s_pump_in = float(fp.entropy("water", P=max(p_pump_in, 1.0), h=max(h_pump_in, 1.0)))
            h_pump_out_s = float(fp.enthalpy("water", P=max(p_pump_out, 1.0), s=s_pump_in))
            w_dot_pump_s = m_dot_pump * (h_pump_out_s - h_pump_in)
            w_dot_pump = w_dot_pump_s / max(eta_pump, 1.0e-9)
            h_pump_out = h_pump_in + w_dot_pump / max(m_dot_pump, 1.0e-9)
        else:
            p_pump_out = p_fw_in
            h_pump_out = h_fw_in
            w_dot_pump = 0.0

        # !!!!Complete mixing tank calcs!!!!
        cross_section = math.pi * (d_tank / 2.0) ** 2
        m_tank_end = m_tank_start + m_dot_bfwh * ts - m_dot_pump * ts
        m_tank_end = max(m_tank_end, 1.0e-9)
        l_tank_end = m_tank_end / rho_water / cross_section
        l_tank_end = min(max(l_tank_end, 0.0), h_tank)
        h_tank_end = (m_tank_start * h_tank_start + m_dot_bfwh * ts * h_bfwh - m_dot_pump * ts * h_tank_start) / m_tank_end
        t_tank_end = float(fp.temperature("water", P=max(p_tank, 1.0), h=max(h_tank_end, 1.0)))

        # !!!!COMPLETE FEEDWATER MIXING CALCS!!!!
        m_dot_fw_out = m_dot_fw_in + m_dot_pump
        p_fw_out = p_fw_in
        if m_dot_fw_out > 0.0:
            h_fw_out = (m_dot_fw_in * h_fw_in + m_dot_pump * h_pump_out) / m_dot_fw_out
        else:
            h_fw_out = h_fw_in

        self.vp_output.v = vp_output
        self.m_dot_fw_out.v = m_dot_fw_out
        self.h_fw_out.v = h_fw_out
        self.p_fw_out.v = p_fw_out
        self.m_dot_pump.v = m_dot_pump
        self.p_pump_out.v = p_pump_out
        self.w_dot_pump.v = w_dot_pump
        self.l_tank_end.v = l_tank_end
        self.t_tank_end.v = t_tank_end
        self.m_tank_end.v = m_tank_end
        self.h_tank_end.v = h_tank_end
