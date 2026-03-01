"""Type 6027 subcooled-water pumps converted from Fortran."""

import numpy as np

from esclab.simulate import Component


class SCWaterPumps(Component):
    """
    TRNSYS Type 6027: SCWaterPumps.

    Fortran uses dynamic sizing (inputs=4+2*N, outputs=9+7*N).

    This implementation keeps scalar plant-level channels and uses array-valued
    ports for the cycle-allocated pump channels.
    """

    MAX_PUMPS = 8

    reserved_parameters = Component.Parameter(np.array([]))
    pump_head_nominal = Component.Parameter()
    pump_head_per_pump = Component.Parameter(np.array([]))
    n_active_pumps = Component.Parameter()

    p_in = Component.Input()
    h_in = Component.Input()
    reserved_input_3 = Component.Input()
    pump_on = Component.Input(np.array([]))
    vp_req = Component.Input(np.array([]))
    reserved_cycle_inputs = Component.Input(np.array([]))

    m_dot_total_out = Component.Output()
    vol_dot_total_out = Component.Output()
    p_out = Component.Output()
    h_out = Component.Output()
    reserved_outputs = Component.Output()
    pump_p_out = Component.Output()
    pump_h_out = Component.Output()
    pump_m_out = Component.Output()
    pump_vp_out = Component.Output()

    @staticmethod
    def _scalar(value, default=0.0):
        return default if value != value else float(value)

    @staticmethod
    def _as_float_array(value):
        if isinstance(value, np.ndarray):
            arr = value.astype(float, copy=False).reshape(-1)
        elif value is None:
            arr = np.array([], dtype=float)
        else:
            arr = np.asarray(value, dtype=float).reshape(-1)
        return np.nan_to_num(arr, nan=0.0)

    def calculate(self):
        rho = 1000.0
        n_active_raw = self._scalar(self.n_active_pumps.v, 1.0)
        n_active = int(max(1, min(self.MAX_PUMPS, round(n_active_raw))))

        p_in = self._scalar(self.p_in.v, 0.0)
        h_in = self._scalar(self.h_in.v, 0.0)

        pump_on_arr = self._as_float_array(self.pump_on.v)
        vp_req_arr = self._as_float_array(self.vp_req.v)

        if pump_on_arr.size < n_active:
            pump_on_arr = np.pad(pump_on_arr, (0, n_active - pump_on_arr.size))
        if vp_req_arr.size < n_active:
            vp_req_arr = np.pad(vp_req_arr, (0, n_active - vp_req_arr.size))

        pump_on_arr = np.clip(pump_on_arr[:n_active], 0.0, 1.0)
        vp_out_arr = np.clip(vp_req_arr[:n_active], 0.0, 1.0)

        head_per_pump = self._as_float_array(self.pump_head_per_pump.v)
        if head_per_pump.size < n_active:
            head_default = self._scalar(self.pump_head_nominal.v, 0.0)
            head_per_pump = np.pad(head_per_pump, (0, n_active - head_per_pump.size), constant_values=head_default)
        head_per_pump = head_per_pump[:n_active]

        pump_is_on = pump_on_arr >= 0.5
        m_dot_arr = np.where(pump_is_on, rho, 0.0)
        p_out_arr = np.where(pump_is_on, head_per_pump * rho * 9.81, p_in)
        h_out_arr = np.full(n_active, h_in)

        m_dot_total = float(np.sum(m_dot_arr))
        h_out = float(np.sum(m_dot_arr * h_out_arr) / max(m_dot_total, 1.0e-9))

        self.m_dot_total_out.v = max(m_dot_total, 1.0e-10)
        self.vol_dot_total_out.v = m_dot_total / rho
        self.p_out.v = p_in
        self.h_out.v = h_out if m_dot_total > 0.0 else h_in

        self.pump_vp_out.v = vp_out_arr
        self.pump_p_out.v = p_out_arr
        self.pump_h_out.v = h_out_arr
        self.pump_m_out.v = m_dot_arr
        self.reserved_outputs.v = np.zeros(max(9 + 7 * n_active - 4 - 4 * n_active, 0), dtype=float)
