"""Type 4035 pipe converted from Fortran."""

import math

import numpy as np

from esclab.components.esol_properties import Incompressible as Inc
from esclab.components.flownetwork.SimplePipe import FricFactor_IC
from esclab.simulate import Component


class Pipe(Component):
    """
    TRNSYS Type 4035: ESOL4035-Pipe.

    Uses a nodal temperature state vector and applies the original end-of-timestep
    RK4 integration structure from the Fortran type.

    Parameters
    ----------
    diameter : float
        Pipe inside diameter [m].
    l_tot : float
        Total pipe length [m].
    n_nodes : int
        Number of thermal nodes (minimum 2, maximum 500).
    fluid_id : str
        Fluid name that maps directly to
        :class:`esclab.components.esol_properties.Incompressible`.
    roughness : float
        Pipe roughness [m].
    n_large_elbows : float
        Number of long elbows.
    n_medium_elbows : float
        Number of medium elbows.
    n_standard_elbows : float
        Number of standard elbows.
    n_contractions : float
        Number of contractions.
    n_expansions : float
        Number of expansions.
    n_gate_valves : float
        Number of gate valves.
    init_temp : float
        Initial nodal temperature [K].
    mc_mult : float
        Thermal capacitance multiplier [-].
    heat_loss : float
        Heat loss per unit length [W/m].

    Inputs
    ------
    temperature : float
        Inlet temperature [K].
    mass_flow : float
        Mass flow [kg/s].
    t_amb : float
        Ambient temperature [K] (reserved in Type4035 body).
    pressure : float
        Inlet pressure [Pa].
    wind : float
        Wind speed [m/s] (reserved in Type4035 body).
    mass_counter : float
        Upstream mass counter [kg] (reserved in Type4035 body).

    Outputs
    -------
    temperature_out : float
        Outlet temperature [K].
    pressure_out : float
        Outlet pressure [Pa].
    mass_flow_out : float
        Outlet mass flow [kg/s].
    mass_counter_out : float
        Pipe hold-up mass [kg].
    """

    diameter = Component.Parameter()
    l_tot = Component.Parameter()
    n_nodes = Component.Parameter()
    fluid_id = Component.Parameter()
    roughness = Component.Parameter()
    n_large_elbows = Component.Parameter()
    n_medium_elbows = Component.Parameter()
    n_standard_elbows = Component.Parameter()
    n_contractions = Component.Parameter()
    n_expansions = Component.Parameter()
    n_gate_valves = Component.Parameter()
    init_temp = Component.Parameter()
    mc_mult = Component.Parameter()
    heat_loss = Component.Parameter()

    temperature = Component.Input()
    mass_flow = Component.Input()
    t_amb = Component.Input()
    pressure = Component.Input()
    wind = Component.Input()
    mass_counter = Component.Input()

    temperature_out = Component.Output()
    pressure_out = Component.Output()
    mass_flow_out = Component.Output()
    mass_counter_out = Component.Output()

    _ff_guess = 0.1
    _is_initialized = False
    _t_nodes = np.array([], dtype=float)
    _props = Inc()

    @staticmethod
    def _safe(value, default):
        return value if value == value else default

    def _density_sf(self, fluid_id, temperature_k, pressure_pa):
        name = str(fluid_id) if fluid_id == fluid_id else "Nitrate Salt"
        p = max(self._safe(pressure_pa, 0.0), 0.0)
        t = max(self._safe(temperature_k, 273.15), 1.0)
        try:
            return max(float(self._props.density(name, t, p)), 1.0)
        except Exception:
            return max(2090.0 - 0.636 * (t - 273.15), 1000.0)

    def _viscosity_sf(self, fluid_id, temperature_k, pressure_pa):
        name = str(fluid_id) if fluid_id == fluid_id else "Nitrate Salt"
        p = max(self._safe(pressure_pa, 0.0), 0.0)
        t = max(self._safe(temperature_k, 273.15), 1.0)
        try:
            return max(float(self._props.viscosity(name, t, p)), 1.0e-9)
        except Exception:
            return 2.5e-3

    def _cp_sf(self, fluid_id, temperature_k, pressure_pa):
        name = str(fluid_id) if fluid_id == fluid_id else "Nitrate Salt"
        p = max(self._safe(pressure_pa, 0.0), 0.0)
        t = max(self._safe(temperature_k, 273.15), 1.0)
        try:
            return max(float(self._props.specheat(name, t, p)) * 1000.0, 100.0)
        except Exception:
            return 1500.0

    def _pressure_drop(self, mass_flow, t_ref):
        # 1. PressureDrop  (Pa)
        m_dot = abs(self._safe(mass_flow, 0.0))
        if m_dot <= 0.0:
            return 0.0

        fluid = self._safe(self.fluid_id.v, "Nitrate Salt")
        d = max(self._safe(self.diameter.v, 0.0), 1.0e-6)
        rough = max(self._safe(self.roughness.v, 0.0), 1.0e-10)
        l_pipe = max(self._safe(self.l_tot.v, 0.0), 0.0)

        rho = self._density_sf(fluid, t_ref, 1.0)
        mu = self._viscosity_sf(fluid, t_ref, 1.0)
        nu = mu / max(rho, 1.0e-9)
        v_dot = m_dot / max(rho, 1.0e-9)
        u_fluid = v_dot / (math.pi * (d / 2.0) * (d / 2.0))

        re = abs(u_fluid * d / max(nu, 1.0e-12))
        if re < 2300.0:
            f = 64.0 / max(re, 1.0)
        else:
            f_val = FricFactor_IC(rough / d, max(re, 1.0), self._ff_guess)
            try:
                f = float(f_val)
            except Exception:
                f = self._ff_guess
            if not math.isfinite(f):
                f = self._ff_guess
        f = max(f, 1.0e-5)
        self._ff_guess = f

        g = 9.80665
        hl_pm = f * u_fluid * u_fluid / (2.0 * d * g)
        dp_pipe = hl_pm * rho * g * l_pipe

        n_exp = max(self._safe(self.n_expansions.v, 0.0), 0.0)
        n_con = max(self._safe(self.n_contractions.v, 0.0), 0.0)
        n_els = max(self._safe(self.n_standard_elbows.v, 0.0), 0.0)
        n_elm = max(self._safe(self.n_medium_elbows.v, 0.0), 0.0)
        n_ell = max(self._safe(self.n_large_elbows.v, 0.0), 0.0)
        n_gav = max(self._safe(self.n_gate_valves.v, 0.0), 0.0)

        d_over_f_hl = d / max(f, 1.0e-12) * hl_pm * rho * g
        dp_exp = 0.25 * rho * u_fluid * u_fluid * n_exp
        dp_con = 0.25 * rho * u_fluid * u_fluid * n_con
        dp_els = 0.9 * d_over_f_hl * n_els
        dp_elm = 0.75 * d_over_f_hl * n_elm
        dp_ell = 0.6 * d_over_f_hl * n_ell
        dp_gav = 0.19 * d_over_f_hl * n_gav

        return max(dp_pipe + dp_exp + dp_con + dp_els + dp_elm + dp_ell + dp_gav, 0.0)

    def _pipe_dtdt(self, t_nodes, vol, mass_flow, mc_mult, fluid_id, heat_loss, l_cv):
        # Loop Through Control Volumes to compute CV temperature rate
        n_nodes = int(t_nodes.size)
        dt = np.zeros(n_nodes, dtype=float)
        if n_nodes < 2:
            return dt

        dtdt_bar = np.zeros(n_nodes - 1, dtype=float)
        for n in range(n_nodes - 1):
            # Heat loss
            q_out = heat_loss * l_cv
            # Compute CV average temp and properties
            t_ave = 0.5 * (t_nodes[n] + t_nodes[n + 1])
            rho = self._density_sf(fluid_id, t_ave, 0.0)
            c = self._cp_sf(fluid_id, t_ave, 0.0)
            # Compute CV temperature rate
            dtdt_bar[n] = (mass_flow * c * (t_nodes[n] - t_nodes[n + 1]) - q_out) / max(vol * rho * c * mc_mult, 1.0e-12)

        # Compute Nodal Temperature Rates
        dt[0] = 0.0
        dt[-1] = dtdt_bar[-1]
        if n_nodes > 2:
            dt[1:-1] = 0.5 * (dtdt_bar[:-1] + dtdt_bar[1:])
        return dt

    def _initialize_state(self):
        n_nodes = int(round(self._safe(self.n_nodes.v, 2.0)))
        n_nodes = max(2, min(n_nodes, 500))
        init_temp = self._safe(self.init_temp.v, self._safe(self.temperature.v, 300.0))

        # Initialize temperatures in nodes of pipe
        self._t_nodes = np.full(n_nodes, init_temp, dtype=float)

        # Compute total mass in the pipe for mass counter
        # Define length of control volumes
        l_cv = max(self._safe(self.l_tot.v, 0.0), 0.0) / max(float(n_nodes - 1), 1.0)
        # Define Volume of CV's
        vol = math.pi * (max(self._safe(self.diameter.v, 0.0), 1.0e-6) / 2.0) ** 2 * l_cv

        mass_counter = 0.0
        fluid = self._safe(self.fluid_id.v, "Nitrate Salt")
        for _ in range(n_nodes - 1):
            t_cv = init_temp
            mass_counter += vol * self._density_sf(fluid, t_cv, 0.0)

        d_p = self._pressure_drop(self._safe(self.mass_flow.v, 0.0), init_temp)

        self.temperature_out.v = init_temp
        self.pressure_out.v = self._safe(self.pressure.v, 0.0) - d_p
        self.mass_flow_out.v = self._safe(self.mass_flow.v, 0.0)
        self.mass_counter_out.v = mass_counter
        self._is_initialized = True

    def _advance_end_of_timestep(self):
        if self._t_nodes.size < 2:
            return

        n_nodes = self._t_nodes.size
        l_cv = max(self._safe(self.l_tot.v, 0.0), 0.0) / max(float(n_nodes - 1), 1.0)
        vol = math.pi * (max(self._safe(self.diameter.v, 0.0), 1.0e-6) / 2.0) ** 2 * l_cv
        mass_flow = self._safe(self.mass_flow.v, 0.0)
        mc_mult = max(self._safe(self.mc_mult.v, 1.0), 1.0e-9)
        fluid_id = self._safe(self.fluid_id.v, "Nitrate Salt")
        heat_loss = max(self._safe(self.heat_loss.v, 0.0), 0.0)

        model = getattr(self, "model", None)
        timestep_hours = self._safe(getattr(getattr(model, "settings", None), "timestep", 1.0), 1.0)
        timestep_s = max(timestep_hours * 3600.0, 1.0e-9)

        # Load in temperatures from last timestep from dynamic array
        t_hat = self._t_nodes.copy()
        # Update Input Temperature
        t_hat[0] = self._safe(self.temperature.v, t_hat[0])
        t_prev = t_hat.copy()

        # Step through time with RK-4
        k1 = self._pipe_dtdt(t_hat, vol, mass_flow, mc_mult, fluid_id, heat_loss, l_cv)
        t_hat = t_prev + k1 * timestep_s / 2.0

        k2 = self._pipe_dtdt(t_hat, vol, mass_flow, mc_mult, fluid_id, heat_loss, l_cv)
        t_hat = t_prev + k2 * timestep_s / 2.0

        k3 = self._pipe_dtdt(t_hat, vol, mass_flow, mc_mult, fluid_id, heat_loss, l_cv)
        t_hat = t_prev + k3 * timestep_s

        k4 = self._pipe_dtdt(t_hat, vol, mass_flow, mc_mult, fluid_id, heat_loss, l_cv)

        self._t_nodes = t_prev + (k1 / 6.0 + k2 / 3.0 + k3 / 3.0 + k4 / 6.0) * timestep_s

    def _mass_counter_pipe(self):
        if self._t_nodes.size < 2:
            return 0.0
        n_nodes = self._t_nodes.size
        l_cv = max(self._safe(self.l_tot.v, 0.0), 0.0) / max(float(n_nodes - 1), 1.0)
        vol = math.pi * (max(self._safe(self.diameter.v, 0.0), 1.0e-6) / 2.0) ** 2 * l_cv
        fluid = self._safe(self.fluid_id.v, "Nitrate Salt")

        # Loop through each control volume
        mass_counter = 0.0
        for n in range(n_nodes - 1):
            # Compute control volume average temperature
            t_cv = 0.5 * (self._t_nodes[n] + self._t_nodes[n + 1])
            # Compute mass in CV
            mass_counter += vol * self._density_sf(fluid, t_cv, 0.0)
        return mass_counter

    def calculate(self):
        model = getattr(self, "model", None)
        is_first_step = bool(getattr(model, "is_first_step", False))
        is_converged = bool(getattr(model, "is_converged", False))

        if is_first_step or not self._is_initialized:
            self._initialize_state()
            return

        # Perform Thermal Calculations at the End of Each Timestep
        if is_converged:
            self._advance_end_of_timestep()
            return

        # Compute total mass in the pipe for mass counter
        mass_counter = self._mass_counter_pipe()

        # Compute pressure drop
        t_ave = 0.5 * (self._t_nodes[0] + self._t_nodes[-1])
        d_p = self._pressure_drop(self._safe(self.mass_flow.v, 0.0), t_ave)

        self.temperature_out.v = float(self._t_nodes[-1])
        self.pressure_out.v = self._safe(self.pressure.v, 0.0) - d_p
        self.mass_flow_out.v = self._safe(self.mass_flow.v, 0.0)
        self.mass_counter_out.v = mass_counter