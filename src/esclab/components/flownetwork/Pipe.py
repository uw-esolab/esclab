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

    def _pressure_drop(self, mass_flow, t_ref):
        # 1. PressureDrop  (Pa)
        m_dot = mass_flow

        fluid_name = str(self.fluid_id.v)
        d = self.diameter.v
        rough = self.roughness.v
        l_pipe = self.l_tot.v

        rho = float(self._props.density(fluid_name, t_ref, 1.0))
        mu = float(self._props.viscosity(fluid_name, t_ref, 1.0))
        nu = mu / rho
        v_dot = m_dot / rho
        u_fluid = v_dot / (math.pi * (d / 2.0) * (d / 2.0))

        re = abs(u_fluid * d / nu)
        if re < 2300.0:
            f = 64.0 / re
        else:
            f = float(FricFactor_IC(rough / d, re, self._ff_guess))
        self._ff_guess = f

        g = 9.80665
        hl_pm = f * u_fluid * u_fluid / (2.0 * d * g)
        dp_pipe = hl_pm * rho * g * l_pipe

        n_exp = self.n_expansions.v
        n_con = self.n_contractions.v
        n_els = self.n_standard_elbows.v
        n_elm = self.n_medium_elbows.v
        n_ell = self.n_large_elbows.v
        n_gav = self.n_gate_valves.v

        d_over_f_hl = d / f * hl_pm * rho * g
        dp_exp = 0.25 * rho * u_fluid * u_fluid * n_exp
        dp_con = 0.25 * rho * u_fluid * u_fluid * n_con
        dp_els = 0.9 * d_over_f_hl * n_els
        dp_elm = 0.75 * d_over_f_hl * n_elm
        dp_ell = 0.6 * d_over_f_hl * n_ell
        dp_gav = 0.19 * d_over_f_hl * n_gav

        return dp_pipe + dp_exp + dp_con + dp_els + dp_elm + dp_ell + dp_gav

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
            fluid_name = str(fluid_id)
            rho = float(self._props.density(fluid_name, t_ave, 0.0))
            c = float(self._props.specheat(fluid_name, t_ave, 0.0)) * 1000.0
            # Compute CV temperature rate
            dtdt_bar[n] = (mass_flow * c * (t_nodes[n] - t_nodes[n + 1]) - q_out) / (vol * rho * c * mc_mult)

        # Compute Nodal Temperature Rates
        dt[0] = 0.0
        dt[-1] = dtdt_bar[-1]
        if n_nodes > 2:
            dt[1:-1] = 0.5 * (dtdt_bar[:-1] + dtdt_bar[1:])
        return dt

    def _initialize_state(self):
        n_nodes = int(round(self.n_nodes.v))
        init_temp = self.init_temp.v

        # Initialize temperatures in nodes of pipe
        self._t_nodes = np.full(n_nodes, init_temp, dtype=float)

        # Compute total mass in the pipe for mass counter
        # Define length of control volumes
        l_cv = self.l_tot.v / float(n_nodes - 1)
        # Define Volume of CV's
        vol = math.pi * (self.diameter.v / 2.0) ** 2 * l_cv

        mass_counter = 0.0
        fluid_name = str(self.fluid_id.v)
        for _ in range(n_nodes - 1):
            t_cv = init_temp
            mass_counter += vol * float(self._props.density(fluid_name, t_cv, 0.0))

        d_p = self._pressure_drop(self.mass_flow.v, init_temp)

        self.temperature_out.v = init_temp
        self.pressure_out.v = self.pressure.v - d_p
        self.mass_flow_out.v = self.mass_flow.v
        self.mass_counter_out.v = mass_counter
        self._is_initialized = True

    def _advance_end_of_timestep(self):
        if self._t_nodes.size < 2:
            return

        n_nodes = self._t_nodes.size
        l_cv = self.l_tot.v / float(n_nodes - 1)
        vol = math.pi * (self.diameter.v / 2.0) ** 2 * l_cv
        mass_flow = self.mass_flow.v
        mc_mult = self.mc_mult.v
        fluid_id = self.fluid_id.v
        heat_loss = self.heat_loss.v

        model = getattr(self, "model", None)
        timestep_hours = getattr(getattr(model, "settings", None), "timestep", 1.0)
        timestep_s = timestep_hours * 3600.0

        # Load in temperatures from last timestep from dynamic array
        t_hat = self._t_nodes.copy()
        # Update Input Temperature
        t_hat[0] = self.temperature.v
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
        l_cv = self.l_tot.v / float(n_nodes - 1)
        vol = math.pi * (self.diameter.v / 2.0) ** 2 * l_cv
        fluid_name = str(self.fluid_id.v)

        # Loop through each control volume
        mass_counter = 0.0
        for n in range(n_nodes - 1):
            # Compute control volume average temperature
            t_cv = 0.5 * (self._t_nodes[n] + self._t_nodes[n + 1])
            # Compute mass in CV
            mass_counter += vol * float(self._props.density(fluid_name, t_cv, 0.0))
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
        d_p = self._pressure_drop(self.mass_flow.v, t_ave)

        self.temperature_out.v = float(self._t_nodes[-1])
        self.pressure_out.v = self.pressure.v - d_p
        self.mass_flow_out.v = self.mass_flow.v
        self.mass_counter_out.v = mass_counter