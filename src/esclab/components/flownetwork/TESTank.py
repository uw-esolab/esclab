"""Type 4100 TES tank converted from Fortran."""

import math

from esclab.components.esol_properties import Incompressible as Inc

from esclab.simulate import Component


class TESTank(Component):
    """
    TRNSYS Type 4100: ESOL4100-TESTank.

    Parameters
    ----------
    d_in, height, ins_th, k_iso, emiss, t0, l0 : float
    id_fluid : str
        Fluid name passed to :class:`esclab.components.esol_properties.Incompressible`.

    Inputs
    ------
    t_in, m_in, m_out, t_env, v_env, v_air : float

    Outputs
    -------
    output_1..output_9 : float
        Pump draw, state, wall temperature, and heat loss diagnostics.
    """

    d_in = Component.Parameter()
    height = Component.Parameter()
    ins_th = Component.Parameter()
    k_iso = Component.Parameter()
    emiss = Component.Parameter()
    t0 = Component.Parameter()
    l0 = Component.Parameter()
    id_fluid = Component.Parameter()

    t_in = Component.Input()
    m_in = Component.Input()
    m_out = Component.Input()
    t_env = Component.Input()
    v_env = Component.Input()
    v_air = Component.Input()

    m_outlet = Component.Output()
    vol_dot_outlet = Component.Output()
    t_outlet = Component.Output()
    p_outlet = Component.Output()
    t_tank_out = Component.Output()
    level_out = Component.Output()
    mass_out = Component.Output()
    wall_temp_out = Component.Output()
    q_loss_out = Component.Output()

    _tank_temp = 300.0
    _tank_level = 0.0
    _wall_temp = 305.0
    _props = Inc()

    @staticmethod
    def _safe(value, default):
        return value if value == value else default

    def calculate(self):
        # Model Parameters
        d_in = max(self._safe(self.d_in.v, 1.0), 1.0e-6)
        height = max(self._safe(self.height.v, 1.0), 1.0e-6)
        ins_th = max(self._safe(self.ins_th.v, 0.0), 0.0)
        k_iso = max(self._safe(self.k_iso.v, 0.04), 1.0e-8)
        emiss = min(max(self._safe(self.emiss.v, 0.8), 0.0), 1.0)
        t0 = self._safe(self.t0.v, 600.0)
        l0 = max(self._safe(self.l0.v, 0.0), 0.0)
        fluid_id = self._safe(self.id_fluid.v, "Nitrate Salt")
        fluid_name = str(fluid_id)

        # Model Inputs
        t_in = self._safe(self.t_in.v, t0)
        m_in = max(self._safe(self.m_in.v, 0.0), 0.0)
        m_out = max(self._safe(self.m_out.v, 0.0), 0.0)
        t_env = self._safe(self.t_env.v, 300.0)
        v_env = max(self._safe(self.v_env.v, 0.0), 0.0)
        v_air = max(self._safe(self.v_air.v, 0.0), 0.0)

        p_salt = 101325.0

        # Constants in the code
        g = 9.81
        mu_air = 184.6e-7
        k_air = 26.3e-3
        rho_air = 1.1614
        visc_air = 15.89e-6
        alpha_air = 22.5e-6
        beta_air = 1.0 / max(t_env, 1.0)
        sigma = 5.67e-8
        pr_ext = 0.7
        t_crit = 10.0

        # Do All of the First Timestep Manipulations Here - There Are No Iterations at the Initial Time
        if self.model.is_first_step:
            rho_salt = max(float(self._props.density(fluid_name, t0, p_salt)), 1.0)
            area_cs = math.pi * d_in**2 / 4.0
            m_tank = rho_salt * (l0 * area_cs)
            tw = 305.0
            p_out = 101325.0 + g * rho_salt * l0

            self._tank_temp = t0
            self._tank_level = l0
            self._wall_temp = tw

            self.m_outlet.v = 0.0
            self.vol_dot_outlet.v = 0.0
            self.t_outlet.v = t0
            self.p_outlet.v = p_out
            self.t_tank_out.v = t0
            self.level_out.v = l0
            self.mass_out.v = m_tank
            self.wall_temp_out.v = tw
            self.q_loss_out.v = 0.0
            return

        # Tank level and temperature at the end of the last timestep
        t_tank_start = self._tank_temp
        l_tank_start = self._tank_level
        tw_start = self._wall_temp

        # BREAK UP TIMESTEP INTO SUB-TIMESTEPS IF NEEDED
        ts = max(self.model.settings.timestep * 3600.0, 1.0e-9)
        if ts > t_crit:
            n_sub = int(math.ceil(ts / t_crit))
            sub_ts = ts / n_sub
        else:
            n_sub = 1
            sub_ts = ts

        area_cs = math.pi * d_in**2 / 4.0
        a_total = math.pi * d_in * height + area_cs
        a_ext = math.pi * (d_in + 2.0 * ins_th) * height
        d_ext = d_in + 2.0 * ins_th

        t_tank = t_tank_start
        l_tank = l_tank_start
        tw = tw_start
        q_loss_total = 0.0

        # SALT Properties based on fluid temperature
        cp_salt_in = max(float(self._props.specheat(fluid_name, t_in, p_salt)), 100.0)

        # Here starts a "for" loop to simulate each "subtime step"
        for _ in range(n_sub):
            rho_salt = max(float(self._props.density(fluid_name, t_tank, p_salt)), 1.0)
            cp_salt = max(float(self._props.specheat(fluid_name, t_tank, p_salt)), 100.0)
            m_tank = rho_salt * (l_tank * area_cs)

            # Here starts the "while" loop to converge wall temperature
            err = 1.0
            iter_count = 0
            while err > 1.0e-6 and iter_count < 200:
                iter_count += 1

                # Thermal Resistances
                # Conductivity resistance
                r_i = ins_th / max(a_total * k_iso, 1.0e-12)

                # Radiative resistance
                h_rad = emiss * sigma * (t_env**2 + tw**2) * (t_env + tw)
                r_ii = 1.0 / max(a_ext * h_rad, 1.0e-12)

                # Convective resistance (natural/forced convection outside tank)
                if v_env < 0.1:
                    ra_ext = (beta_air * (tw - t_env) * g * height**3) / max(visc_air * alpha_air, 1.0e-12)
                    if ra_ext < 0.0:
                        nu_ext = 0.68
                    else:
                        nu_ext = 0.68 + (0.67 * ra_ext**0.25) / (1.0 + (0.492 / pr_ext) ** (9.0 / 16.0))
                else:
                    re_ext = rho_air * v_env * height / max(mu_air, 1.0e-12)
                    if re_ext < 3000.0:
                        nu_ext = 0.664 * re_ext**0.5 * pr_ext**0.333
                    else:
                        nu_ext = 0.0296 * re_ext**0.8 * pr_ext**0.333
                h_ext = k_air * nu_ext / max(height, 1.0e-12)
                r_iii = 1.0 / max(a_ext * h_ext, 1.0e-12)

                # Convective resistance (forced air around tank)
                re_air = rho_air * v_air * d_ext / max(mu_air, 1.0e-12)
                nu_air = 0.0296 * max(re_air, 0.0) ** 0.8 * pr_ext**0.333
                h_air = k_air * nu_air / max(d_ext, 1.0e-12)
                r_iv = 1.0 / max(d_ext * h_air, 1.0e-12)

                # Parallel external resistances
                r_v_num = r_ii * r_iii * r_iv
                r_v_den = r_iii * r_iv + r_ii * r_ii + r_ii * r_iv
                r_v = r_v_num / max(r_v_den, 1.0e-12)

                m_new = m_tank + (m_in - m_out) * sub_ts
                h_in = cp_salt_in * t_in
                div = (
                    m_new * cp_salt / sub_ts
                    + cp_salt * (m_in - m_out)
                    + 1.0 / max(r_i + r_v, 1.0e-12)
                    + cp_salt * m_out
                )
                t_new = (
                    (m_new * cp_salt * t_tank / sub_ts)
                    + (h_in * m_in)
                    + (t_env / max(r_i + r_v, 1.0e-12))
                ) / max(div, 1.0e-12)

                hl = (t_new - t_env) / max(r_i + r_v, 1.0e-12)

                tw_prev = tw
                tw = hl * r_v + t_env
                err = (tw - tw_prev) ** 2

            # Once wall-temperature error is small, move to next subtime step
            t_tank = t_new
            m_tank = m_new
            q_loss_total += hl * sub_ts

            rho_salt_next = max(float(self._props.density(fluid_name, t_tank, p_salt)), 1.0)
            l_tank = (m_tank / max(rho_salt_next, 1.0e-9)) / max(area_cs, 1.0e-12)

        vol_dot_out = m_out / max(float(self._props.density(fluid_name, t_tank, p_salt)), 1.0e-9)

        self.m_outlet.v = m_out
        self.vol_dot_outlet.v = vol_dot_out
        self.t_outlet.v = t_tank
        self.p_outlet.v = p_salt
        self.t_tank_out.v = t_tank
        self.level_out.v = l_tank
        self.mass_out.v = m_tank
        self.wall_temp_out.v = tw
        self.q_loss_out.v = q_loss_total

        # Store dynamic values when the timestep is converged
        if self.model.is_converged:
            self._tank_temp = t_tank
            self._tank_level = l_tank
            self._wall_temp = tw

