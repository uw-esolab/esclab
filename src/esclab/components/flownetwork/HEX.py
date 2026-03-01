"""Type 4101 shell-and-tube HX converted from Fortran."""

import math

import numpy as np
from esclab.components.esol_properties import Incompressible as Inc

from esclab.simulate import Component


class HEX(Component):
    """
    TRNSYS Type 4101: ESOL4101-HEX.

    Parameters
    ----------
    n_baffles, n_hex, tube_passes, shell_passes, d_shell, l_shell, r_in, r_out,
    l_baffle, th_baffle, th_ins, s_t, s_l, n_tubes, config, emiss, k_ins,
    k_steel : float
    fluid_tube, fluid_shell : str

    Inputs
    ------
    v_out, t_env,
    m_dot_shell_c, t_shell_in_c, p_shell_in_c,
    m_dot_shell_d, t_shell_in_d, p_shell_in_d,
    m_dot_tube_c, t_tube_in_c, p_tube_in_c,
    m_dot_tube_d, t_tube_in_d, p_tube_in_d,
    mass_counter_c, mass_counter_d : float

    Outputs
    -------
    mode,
    m_dot_shell_out_c, t_shell_out_c, p_shell_out_c,
    m_dot_shell_out_d, t_shell_out_d, p_shell_out_d,
    m_dot_tube_out_c, t_tube_out_c, p_tube_out_c,
    m_dot_tube_out_d, t_tube_out_d, p_tube_out_d,
    mass_counter_out_c, mass_counter_out_d : float
    """

    n_baffles = Component.Parameter()
    n_hex = Component.Parameter()
    tube_passes = Component.Parameter()
    shell_passes = Component.Parameter()
    d_shell = Component.Parameter()
    l_shell = Component.Parameter()
    r_in = Component.Parameter()
    r_out = Component.Parameter()
    l_baffle = Component.Parameter()
    th_baffle = Component.Parameter()
    th_ins = Component.Parameter()
    s_t = Component.Parameter()
    s_l = Component.Parameter()
    n_tubes = Component.Parameter()
    config = Component.Parameter()
    emiss = Component.Parameter()
    k_ins = Component.Parameter()
    k_steel = Component.Parameter()
    fluid_tube = Component.Parameter()
    fluid_shell = Component.Parameter()

    v_out = Component.Input()
    t_env = Component.Input()
    m_dot_shell_c = Component.Input()
    t_shell_in_c = Component.Input()
    p_shell_in_c = Component.Input()
    m_dot_shell_d = Component.Input()
    t_shell_in_d = Component.Input()
    p_shell_in_d = Component.Input()
    m_dot_tube_c = Component.Input()
    t_tube_in_c = Component.Input()
    p_tube_in_c = Component.Input()
    m_dot_tube_d = Component.Input()
    t_tube_in_d = Component.Input()
    p_tube_in_d = Component.Input()
    mass_counter_c = Component.Input()
    mass_counter_d = Component.Input()

    mode = Component.Output()
    m_dot_shell_out_c = Component.Output()
    t_shell_out_c = Component.Output()
    p_shell_out_c = Component.Output()
    m_dot_shell_out_d = Component.Output()
    t_shell_out_d = Component.Output()
    p_shell_out_d = Component.Output()
    m_dot_tube_out_c = Component.Output()
    t_tube_out_c = Component.Output()
    p_tube_out_c = Component.Output()
    m_dot_tube_out_d = Component.Output()
    t_tube_out_d = Component.Output()
    p_tube_out_d = Component.Output()
    mass_counter_out_c = Component.Output()
    mass_counter_out_d = Component.Output()

    _initialized = False
    _n_nodes = 1
    _mode_prev = 1.0
    _t_shell_nodes = np.array([400.0])
    _t_tube_nodes = np.array([650.0])
    _h_shell_nodes = np.array([400.0 * 1500.0])
    _h_tube_nodes = np.array([650.0 * 2200.0])
    _props = Inc()

    @staticmethod
    def _is_num(value):
        return value == value and np.isfinite(value)

    @staticmethod
    def _safe(value, default):
        return value if HEX._is_num(value) else default

    @staticmethod
    def _diag(vector):
        return np.diag(vector)

    def _node_count(self):
        n_baffles = max(int(round(self._safe(self.n_baffles.v, 1.0))), 1)
        shell_p = max(int(round(self._safe(self.shell_passes.v, 1.0))), 1)
        n_hex = max(int(round(self._safe(self.n_hex.v, 1.0))), 1)
        return n_baffles * shell_p * n_hex + 1

    @staticmethod
    def _h_from_tube_temp_c(temp_c):
        return (
            1.181138093142442e-11 * temp_c**5
            - 1.340092005730253e-08 * temp_c**4
            + 5.872766536462754e-06 * temp_c**3
            + 1.537140422596891e-04 * temp_c**2
            + 1.624749057644918 * temp_c
            - 19.52
        ) * 1.0e3

    @staticmethod
    def _vspec_from_h(h_j_kg):
        hk = h_j_kg * 1.0e-3
        return (
            5.161496408907718e-19 * hk**5
            - 4.960573299627946e-16 * hk**4
            + 3.134480901443810e-13 * hk**3
            + 7.807498873110475e-11 * hk**2
            + 4.315796357756615e-07 * hk
            + 9.381445088449658e-04
        )

    @staticmethod
    def _dv_dh(h_j_kg):
        hk = h_j_kg * 1.0e-3
        return (
            2.580748204453859e-18 * hk**4
            - 1.984229319851178e-15 * hk**3
            + 9.403442704331430e-13 * hk**2
            + 1.561499774622095e-10 * hk
            + 4.315796357756615e-07
        ) * 1.0e-3

    @staticmethod
    def _tube_temp_from_h(h_j_kg):
        hk = h_j_kg * 1.0e-3
        return (
            -1.629510314864758e-13 * hk**5
            + 3.147528019331875e-10 * hk**4
            - 1.359541381535188e-07 * hk**3
            - 1.822021668868121e-04 * hk**2
            + 0.619357049700913 * hk
            + 11.960141056529313
        ) + 273.15

    def _ensure_state(self, n_nodes):
        if n_nodes != self._n_nodes or len(self._t_shell_nodes) != n_nodes:
            self._n_nodes = n_nodes
            self._t_shell_nodes = np.full(n_nodes, 400.0, dtype=float)
            self._t_tube_nodes = np.full(n_nodes, 650.0, dtype=float)
            self._h_shell_nodes = np.full(n_nodes, 400.0 * 1500.0, dtype=float)
            self._h_tube_nodes = np.full(n_nodes, 650.0 * 2200.0, dtype=float)

    def _build_d1dz(self, n_nodes, dz):
        if n_nodes < 5:
            mat = np.zeros((n_nodes, n_nodes), dtype=float)
            for idx in range(n_nodes):
                if idx == 0 and n_nodes > 1:
                    mat[idx, idx] = -1.0 / dz
                    mat[idx, idx + 1] = 1.0 / dz
                elif idx == n_nodes - 1 and n_nodes > 1:
                    mat[idx, idx - 1] = -1.0 / dz
                    mat[idx, idx] = 1.0 / dz
                elif 0 < idx < n_nodes - 1:
                    mat[idx, idx - 1] = -0.5 / dz
                    mat[idx, idx + 1] = 0.5 / dz
            return mat

        mat = np.zeros((n_nodes, n_nodes), dtype=float)
        for idx in range(2, n_nodes - 2):
            mat[idx, idx - 2] = 1.0
            mat[idx, idx - 1] = -8.0
            mat[idx, idx + 1] = 8.0
            mat[idx, idx + 2] = -1.0
        mat /= 12.0 * dz

        mat[0, 0:5] = np.array([-25.0, 48.0, -36.0, 16.0, -3.0]) / (12.0 * dz)
        mat[1, 0:5] = np.array([-3.0, -10.0, 18.0, -6.0, 1.0]) / (12.0 * dz)
        mat[-1, -5:] = np.array([3.0, -16.0, 36.0, -48.0, 25.0]) / (12.0 * dz)
        mat[-2, -5:] = np.array([-1.0, 6.0, -18.0, 10.0, 3.0]) / (12.0 * dz)
        return mat

    def _gnielinski(self, t_fluid, m_tube, n_tube, d_tube, fluid_number, pressure):
        fluid_name = str(fluid_number) if fluid_number == fluid_number else "Nitrate Salt"
        rho_tube = max(float(self._props.density(fluid_name, t_fluid, pressure)), 0.1)
        visc_tube = max(float(self._props.viscosity(fluid_name, t_fluid, pressure)), 1.0e-9)
        cp_tube = max(float(self._props.specheat(fluid_name, t_fluid, pressure)) * 1000.0, 50.0)
        k_film = 0.1
        pr_tube = max(visc_tube * cp_tube / k_film, 1.0e-8)

        m_dot_tube = max(m_tube / max(n_tube, 1.0), 1.0e-12)
        area = max((math.pi / 4.0) * d_tube**2, 1.0e-12)
        vel_tube = m_dot_tube / max(rho_tube * area, 1.0e-12)
        re_tube = rho_tube * vel_tube * d_tube / max(visc_tube, 1.0e-12)

        if re_tube > 4000.0:
            f = 0.035
            nu_tube = ((f / 8.0) * (re_tube - 1000.0) * pr_tube) / (
                1.0 + 12.7 * math.sqrt(f / 8.0) * (pr_tube ** (2.0 / 3.0) - 1.0)
            )
        else:
            nu_tube = 4.36

        h_tube = max(nu_tube * k_film / max(d_tube, 1.0e-9), 0.01)
        return h_tube, vel_tube

    def _zukauskas(self, fluid_shell, config, d_shell, s_t, s_l, n_tube, d_tube, t_shell, m_shell, l_shell, n_baffles, pressure):
        fluid_name = str(fluid_shell) if fluid_shell == fluid_shell else "Nitrate Salt"
        rho_shell = max(float(self._props.density(fluid_name, t_shell, pressure)), 0.1)
        visc_shell = max(float(self._props.viscosity(fluid_name, t_shell, pressure)), 1.0e-9)
        cp_shell = max(float(self._props.specheat(fluid_name, t_shell, pressure)) * 1000.0, 50.0)

        k_film = 0.5
        pr_shell = max(visc_shell * cp_shell / k_film, 1.0e-8)
        pr_film = pr_shell

        gap_baffle = l_shell / max(n_baffles + 1, 1)
        vel_shell = 3.0 * m_shell / max(rho_shell * (math.pi * d_shell * gap_baffle / 2.0), 1.0e-12)
        re_shell = rho_shell * vel_shell * d_tube / max(visc_shell, 1.0e-12)

        config_aligned = bool(round(config))
        if config_aligned:
            if re_shell < 100.0:
                c_coef, m_coef = 0.8, 0.4
            elif re_shell < 1000.0:
                c_coef, m_coef = 0.51, 0.50
            elif re_shell < 2.0e5:
                c_coef, m_coef = 0.27, 0.63
            else:
                c_coef, m_coef = 0.021, 0.84
        else:
            if re_shell < 100.0:
                c_coef, m_coef = 0.9, 0.4
            elif re_shell < 1000.0:
                c_coef, m_coef = 0.51, 0.5
            elif re_shell < 2.0e5:
                if (s_t / max(s_l, 1.0e-12)) > 2.0:
                    c_coef, m_coef = 0.4, 0.60
                else:
                    c_coef, m_coef = 0.35 * (s_t / max(s_l, 1.0e-12)) ** 0.2, 0.60
            else:
                c_coef, m_coef = 0.022, 0.84

        if re_shell > 1000.0 and n_tube < 16.0:
            if config_aligned:
                c2 = 0.6233 + 0.1007 * n_tube - 0.0095 * n_tube**2 + 2.8849e-4 * n_tube**3
            else:
                c2 = 0.5385 + 0.1292 * n_tube - 0.0123 * n_tube**2 + 3.7849e-4 * n_tube**3
        else:
            c2 = 1.0

        nu_shell = c_coef * max(re_shell, 1.0e-9) ** m_coef * pr_shell**0.36 * (pr_shell / pr_film) ** 0.25
        nu_shell *= c2
        h_shell = max(nu_shell * k_film / max(d_tube, 1.0e-9), 0.01)
        return h_shell, vel_shell

    def _zero_d_eq(self, n_nodes, h_tube, h_shell, t_shell, t_tube, k_tube, delta_x, d_in, d_out):
        a_out = d_out * math.pi * delta_x
        a_in = d_in * math.pi * delta_x

        r_conv_shell = 1.0 / np.maximum(a_out * h_shell, 1.0e-12)
        r_conv_tube = 1.0 / np.maximum(a_in * h_tube, 1.0e-12)

        r_in = d_in / 2.0
        r_out = d_out / 2.0
        r_cond = math.log(max(r_out / max(r_in, 1.0e-12), 1.0 + 1.0e-12)) / (
            2.0 * math.pi * max(k_tube, 1.0e-12) * max(delta_x, 1.0e-12)
        )

        r_fouling_out = (0.0005 * 0.3048**2 * 3.41 * 0.5556) / max(a_out, 1.0e-12)
        r_fouling_in = (0.0005 * 0.3048**2 * 3.41 * 0.5556) / max(a_in, 1.0e-12)

        q_htf_tube = np.zeros(n_nodes, dtype=float)
        q_htf_shell = np.zeros(n_nodes, dtype=float)

        for idx in range(n_nodes):
            r_idx = n_nodes - 1 - idx
            denom = r_conv_shell[idx] + r_conv_tube[r_idx] + r_cond + r_fouling_in + r_fouling_out
            q = (t_shell[idx] - t_tube[r_idx]) / max(denom, 1.0e-12)
            q_htf_shell[idx] = q
            q_htf_tube[r_idx] = q

        return q_htf_tube, q_htf_shell

    def _shell_subroutine(self, config, shell_p, t_env, p_in_shell, gap_baffle, n_nodes, k_steel, k_ins, s_t, s_l,
                          m_shell, n_baffles, n_hex, th_baffle, th_ins, t_shell, d_shell, l_baffle, a_baffle,
                          shell_passes, l_shell, n_tubes, d_tube, q_htf_shell, v_out, delta_t, fluid_shell):
        v_tubes = ((math.pi * d_tube**2) / 4.0) * gap_baffle * n_tubes
        vc = ((math.pi * d_shell**2) / 4.0) * gap_baffle / max(shell_p, 1) - v_tubes
        vc = max(vc, 1.0e-9)

        g = 9.8
        n_regions = 2 * shell_passes - 1
        u_region = 1

        cp_shell = np.zeros(n_nodes, dtype=float)
        rho_shell = np.zeros(n_nodes, dtype=float)
        mu_shell = np.zeros(n_nodes, dtype=float)
        k_shell = np.full(n_nodes, 0.5, dtype=float)
        pr_shell = np.zeros(n_nodes, dtype=float)

        tw = np.full(n_nodes, t_env + 5.0, dtype=float)
        tww = tw.copy()

        for idx in range(n_nodes):
            fluid_name = str(fluid_shell) if fluid_shell == fluid_shell else "Nitrate Salt"
            cp_shell[idx] = max(float(self._props.specheat(fluid_name, t_shell[idx], p_in_shell)) * 1000.0, 50.0)
            rho_shell[idx] = max(float(self._props.density(fluid_name, t_shell[idx], p_in_shell)), 0.1)
            mu_shell[idx] = max(float(self._props.viscosity(fluid_name, t_shell[idx], p_in_shell)), 1.0e-9)
            pr_shell[idx] = cp_shell[idx] * mu_shell[idx] / max(k_shell[idx], 1.0e-12)

        vel_shell = 2.0 * m_shell / np.maximum(rho_shell * (d_shell * gap_baffle), 1.0e-12)

        def baffle_resistance(k_idx, k1_idx):
            re_k = rho_shell[k_idx] * l_baffle * vel_shell[k_idx] / max(mu_shell[k_idx], 1.0e-12)
            nu_plate_k = 0.0296 * max(re_k, 1.0e-12) ** 0.8 * max(pr_shell[k_idx], 1.0e-12) ** 0.333
            h_k = nu_plate_k * k_shell[k_idx] / max(l_baffle, 1.0e-12)
            r1conv_i = 1.0 / max(a_baffle * h_k, 1.0e-12)

            if k1_idx is None:
                r1conv_i1 = 1.0e9
            else:
                re_k1 = rho_shell[k1_idx] * l_baffle * vel_shell[k1_idx] / max(mu_shell[k1_idx], 1.0e-12)
                nu_plate_k1 = 0.0296 * max(re_k1, 1.0e-12) ** 0.8 * max(pr_shell[k1_idx], 1.0e-12) ** 0.333
                h_k1 = nu_plate_k1 * k_shell[k1_idx] / max(l_baffle, 1.0e-12)
                r1conv_i1 = 1.0 / max(a_baffle * h_k1, 1.0e-12)

            r1cond_i = th_baffle / max(a_baffle * k_steel, 1.0e-12)
            return r1cond_i + r1conv_i1 + r1conv_i

        rt_hl1 = baffle_resistance(0, 1 if n_nodes > 1 else None)
        t_shell_new = np.array(t_shell, dtype=float)

        error = 1.0
        while error > 1.0e-5:
            k_idx = 1
            for _ in range(n_regions):
                if u_region > 0:
                    vol = vc
                    a_gap = math.pi * d_shell * gap_baffle / 2.0
                    nodes_region = max(n_baffles - 1, 0)
                else:
                    vol = 2.0 * vc
                    a_gap = math.pi * d_shell * gap_baffle
                    nodes_region = 1

                for _ in range(nodes_region):
                    if k_idx >= n_nodes:
                        break

                    rt_hl4 = rt_hl1
                    rt_hl1 = baffle_resistance(k_idx, k_idx + 1 if k_idx < (n_nodes - 1) else None)

                    mu_air = 184.6e-7
                    k_air = 26.3e-3
                    rho_air = 1.1614
                    visc_air = 15.89e-6
                    alpha_air = 22.5e-6
                    beta_air = 1.0 / max(t_env, 1.0)

                    r_i = th_ins / max(a_gap * k_ins, 1.0e-12)

                    emis = 0.9
                    sigma = 5.67e-8
                    h_rad = emis * sigma * (t_env**2 + tw[k_idx] ** 2) * (t_env + tw[k_idx])
                    r_ii = 1.0 / max(a_gap * h_rad, 1.0e-12)

                    pr_ext = 0.7
                    if v_out < 0.1:
                        ra_ext = beta_air * (tw[k_idx] - t_env) * g * d_shell**3 / max(visc_air * alpha_air, 1.0e-12)
                        nu_ext = (0.60 + 0.387 * max(ra_ext, 0.0) ** 0.16 / (1.0 + (0.559 / pr_ext) ** (9.0 / 16.0)) ** 0.2963) ** 2
                    else:
                        re_ext = rho_air * v_out * d_shell / max(mu_air, 1.0e-12)
                        nu_ext = (
                            0.3
                            + (0.62 * max(re_ext, 0.0) ** 0.5 * pr_ext**0.333) / (1.0 + (0.4 / pr_ext) ** 0.666) ** 0.25
                        ) * (1.0 + (max(re_ext, 0.0) / 282000.0) ** 0.625) ** 0.8

                    h_ext = k_air * nu_ext / max(d_shell, 1.0e-12)
                    r_iii = 1.0 / max(a_gap * h_ext, 1.0e-12)

                    rt_parallel = (r_ii * r_iii) / max(r_iii + r_ii, 1.0e-12)
                    rt_hl2 = rt_parallel + r_i

                    hl_3 = q_htf_shell[k_idx] * n_tubes

                    km1 = max(k_idx - 1, 0)
                    kp1 = min(k_idx + 1, n_nodes - 1)
                    num = (
                        rho_shell[k_idx] * vol * cp_shell[k_idx] * t_shell[k_idx] / max(delta_t, 1.0e-12)
                        + m_shell * cp_shell[km1] * t_shell_new[km1]
                        - hl_3
                        + t_shell[kp1] / max(rt_hl1, 1.0e-12)
                        + t_env / max(rt_hl2, 1.0e-12)
                        + t_shell_new[km1] / max(rt_hl4, 1.0e-12)
                    )
                    den = (
                        rho_shell[k_idx] * vol * cp_shell[k_idx] / max(delta_t, 1.0e-12)
                        + m_shell * cp_shell[k_idx]
                        + 1.0 / max(rt_hl1, 1.0e-12)
                        + 1.0 / max(rt_hl2, 1.0e-12)
                        + 1.0 / max(rt_hl4, 1.0e-12)
                    )
                    t_shell_new[k_idx] = num / max(den, 1.0e-12)

                    q_loss = (t_shell_new[k_idx] - t_env) / max(rt_hl2, 1.0e-12)
                    tww[k_idx] = q_loss * rt_parallel + t_env
                    k_idx += 1

                u_region = -u_region

            error = 1.0e-9
            tw = tww.copy()
            u_region = 1

        if n_nodes > 2:
            t_shell_new[-1] = t_shell_new[-2] - (t_shell_new[-3] - t_shell_new[-2])

        h_shell = np.zeros(n_nodes, dtype=float)
        vel_shell_out = np.zeros(n_nodes, dtype=float)
        for idx in range(n_nodes):
            h_shell_i, vel_shell_i = self._zukauskas(
                fluid_shell=fluid_shell,
                config=config,
                d_shell=d_shell,
                s_t=s_t,
                s_l=s_l,
                n_tube=n_tubes,
                d_tube=d_tube,
                t_shell=t_shell_new[idx],
                m_shell=m_shell,
                l_shell=l_shell,
                n_baffles=n_baffles,
                pressure=p_in_shell,
            )
            h_shell[idx] = h_shell_i
            vel_shell_out[idx] = vel_shell_i

        p_shell_out = p_in_shell
        return t_shell_new, h_shell, vel_shell_out, p_shell_out

    def _tube_one_d_inc(self, q_htf_tube, t_tube, n_nodes, gap_baffle, d_in, p_tube_in, m_tube, n_tube, fluid_tube, delta_t):
        j = n_nodes
        d = d_in
        rh = d / 4.0
        lambda_f = 0.003

        dz = gap_baffle if n_nodes <= 1 else gap_baffle
        d1_d1z = self._build_d1dz(n_nodes=j, dz=max(dz, 1.0e-12))

        t_in_c = t_tube[0] - 273.15
        h_in = self._h_from_tube_temp_c(t_in_c)
        rho_in = 1.0 / max(self._vspec_from_h(h_in), 1.0e-12)
        p_in = p_tube_in
        area = max(math.pi * d**2 / 4.0, 1.0e-12)
        v_in = m_tube / max(rho_in * area, 1.0e-12)

        q_term = q_htf_tube / max(dz * area, 1.0e-12)

        h_old = self._h_from_tube_temp_c(t_tube - 273.15)
        rho_old = 1.0 / np.maximum(self._vspec_from_h(h_old), 1.0e-12)
        p_old = np.full(j, p_tube_in, dtype=float)
        v_old = m_tube / np.maximum(rho_old * area, 1.0e-12)

        h = h_old.copy()
        v = v_old.copy()
        p = p_old.copy()

        lim = 1.0
        max_iter = 30
        iter_count = 0

        while lim > 1.0e-5 and iter_count < max_iter:
            iter_count += 1

            dv_dt = (v - v_old) / max(delta_t, 1.0e-12)
            dh_dt = (h - h_old) / max(delta_t, 1.0e-12)

            d_v_dh = self._dv_dh(h)
            v_spec = self._vspec_from_h(h)
            rho = 1.0 / np.maximum(v_spec, 1.0e-12)
            d_rho_dh = -(rho**2) * d_v_dh

            dh_dz = d1_d1z @ h
            dv_dz = d1_d1z @ v

            f1 = d_rho_dh * dh_dt + rho * (d1_d1z @ v) + v * d_rho_dh * dh_dz
            j11 = self._diag(rho) @ d1_d1z + self._diag(d_rho_dh * dh_dz)
            j13 = self._diag(d_rho_dh) / max(delta_t, 1.0e-12) + self._diag(d_rho_dh * dv_dz) + self._diag(v * d_rho_dh) @ d1_d1z

            dp_dz = d1_d1z @ p
            f2 = v_spec * dp_dz + dv_dt + v * dv_dz + lambda_f * (v**2) / max(8.0 * rh, 1.0e-12)

            j23 = self._diag(d_v_dh * dp_dz)
            j21 = np.eye(j) / max(delta_t, 1.0e-12)
            for ii in range(j):
                j21[ii, ii] += dv_dz[ii] + 2.0 * lambda_f * v[ii] / max(8.0 * rh, 1.0e-12)
            j21 += self._diag(v) @ d1_d1z
            j22 = self._diag(v_spec) @ d1_d1z

            dp_dt = (p - p_old) / max(delta_t, 1.0e-12)
            f3 = dh_dt + v * dh_dz + v_spec * (-q_term - dp_dt)

            j31 = self._diag(dh_dz)
            j32 = -self._diag(v_spec / max(delta_t, 1.0e-12))
            j33 = np.eye(j) / max(delta_t, 1.0e-12) + self._diag(v) @ d1_d1z + self._diag(d_v_dh * (-q_term - dp_dt))

            j11[0, :] = 0.0
            j11[0, 0] = 1.0
            j13[0, :] = 0.0
            f1[0] = v[0] - v_in

            j21[0, :] = 0.0
            j22[0, :] = 0.0
            j22[0, 0] = 1.0
            j23[0, :] = 0.0
            f2[0] = p[0] - p_in

            j31[0, :] = 0.0
            j32[0, :] = 0.0
            j33[0, :] = 0.0
            j33[0, 0] = 1.0
            f3[0] = h[0] - h_in

            f = np.concatenate([f1, f2, f3])
            jacob = np.zeros((3 * j, 3 * j), dtype=float)

            jacob[0:j, 0:j] = j11
            jacob[0:j, 2 * j:3 * j] = j13
            jacob[j:2 * j, 0:j] = j21
            jacob[j:2 * j, j:2 * j] = j22
            jacob[j:2 * j, 2 * j:3 * j] = j23
            jacob[2 * j:3 * j, 0:j] = j31
            jacob[2 * j:3 * j, j:2 * j] = j32
            jacob[2 * j:3 * j, 2 * j:3 * j] = j33

            try:
                dx = np.linalg.solve(jacob, f)
            except np.linalg.LinAlgError:
                dx = np.linalg.lstsq(jacob, f, rcond=None)[0]

            v = v - dx[0:j]
            p = p - dx[j:2 * j]
            h = h - dx[2 * j:3 * j]

            lim = math.sqrt(float(np.max(dx**2)))

        t_tube_new = self._tube_temp_from_h(h)

        h_tube = np.array(t_tube_new, dtype=float)
        vel_tube = np.zeros(j, dtype=float)
        for idx in range(j):
            h_i, vel_i = self._gnielinski(
                t_fluid=t_tube_new[idx],
                m_tube=m_tube,
                n_tube=n_tube,
                d_tube=d_in,
                fluid_number=fluid_tube,
                pressure=p[idx],
            )
            h_tube[idx] = h_i
            vel_tube[idx] = vel_i

        p_tube_out = p[-1]
        return t_tube_new, h_tube, vel_tube, p_tube_out

    def _initialize_profiles(self):
        n_nodes = self._node_count()
        self._ensure_state(n_nodes)

        p_in_shell = self.p_shell_in_c.v if self.p_shell_in_c.v > 0.0 else 200000.0
        p_in_tube = self.p_tube_in_c.v if self.p_tube_in_c.v > 0.0 else 1000000.0
        t_tube_in = self.t_tube_in_c.v if self.t_tube_in_c.v > 350.0 else 650.0
        t_shell_in = self.t_shell_in_c.v if self.t_shell_in_c.v > 200.0 else 400.0

        m_tube = max(self._safe(self.m_dot_tube_c.v, 0.0), 0.0)
        m_shell = max(self._safe(self.m_dot_shell_c.v, 0.0), 0.0)
        cp_shell = max(float(self._props.specheat(str(self.fluid_shell.v), t_shell_in, p_in_shell)) * 1000.0, 50.0)
        cp_tube = max(float(self._props.specheat(str(self.fluid_tube.v), t_tube_in, p_in_tube)) * 1000.0, 50.0)

        c_shell = max(m_shell * cp_shell, 1.0)
        c_tube = max(m_tube * cp_tube, 1.0)
        c_min = min(c_shell, c_tube)

        l_shell = max(self._safe(self.l_shell.v, 1.0), 1.0e-6)
        d_shell = max(self._safe(self.d_shell.v, 1.0), 1.0e-6)
        k_steel = max(self._safe(self.k_steel.v, 15.0), 1.0e-6)
        n_tubes = max(int(round(self._safe(self.n_tubes.v, 1.0))), 1)
        d_out = max(2.0 * self._safe(self.r_out.v, 0.01), 1.0e-5)
        area = math.pi * d_out * l_shell * n_tubes
        ua = max(k_steel * area / max(d_shell, 1.0e-6), 25.0)
        ntu = ua / c_min
        eff = 1.0 - math.exp(-ntu)

        q_dot = eff * c_min * (t_shell_in - t_tube_in)
        t_shell_out = t_shell_in - q_dot / c_shell
        t_tube_out = t_tube_in + q_dot / c_tube

        self._t_shell_nodes = np.linspace(t_shell_in, t_shell_out, n_nodes)
        self._t_tube_nodes = np.linspace(t_tube_in, t_tube_out, n_nodes)
        self._h_shell_nodes = self._t_shell_nodes * cp_shell
        self._h_tube_nodes = self._t_tube_nodes * cp_tube
        self._mode_prev = 1.0
        self._initialized = True

        self._set_outputs(
            mode_iter=1.0,
            t_shell_out=t_shell_out,
            p_shell_out=self.p_shell_in_c.v,
            t_tube_out=t_tube_out,
            p_tube_out=self.p_tube_in_c.v,
        )

    def _set_outputs(self, mode_iter, t_shell_out, p_shell_out, t_tube_out, p_tube_out):
        p_shell_out_c = p_shell_out
        p_tube_out_c = p_tube_out
        p_shell_out_d = p_shell_out
        p_tube_out_d = p_tube_out

        if mode_iter == 1.0:
            t_shell_out_c = t_shell_out
            if self.p_shell_in_c.v < 0.0:
                p_shell_out_c = self.p_shell_in_c.v

            t_tube_out_c = t_tube_out
            if self.p_tube_in_c.v < 0.0:
                p_tube_out_c = self.p_tube_in_c.v

            t_shell_out_d = self.t_shell_in_d.v
            p_shell_out_d = self.p_shell_in_d.v
            t_tube_out_d = self.t_tube_in_d.v
            p_tube_out_d = self.p_tube_in_d.v
        else:
            t_shell_out_d = t_shell_out
            if self.p_shell_in_d.v < 0.0:
                p_shell_out_d = self.p_shell_in_d.v

            t_tube_out_d = t_tube_out
            if self.p_tube_in_d.v < 0.0:
                p_tube_out_d = self.p_tube_in_d.v

            t_shell_out_c = self.t_shell_in_c.v
            p_shell_out_c = self.p_shell_in_c.v
            t_tube_out_c = self.t_tube_in_c.v
            p_tube_out_c = self.p_tube_in_c.v

        self.mode.v = mode_iter
        self.m_dot_shell_out_c.v = self.m_dot_shell_c.v
        self.t_shell_out_c.v = t_shell_out_c
        self.p_shell_out_c.v = p_shell_out_c
        self.m_dot_shell_out_d.v = self.m_dot_shell_d.v
        self.t_shell_out_d.v = t_shell_out_d
        self.p_shell_out_d.v = p_shell_out_d
        self.m_dot_tube_out_c.v = self.m_dot_tube_c.v
        self.t_tube_out_c.v = t_tube_out_c
        self.p_tube_out_c.v = p_tube_out_c
        self.m_dot_tube_out_d.v = self.m_dot_tube_d.v
        self.t_tube_out_d.v = t_tube_out_d
        self.p_tube_out_d.v = p_tube_out_d
        self.mass_counter_out_c.v = self.mass_counter_c.v
        self.mass_counter_out_d.v = self.mass_counter_d.v

    def calculate(self):
        if self.model.is_first_step:
            if not self._initialized:
                self._initialize_profiles()
            else:
                self._set_outputs(
                    mode_iter=self._mode_prev,
                    t_shell_out=float(self._t_shell_nodes[-1]),
                    p_shell_out=self.p_shell_in_c.v if self._mode_prev == 1.0 else self.p_shell_in_d.v,
                    t_tube_out=float(self._t_tube_nodes[-1]),
                    p_tube_out=self.p_tube_in_c.v if self._mode_prev == 1.0 else self.p_tube_in_d.v,
                )
            return

        n_nodes = self._node_count()
        self._ensure_state(n_nodes)

        shell_passes_total = max(int(round(self.shell_passes.v * self.n_hex.v)), 1)

        if self.m_dot_tube_c.v > self.m_dot_tube_d.v:
            mode_iter = 1.0
            m_shell = self.m_dot_shell_c.v if self.m_dot_shell_c.v > 0.0 else 1.0
            t_shell_in = self.t_shell_in_c.v if self.t_shell_in_c.v > 0.0 else 600.0
            p_in_shell = self.p_shell_in_c.v if self.p_shell_in_c.v > 0.0 else 101325.0
            p_in_tube = self.p_tube_in_c.v if self.p_tube_in_c.v > 0.0 else 101325.0
            m_tubes = self.m_dot_tube_c.v
            t_tube_in = self.t_tube_in_c.v
        else:
            mode_iter = 0.0
            m_shell = self.m_dot_shell_d.v if self.m_dot_shell_d.v > 0.0 else 1.0
            t_shell_in = self.t_shell_in_d.v if self.t_shell_in_d.v > 0.0 else 600.0
            p_in_shell = self.p_shell_in_d.v if self.p_shell_in_d.v > 0.0 else 101325.0
            p_in_tube = self.p_tube_in_d.v if self.p_tube_in_d.v > 0.0 else 101325.0
            m_tubes = self.m_dot_tube_d.v
            t_tube_in = self.t_tube_in_d.v

        if mode_iter != self._mode_prev:
            self._t_shell_nodes = self._t_shell_nodes[::-1].copy()
            self._t_tube_nodes = self._t_tube_nodes[::-1].copy()
            self._h_shell_nodes = self._h_shell_nodes[::-1].copy()
            self._h_tube_nodes = self._h_tube_nodes[::-1].copy()

        gap_baffle = self.l_shell.v / max(self.n_baffles.v + 1.0, 1.0)
        n_tubes = max(int(round(self._safe(self.n_tubes.v, 1.0))), 1)
        m_tube = max(m_tubes / n_tubes, 1.0e-9)
        d_in = max(2.0 * self._safe(self.r_in.v, 0.01), 1.0e-5)
        d_out = max(2.0 * self._safe(self.r_out.v, 0.012), d_in + 1.0e-6)

        delta_t = max(float(self._safe(self.model.settings.timestep, 1.0)), 1.0e-9)

        q_htf_tube, q_htf_shell = self._zero_d_eq(
            n_nodes=n_nodes,
            h_tube=self._h_tube_nodes,
            h_shell=self._h_shell_nodes,
            t_shell=self._t_shell_nodes,
            t_tube=self._t_tube_nodes,
            k_tube=self.k_steel.v,
            delta_x=gap_baffle,
            d_in=d_in,
            d_out=d_out,
        )

        a_baffle = (2.0 / 5.0) * math.pi * self.d_shell.v**2 / 4.0

        t_shell_new, h_shell, _vel_shell, p_shell_out = self._shell_subroutine(
            config=self.config.v,
            shell_p=max(int(round(self.shell_passes.v)), 1),
            t_env=self.t_env.v,
            p_in_shell=p_in_shell,
            gap_baffle=gap_baffle,
            n_nodes=n_nodes,
            k_steel=max(self.k_steel.v, 1.0e-9),
            k_ins=max(self.k_ins.v, 1.0e-9),
            s_t=max(self.s_t.v, 1.0e-9),
            s_l=max(self.s_l.v, 1.0e-9),
            m_shell=m_shell,
            n_baffles=max(int(round(self.n_baffles.v)), 1),
            n_hex=max(int(round(self.n_hex.v)), 1),
            th_baffle=max(self.th_baffle.v, 1.0e-9),
            th_ins=max(self.th_ins.v, 0.0),
            t_shell=self._t_shell_nodes,
            d_shell=max(self.d_shell.v, 1.0e-6),
            l_baffle=max(self.l_baffle.v, 1.0e-6),
            a_baffle=max(a_baffle, 1.0e-9),
            shell_passes=shell_passes_total,
            l_shell=max(self.l_shell.v, 1.0e-6),
            n_tubes=max(self.n_tubes.v, 1.0),
            d_tube=d_out,
            q_htf_shell=q_htf_shell,
            v_out=max(self.v_out.v, 0.0),
            delta_t=delta_t,
            fluid_shell=self.fluid_shell.v,
        )

        t_tube_new, h_tube, _vel_tube, p_tube_out = self._tube_one_d_inc(
            q_htf_tube=q_htf_tube,
            t_tube=self._t_tube_nodes,
            n_nodes=n_nodes,
            gap_baffle=gap_baffle,
            d_in=d_in,
            p_tube_in=p_in_tube,
            m_tube=max(m_tube, 1.0e-9),
            n_tube=max(self.n_tubes.v, 1.0),
            fluid_tube=self.fluid_tube.v,
            delta_t=delta_t,
        )

        t_shell_out = float(t_shell_new[-1])
        t_tube_out = float(t_tube_new[-1])

        self._set_outputs(mode_iter, t_shell_out, p_shell_out, t_tube_out, p_tube_out)

        if self.model.is_converged:
            self._t_shell_nodes = t_shell_new
            self._t_tube_nodes = t_tube_new
            self._h_shell_nodes = h_shell
            self._h_tube_nodes = h_tube
            self._mode_prev = mode_iter
