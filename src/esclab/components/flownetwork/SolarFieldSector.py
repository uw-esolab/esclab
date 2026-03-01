"""Type 4034 solar-field sector component converted from Fortran."""

import math
import re
from pathlib import Path

import numpy as np
from eeslib import fluid_properties as fp

from esclab.components.esol_properties import Incompressible
from esclab.components.flownetwork.SimplePipe import FricFactor_IC
from esclab.simulate import Component


class SolarFieldSector(Component):
    """
    TRNSYS Type 4034: ESOL4034-SolarFieldSector.

    Parameters
    ----------
    eta_defocus_1, eta_defocus_2, eta_defocus_3, eta_defocus_4 : float
        Defocusing setpoints from the original Fortran type.
    iam_a0, iam_a1, iam_a2 : float
        Incidence-angle modifier coefficients.
    eta_tracking, eta_soil, eta_reflect, sf_avail : float
        Optical efficiency factors.
    distance_sca, ave_focal_length, w_ap, row_distance, l_exp_loop : float
        Solar-field geometry terms.
    n_sca, n_loop : float
        Number of SCA elements and loops in the sector.
    fluid_id : float|str
        HTF identifier.
    mc_receiver_mult, mc_header_mult : float
        Thermal mass multipliers.
    t_nominal, time_lim_df : float
        Defocusing/tracking reference temperatures and time limit.
    sf_label, max_loop, n_sectors : float
        Sector indexing and global solar-field dimensions.
    l_tot, n_nodes_per_loop : float
        Receiver loop length and discretization.
    t_init_sf, t_init_in_header, t_init_return_header : float
        Initial temperatures for loop and header states [K].
    d_receiver, roughness_pipe : float
        Pipe geometry and roughness.
    prob_pristine, prob_broken_glass, prob_lost_vacuum : float
        Degradation state probabilities used by loop-file generation in Fortran.
    m_dot_std, n_h2, new_loops : float
        Statistical loop variability and hydrogen-state controls.
    inlet_header_heat_loss, return_header_heat_loss : float
        Header heat-loss multipliers used in Fortran RK4 equations.
    geom_file, nn_base : str
        Fortran label equivalents for header geometry data and NN base path.

    Inputs
    ------
    mass_flow, pressure, temperature, control_signal, mass_counter : float
        Inlet hydraulic and inventory states.
    ani, t_amb, t_sky, wind, theta, phi, t_tracking : float
        Weather and sun-angle conditions.

    Outputs
    -------
    mass_flow_out, pressure_out, temperature_out, t4_ave, mass_counter_out : float
        Main outlet state in Fortran output order 1..5.
    defocusing : float
        Aggregate defocusing metric (Fortran output 6).
    defocus_group_1..defocus_group_8 : float
        Grouped defocusing diagnostics (Fortran outputs 7..14).
    temperature_group_1..temperature_group_8 : float
        Grouped temperature diagnostics (Fortran outputs 15..22).

    Notes
    -----
    This pass maps the Type 4034 interface, simulation-phase structure, and
    output ordering while using a reduced thermal/hydraulic surrogate.
    """

    eta_defocus_1 = Component.Parameter()
    eta_defocus_2 = Component.Parameter()
    eta_defocus_3 = Component.Parameter()
    eta_defocus_4 = Component.Parameter()
    iam_a0 = Component.Parameter()
    iam_a1 = Component.Parameter()
    iam_a2 = Component.Parameter()
    eta_tracking = Component.Parameter()
    eta_soil = Component.Parameter()
    eta_reflect = Component.Parameter()
    sf_avail = Component.Parameter()
    distance_sca = Component.Parameter()
    ave_focal_length = Component.Parameter()
    w_ap = Component.Parameter()
    row_distance = Component.Parameter()
    l_exp_loop = Component.Parameter()
    n_sca = Component.Parameter()
    n_loop = Component.Parameter()
    fluid_id = Component.Parameter()
    mc_receiver_mult = Component.Parameter()
    mc_header_mult = Component.Parameter()
    t_nominal = Component.Parameter()
    time_lim_df = Component.Parameter()
    sf_label = Component.Parameter()
    max_loop = Component.Parameter()
    n_sectors = Component.Parameter()
    l_tot = Component.Parameter()
    n_nodes_per_loop = Component.Parameter()
    t_init_sf = Component.Parameter()
    t_init_in_header = Component.Parameter()
    t_init_return_header = Component.Parameter()
    d_receiver = Component.Parameter()
    roughness_pipe = Component.Parameter()
    prob_pristine = Component.Parameter()
    prob_broken_glass = Component.Parameter()
    prob_lost_vacuum = Component.Parameter()
    m_dot_std = Component.Parameter()
    n_h2 = Component.Parameter()
    new_loops = Component.Parameter()
    inlet_header_heat_loss = Component.Parameter()
    return_header_heat_loss = Component.Parameter()
    geom_file = Component.Parameter("")
    nn_base = Component.Parameter("")

    mass_flow = Component.Input()
    pressure = Component.Input()
    temperature = Component.Input()
    control_signal = Component.Input()
    mass_counter = Component.Input()
    ani = Component.Input()
    t_amb = Component.Input()
    t_sky = Component.Input()
    wind = Component.Input()
    theta = Component.Input()
    phi = Component.Input()
    t_tracking = Component.Input()

    mass_flow_out = Component.Output()
    pressure_out = Component.Output()
    temperature_out = Component.Output()
    t4_ave = Component.Output()
    mass_counter_out = Component.Output()
    defocusing = Component.Output()
    defocus_group_1 = Component.Output()
    defocus_group_2 = Component.Output()
    defocus_group_3 = Component.Output()
    defocus_group_4 = Component.Output()
    defocus_group_5 = Component.Output()
    defocus_group_6 = Component.Output()
    defocus_group_7 = Component.Output()
    defocus_group_8 = Component.Output()
    temperature_group_1 = Component.Output()
    temperature_group_2 = Component.Output()
    temperature_group_3 = Component.Output()
    temperature_group_4 = Component.Output()
    temperature_group_5 = Component.Output()
    temperature_group_6 = Component.Output()
    temperature_group_7 = Component.Output()
    temperature_group_8 = Component.Output()

    _props = Incompressible()
    _is_initialized = False
    _n_loop_i = 1
    _n_nodes_i = 2
    _vol_loop_cv = 1.0
    _loop_t = np.array([[300.0, 300.0]])
    _t_header_inlet = np.array([300.0, 300.0])
    _t_header_return = np.array([300.0, 300.0])
    _n_cv_header = 1
    _n_node_header = 2
    _inds_header_in = np.array([1])
    _inds_left = np.array([1])
    _inds_right = np.array([2])
    _l_cv_inlet = np.array([1.0])
    _l_cv_return = np.array([1.0, 1.0])
    _vol_inlet = np.array([1.0])
    _vol_return = np.array([1.0, 1.0])
    _d_inlet = np.array([0.1])
    _d_return = np.array([0.1, 0.1])
    _m_dots_in = np.array([0.0])
    _m_dots_return = np.array([0.0, 0.0])
    _m_dot_var = np.array([1.0])
    _defocus_mode = np.array([1.0])
    _defocus_groups = np.zeros(8)
    _temperature_groups = np.zeros(8)
    _t4_ave_hold = 300.0
    _mass_htf_hold = 0.0
    _ff_guess = 0.1
    _time_df = np.array([0.0])
    _t_bar_sf = np.array([300.0])
    _ncv_state = np.zeros((4, 1), dtype=int)
    _state_indices = [[np.array([], dtype=int)] for _ in range(4)]
    _h2_pressure = [np.array([], dtype=float)]
    _features = np.zeros((1, 7), dtype=float)
    _nn_layers = {}
    _nn_loaded = False
    _minmax2 = np.array(
        [
            [0.0, 0.0, 0.0, 0.0],
            [1.0, 1.0, 1.0, 1.0],
            [300.0, 300.0, 300.0, 300.0],
            [800.0, 800.0, 800.0, 800.0],
            [0.0, 0.0, 0.0, 0.0],
            [2.0e6, 2.0e6, 2.0e6, 2.0e6],
            [250.0, 250.0, 250.0, 250.0],
            [340.0, 340.0, 340.0, 340.0],
            [220.0, 220.0, 220.0, 220.0],
            [340.0, 340.0, 340.0, 340.0],
            [0.0, 0.0, 0.0, 0.0],
            [30.0, 30.0, 30.0, 30.0],
            [0.0, 0.0, 0.0, 0.0],
            [1500.0, 1500.0, 1500.0, 1500.0],
            [-1000.0, -1000.0, -1000.0, -1000.0],
            [1000.0, 1000.0, 1000.0, 1000.0],
        ],
        dtype=float,
    )

    @staticmethod
    def _debugging_root():
        return Path(__file__).resolve().parent / "fortran-source" / "debugging"

    @staticmethod
    def _safe(value, default=0.0):
        return value if value == value else default

    @staticmethod
    def _clamp(value, low=0.0, high=1.0):
        return max(low, min(high, value))

    def _fluid_density(self, fluid_name, temperature_k):
        name = str(fluid_name) if fluid_name == fluid_name else "Nitrate Salt"
        try:
            return max(self._props.density(name, temperature_k), 1.0)
        except Exception:
            return max(self._props.density("Nitrate Salt", temperature_k), 1.0)

    def _fluid_cp(self, fluid_name, temperature_k):
        name = str(fluid_name) if fluid_name == fluid_name else "Nitrate Salt"
        try:
            return max(self._props.specheat(name, temperature_k) * 1000.0, 1000.0)
        except Exception:
            return max(self._props.specheat("Nitrate Salt", temperature_k) * 1000.0, 1000.0)

    @staticmethod
    def _h_dowtherm_a(temperature_k):
        td = temperature_k - 273.15
        return (-12.7078 + 1.481714 * td + 0.0014292857 * td**2) * 1000.0

    def _fluid_viscosity(self, fluid_name, temperature_k, pressure_pa):
        name = str(fluid_name) if fluid_name == fluid_name else "Nitrate Salt"
        try:
            mu = fp.viscosity(name, T=temperature_k, P=pressure_pa)
            return max(float(mu), 1.0e-6)
        except Exception:
            try:
                mu = self._props.viscosity(name, temperature_k)
                return max(float(mu), 1.0e-6)
            except Exception:
                return 2.5e-3

    def _fluid_density_ees(self, fluid_name, temperature_k, pressure_pa):
        name = str(fluid_name) if fluid_name == fluid_name else "Nitrate Salt"
        try:
            rho = fp.density(name, T=temperature_k, P=pressure_pa)
            return max(float(rho), 1.0)
        except Exception:
            return self._fluid_density(name, temperature_k)

    def _dp_segment(
        self,
        mass_flow,
        diameter,
        length,
        temperature_k,
        pressure_pa,
        n_exp=0.0,
        n_con=0.0,
        n_els=0.0,
        n_elm=0.0,
        n_ell=0.0,
        n_gav=0.0,
        n_glv=0.0,
        n_chv=0.0,
        n_lw=0.0,
        n_lcv=0.0,
        n_bja=0.0,
    ):
        if mass_flow <= 0.0:
            return 0.0

        d = max(diameter, 1.0e-5)
        l = max(length, 0.0)
        rho = self._fluid_density_ees(self.fluid_id.v, temperature_k, pressure_pa)
        mu = self._fluid_viscosity(self.fluid_id.v, temperature_k, pressure_pa)
        area = math.pi * d**2 / 4.0
        vel = mass_flow / max(rho * area, 1.0e-9)
        re = abs(rho * vel * d / max(mu, 1.0e-9))

        self._ff_guess = max(self._ff_guess, 0.05)
        rel_rough = max(self._safe(self.roughness_pipe.v, 1.0e-5), 1.0e-8) / d
        ff = FricFactor_IC(rel_rough, max(re, 1.0), self._ff_guess)
        if ff is None:
            ff = self._ff_guess
        try:
            ff = float(ff)
        except Exception:
            ff = self._ff_guess
        if not math.isfinite(ff):
            ff = self._ff_guess
        ff = max(ff, 0.01)
        self._ff_guess = ff

        g = 9.81
        hl_pm = ff * vel * vel / (2.0 * d * g)
        dp_pipe = hl_pm * rho * g * l

        dp_exp = 0.25 * rho * vel * vel * max(n_exp, 0.0)
        dp_con = 0.25 * rho * vel * vel * max(n_con, 0.0)
        d_over_f_hl = (d / max(ff, 1.0e-9)) * hl_pm * rho * g
        dp_els = 0.9 * d_over_f_hl * max(n_els, 0.0)
        dp_elm = 0.75 * d_over_f_hl * max(n_elm, 0.0)
        dp_ell = 0.6 * d_over_f_hl * max(n_ell, 0.0)
        dp_gav = 0.19 * d_over_f_hl * max(n_gav, 0.0)
        dp_glv = 10.0 * d_over_f_hl * max(n_glv, 0.0)
        dp_chv = 2.5 * d_over_f_hl * max(n_chv, 0.0)
        dp_lw = 1.8 * d_over_f_hl * max(n_lw, 0.0)
        dp_lcv = 10.0 * d_over_f_hl * max(n_lcv, 0.0)
        dp_bja = 8.69 * d_over_f_hl * max(n_bja, 0.0)

        return max(dp_pipe + dp_exp + dp_con + dp_els + dp_elm + dp_ell + dp_gav + dp_glv + dp_chv + dp_lw + dp_lcv + dp_bja, 0.0)

    def _group_index(self, loop_number):
        # Type4034 Fortran grouping: loops 1-13, 14-26, ..., 92-104.
        if loop_number <= 13:
            return 0
        if loop_number <= 26:
            return 1
        if loop_number <= 39:
            return 2
        if loop_number <= 52:
            return 3
        if loop_number <= 65:
            return 4
        if loop_number <= 78:
            return 5
        if loop_number <= 91:
            return 6
        return 7

    @staticmethod
    def _row_shadow(phi, row_distance, w_ap):
        eta_row = abs(math.cos(phi)) * row_distance / max(w_ap, 1.0e-12)
        if eta_row > 1.0:
            eta_row = 1.0
        elif eta_row < 0.0:
            eta_row = 0.0
        return eta_row

    def _set_group_outputs(self):
        self.defocus_group_1.v = self._defocus_groups[0]
        self.defocus_group_2.v = self._defocus_groups[1]
        self.defocus_group_3.v = self._defocus_groups[2]
        self.defocus_group_4.v = self._defocus_groups[3]
        self.defocus_group_5.v = self._defocus_groups[4]
        self.defocus_group_6.v = self._defocus_groups[5]
        self.defocus_group_7.v = self._defocus_groups[6]
        self.defocus_group_8.v = self._defocus_groups[7]

        self.temperature_group_1.v = self._temperature_groups[0]
        self.temperature_group_2.v = self._temperature_groups[1]
        self.temperature_group_3.v = self._temperature_groups[2]
        self.temperature_group_4.v = self._temperature_groups[3]
        self.temperature_group_5.v = self._temperature_groups[4]
        self.temperature_group_6.v = self._temperature_groups[5]
        self.temperature_group_7.v = self._temperature_groups[6]
        self.temperature_group_8.v = self._temperature_groups[7]

    def _read_geometry_pairs(self, geom_file_path, n_half_loop):
        diam_in = []
        diam_ret = []
        geom_text = str(geom_file_path).strip() if geom_file_path == geom_file_path else ""
        path = Path(geom_text) if geom_text else Path("__missing_geom__")
        if not geom_text or not path.exists() or path.is_dir():
            path = self._debugging_root() / "header_geom_PID.txt"
        if not path.exists():
            return None, None

        try:
            with path.open("r", encoding="utf-8") as f:
                for line in f:
                    vals = self._extract_numeric_values(line)
                    if len(vals) < 2:
                        continue
                    first = vals[0]
                    second = vals[1]
                    diam_in.append(first)
                    diam_ret.append(second)
                    if len(diam_in) >= n_half_loop:
                        break
        except Exception:
            return None, None

        if len(diam_in) < n_half_loop:
            return None, None

        return np.array(diam_in[:n_half_loop], dtype=float), np.array(diam_ret[:n_half_loop], dtype=float)

    @staticmethod
    def _extract_numeric_values(line):
        matches = re.findall(r"[-+]?\d*\.?\d+(?:[eEdD][-+]?\d+)?", line)
        values = []
        for token in matches:
            values.append(float(token.replace("D", "E").replace("d", "e")))
        return values

    def _next_numeric_line(self, lines, start_idx, min_count=1):
        idx = start_idx
        while idx < len(lines):
            values = self._extract_numeric_values(lines[idx])
            if len(values) >= min_count:
                return values, idx + 1
            idx += 1
        return [0.0] * min_count, idx

    def _load_loop_file_data(self):
        n_loop = self._n_loop_i
        self._ncv_state = np.zeros((4, n_loop), dtype=int)
        self._state_indices = [[np.array([], dtype=int) for _ in range(n_loop)] for _ in range(4)]
        self._h2_pressure = [np.array([], dtype=float) for _ in range(n_loop)]
        self._m_dot_var = np.ones(n_loop, dtype=float)

        sf_label_int = int(round(self._safe(self.sf_label.v, 1.0)))
        loop_file = Path(f"LoopDef_Sector{sf_label_int}") / "loops.txt"
        if not loop_file.exists():
            loop_file = self._debugging_root() / f"LoopDef_Sector{sf_label_int}" / "loops.txt"
        if not loop_file.exists():
            return

        try:
            lines = loop_file.read_text(encoding="utf-8", errors="ignore").splitlines()
        except Exception:
            return

        has_labeled_blocks = any("Loop Number" in line for line in lines)
        if has_labeled_blocks:
            curr_loop = -1
            line_idx = 0
            while line_idx < len(lines) and curr_loop < n_loop - 1:
                line = lines[line_idx]
                if "Loop Number" in line:
                    curr_loop += 1
                    line_idx += 1
                    continue

                def _next_vals(start_idx, min_count=1):
                    j = start_idx
                    while j < len(lines):
                        vals = self._extract_numeric_values(lines[j])
                        if len(vals) >= min_count:
                            return vals, j + 1
                        j += 1
                    return [0.0] * min_count, j

                if curr_loop >= 0 and "Pristine Mapping" in line:
                    vals, line_idx = _next_vals(line_idx + 1, 1)
                    n_state = max(int(round(vals[0])), 0)
                    self._ncv_state[0, curr_loop] = n_state
                    if n_state > 0:
                        idx_vals, line_idx = _next_vals(line_idx, n_state)
                        inds = np.array([int(round(v)) for v in idx_vals[:n_state]], dtype=int)
                        self._state_indices[0][curr_loop] = inds
                    continue

                if curr_loop >= 0 and "Lost Vacuum Mapping" in line:
                    vals, line_idx = _next_vals(line_idx + 1, 1)
                    n_state = max(int(round(vals[0])), 0)
                    self._ncv_state[1, curr_loop] = n_state
                    if n_state > 0:
                        idx_vals, line_idx = _next_vals(line_idx, n_state)
                        inds = np.array([int(round(v)) for v in idx_vals[:n_state]], dtype=int)
                        self._state_indices[1][curr_loop] = inds
                    continue

                if curr_loop >= 0 and "Broken Glass Mapping" in line:
                    vals, line_idx = _next_vals(line_idx + 1, 1)
                    n_state = max(int(round(vals[0])), 0)
                    self._ncv_state[2, curr_loop] = n_state
                    if n_state > 0:
                        idx_vals, line_idx = _next_vals(line_idx, n_state)
                        inds = np.array([int(round(v)) for v in idx_vals[:n_state]], dtype=int)
                        self._state_indices[2][curr_loop] = inds
                    continue

                if curr_loop >= 0 and "H2 Mapping" in line:
                    vals, line_idx = _next_vals(line_idx + 1, 1)
                    n_state = max(int(round(vals[0])), 0)
                    self._ncv_state[3, curr_loop] = n_state
                    if n_state > 0:
                        idx_vals, line_idx = _next_vals(line_idx, n_state)
                        inds = np.array([int(round(v)) for v in idx_vals[:n_state]], dtype=int)
                        self._state_indices[3][curr_loop] = inds
                    continue

                if curr_loop >= 0 and "Hydrogen Annulus Pressure" in line:
                    n_h2_cv = self._ncv_state[3, curr_loop]
                    if n_h2_cv > 0:
                        h2_vals, line_idx = _next_vals(line_idx + 1, n_h2_cv)
                        self._h2_pressure[curr_loop] = np.array(h2_vals[:n_h2_cv], dtype=float)
                    else:
                        line_idx += 1
                    continue

                if curr_loop >= 0 and "Mass flow variability multiplier" in line:
                    m_vals, line_idx = _next_vals(line_idx + 1, 1)
                    self._m_dot_var[curr_loop] = max(float(m_vals[0]), 0.0)
                    continue

                line_idx += 1
            return

        idx = 0
        for loop in range(n_loop):
            # skip Loop Number / Label, then read total nodes and length
            _, idx = self._next_numeric_line(lines, idx, min_count=1)
            _, idx = self._next_numeric_line(lines, idx, min_count=1)
            _, idx = self._next_numeric_line(lines, idx, min_count=2)

            # Pristine / Lost Vacuum / Broken Glass / Hydrogen mappings
            for state in range(4):
                vals, idx = self._next_numeric_line(lines, idx, min_count=1)
                n_state = max(int(round(vals[0])), 0)
                self._ncv_state[state, loop] = n_state
                if n_state > 0:
                    idx_vals, idx = self._next_numeric_line(lines, idx, min_count=n_state)
                    inds = np.array([int(round(v)) for v in idx_vals[:n_state]], dtype=int)
                    self._state_indices[state][loop] = inds

            # H2 pressure values sized by hydrogen state count
            n_h2_cv = self._ncv_state[3, loop]
            if n_h2_cv > 0:
                h2_vals, idx = self._next_numeric_line(lines, idx, min_count=n_h2_cv)
                self._h2_pressure[loop] = np.array(h2_vals[:n_h2_cv], dtype=float)

            # Read in Mass flow variation
            m_vals, idx = self._next_numeric_line(lines, idx, min_count=1)
            self._m_dot_var[loop] = max(float(m_vals[0]), 0.0)

    def _load_minmax2_from_nn_base(self):
        nn_base = str(self.nn_base.v).strip() if self.nn_base.v == self.nn_base.v else ""
        keys = ["pristine_NN", "vacuumLost_NN", "brokenGlass_NN", "H2_NN"]

        loaded_states = 0
        for state_idx, key in enumerate(keys):
            candidates = [Path(f"{key}{nn_base}.txt")]
            debug_root = self._debugging_root()
            if debug_root.exists():
                candidates.append(debug_root / f"{key}{nn_base}.txt")
            if nn_base:
                base_path = Path(nn_base)
                if base_path.is_dir():
                    candidates.insert(0, base_path / f"{key}.txt")
                debug_base = debug_root / nn_base
                if debug_base.is_dir():
                    candidates.insert(0, debug_base / f"{key}.txt")

            file_path = next((candidate for candidate in candidates if candidate.exists()), None)
            if file_path is None:
                continue

            values = []
            try:
                for line in file_path.read_text(encoding="utf-8", errors="ignore").splitlines():
                    line_vals = self._extract_numeric_values(line)
                    if len(line_vals) == 0:
                        continue
                    values.extend(line_vals)
                    if len(values) >= 16:
                        break
            except Exception:
                continue

            if len(values) >= 16:
                self._minmax2[:, state_idx] = np.array(values[:16], dtype=float)
                loaded_states += 1

        return loaded_states > 0

    def _load_nn_models_from_base(self):
        nn_base = str(self.nn_base.v).strip() if self.nn_base.v == self.nn_base.v else ""
        keys = ["pristine_NN", "vacuumLost_NN", "brokenGlass_NN", "H2_NN"]

        self._nn_layers = {}
        self._nn_loaded = False

        loaded_any = False
        for state_idx, key in enumerate(keys, start=1):
            candidates = [Path(f"{key}{nn_base}.txt")]
            debug_root = self._debugging_root()
            if debug_root.exists():
                candidates.append(debug_root / f"{key}{nn_base}.txt")
            if nn_base:
                base_path = Path(nn_base)
                if base_path.is_dir():
                    candidates.insert(0, base_path / f"{key}.txt")
                debug_base = debug_root / nn_base
                if debug_base.is_dir():
                    candidates.insert(0, debug_base / f"{key}.txt")
            else:
                if debug_root.exists():
                    matches = sorted(debug_root.glob(f"{key}*.txt"))
                    if len(matches) > 0:
                        candidates.insert(0, matches[0])

            file_path = next((candidate for candidate in candidates if candidate.exists()), None)
            if file_path is None:
                continue

            try:
                lines = file_path.read_text(encoding="utf-8", errors="ignore").splitlines()
            except Exception:
                continue

            idx = 0
            minmax_vals, idx = self._next_numeric_line(lines, idx, min_count=16)
            self._minmax2[:, state_idx - 1] = np.array(minmax_vals[:16], dtype=float)

            layers = []
            while idx < len(lines):
                dim_vals, idx = self._next_numeric_line(lines, idx, min_count=2)
                row = int(round(dim_vals[0]))
                col = int(round(dim_vals[1]))
                if row <= 0 or col <= 0:
                    break

                w = np.zeros((col, row), dtype=float)
                ok = True
                for nn in range(row):
                    w_vals, idx = self._next_numeric_line(lines, idx, min_count=col)
                    if len(w_vals) < col:
                        ok = False
                        break
                    w[:, nn] = np.array(w_vals[:col], dtype=float)
                if not ok:
                    break

                b_vals, idx = self._next_numeric_line(lines, idx, min_count=row)
                if len(b_vals) < row:
                    break
                b = np.array(b_vals[:row], dtype=float)
                layers.append((w, b))

            if len(layers) > 0:
                self._nn_layers[state_idx] = layers
                loaded_any = True

        self._nn_loaded = loaded_any
        return loaded_any

    @staticmethod
    def _relu(x):
        return np.maximum(x, 0.0)

    def _nn_forward_prop(self, feat_in, state_1b):
        layers = self._nn_layers.get(state_1b)
        if layers is None or len(layers) == 0:
            return np.zeros(feat_in.shape[0], dtype=float)

        curr = feat_in
        for layer_idx, (w, b) in enumerate(layers):
            curr = curr @ w + b
            if layer_idx < len(layers) - 1:
                curr = self._relu(curr)

        heat = curr[:, 0].astype(float)
        state_col = state_1b - 1
        heat = heat - 0.1
        heat = heat * (self._minmax2[14, state_col] - self._minmax2[15, state_col]) + self._minmax2[15, state_col]
        return heat

    def _dt_dtime_nn(self, loop_idx, t_nodes, t_bar, features, m_dot, mc_sf, l_segment, vol):
        n_nodes = len(t_nodes)
        n_cv = n_nodes - 1

        q_in = np.zeros(n_cv, dtype=float)
        for state in range(1, 5):
            inds = self._get_loop_state_indices(loop_idx, state)
            if inds.size == 0:
                continue
            feat_subset = features[inds, :]
            q_in[inds] = self._nn_forward_prop(feat_subset, state) * l_segment

        rho = np.zeros(n_cv, dtype=float)
        c_bar = np.zeros(n_cv, dtype=float)
        for n in range(n_cv):
            rho[n] = self._fluid_density_ees(self.fluid_id.v, t_bar[n], 0.0)
            c_bar[n] = self._fluid_cp(self.fluid_id.v, t_bar[n])

        dt_bar = (q_in + m_dot * c_bar * (t_nodes[:-1] - t_nodes[1:])) / (mc_sf * rho * vol * c_bar)

        dt_nodes = np.zeros(n_nodes, dtype=float)
        dt_nodes[0] = 0.0
        if n_nodes > 2:
            dt_nodes[1:-1] = 0.5 * (dt_bar[:-1] + dt_bar[1:])
        dt_nodes[-1] = dt_bar[-1]
        return dt_nodes

    def _diams_inlet(self, n_cv, n_loop, l_row, l_exp, diam_file, d_default):
        half_loop = max(n_loop // 2, 1)
        diam_half, _ = self._read_geometry_pairs(diam_file, half_loop)
        if diam_half is None:
            return np.full(n_cv, max(d_default, 1.0e-4), dtype=float)

        diam_cv = np.zeros(n_cv, dtype=float)
        cc = 1
        loop_count = 0
        for n in range(n_cv):
            d_curr = diam_half[min(loop_count, len(diam_half) - 1)]
            diam_cv[n] = max(d_curr, 1.0e-4)
            if cc > 1:
                cc += 1
                if cc == 4:
                    cc = 1
                    loop_count += 1
            else:
                cc += 1
                loop_count += 1
        return diam_cv

    def _vols_inlet(self, n_cv, n_loop, l_row, l_exp, diam_file, d_default):
        diam_cv = self._diams_inlet(n_cv, n_loop, l_row, l_exp, diam_file, d_default)
        vols = np.zeros(n_cv, dtype=float)
        cc = 1
        for n in range(n_cv):
            d_curr = diam_cv[n]
            if cc > 1:
                length = (l_row * 2.0 + l_exp * 2.0) / 2.0
                cc += 1
                if cc == 4:
                    cc = 1
            else:
                length = l_row * 2.0
                cc += 1
            vols[n] = math.pi * (d_curr / 2.0) ** 2 * length
        return vols

    def _diams_return(self, n_cv, n_loop, l_row, l_exp, diam_file, d_default):
        half_loop = max(n_loop // 2, 1)
        _, diam_half = self._read_geometry_pairs(diam_file, half_loop)
        if diam_half is None:
            return np.full(n_cv, max(d_default, 1.0e-4), dtype=float)

        diam_cv = np.zeros(n_cv, dtype=float)
        cc = 1
        loop_count = 0
        for n in range(n_cv):
            d_curr = diam_half[min(loop_count, len(diam_half) - 1)]
            diam_cv[n] = max(d_curr, 1.0e-4)
            if cc > 1:
                cc += 1
                if cc == 4:
                    cc = 1
                    loop_count += 1
            else:
                cc += 1
                loop_count += 1
        return diam_cv

    def _vols_return(self, n_cv, n_loop, l_row, l_exp, diam_file, d_default):
        diam_cv = self._diams_return(n_cv, n_loop, l_row, l_exp, diam_file, d_default)
        vols = np.zeros(n_cv, dtype=float)
        cc = 1
        for n in range(n_cv):
            d_curr = diam_cv[n]
            if cc > 1:
                length = (l_row * 2.0 + l_exp * 2.0) / 2.0
                cc += 1
                if cc == 4:
                    cc = 1
            else:
                length = l_row * 2.0
                cc += 1
            vols[n] = math.pi * (d_curr / 2.0) ** 2 * length
        return vols

    def _configure_header_topology(self):
        # Determine control volumes of inlet and return headers
        half_loop = max(self._n_loop_i // 2, 1)
        cc = 0
        for k in range(1, half_loop + 1):
            if k == 1:
                cc += 1
            elif (k % 2) == 0:
                cc += 1
            else:
                cc += 2

        self._n_node_header = max(cc, 2)
        self._n_cv_header = max(self._n_node_header - 1, 1)

        row_distance = max(self._safe(self.row_distance.v, 0.0), 0.0)
        l_exp_loop = max(self._safe(self.l_exp_loop.v, 0.0), 0.0)
        l_base = row_distance * 2.0
        l_exp = (row_distance * 2.0 + l_exp_loop * 2.0) / 2.0

        self._l_cv_inlet = np.zeros(self._n_cv_header, dtype=float)
        self._l_cv_return = np.zeros(self._n_cv_header + 1, dtype=float)
        cc = 1
        for n in range(self._n_cv_header):
            if cc > 1:
                self._l_cv_return[n] = l_exp
                self._l_cv_inlet[n] = l_exp
                cc += 1
                if cc == 4:
                    cc = 1
            else:
                self._l_cv_return[n] = l_base
                self._l_cv_inlet[n] = l_base
                cc += 1
        self._l_cv_return[self._n_cv_header] = l_base

        d_header = max(1.5 * max(self._safe(self.d_receiver.v, 0.07), 1.0e-4), 1.0e-4)
        self._d_inlet = self._diams_inlet(self._n_cv_header, self._n_loop_i, row_distance, l_exp_loop, self.geom_file.v, d_header)
        self._d_return = self._diams_return(self._n_cv_header + 1, self._n_loop_i, row_distance, l_exp_loop, self.geom_file.v, d_header)
        self._vol_inlet = self._vols_inlet(self._n_cv_header, self._n_loop_i, row_distance, l_exp_loop, self.geom_file.v, d_header)
        self._vol_return = self._vols_return(self._n_cv_header + 1, self._n_loop_i, row_distance, l_exp_loop, self.geom_file.v, d_header)

        # Compute indices where mass flow is leaving inlet header
        inds = []
        cc = 1
        for n in range(1, self._n_node_header + 1):
            if cc == 1:
                inds.append(n)
                cc += 1
            elif cc == 2:
                inds.append(n)
                cc += 1
            else:
                cc = 1
        self._inds_header_in = np.array(inds[:half_loop], dtype=int)

        # Determine index mapping (solar field loop to return header)
        inds_left = []
        inds_right = []
        for n in range(self._n_loop_i, 0, -1):
            if (n % 2) == 0:
                inds_right.append(n)
            else:
                inds_left.append(n)
        self._inds_left = np.array(inds_left, dtype=int)
        self._inds_right = np.array(inds_right, dtype=int)

        self._m_dots_in = np.zeros(self._n_cv_header, dtype=float)
        self._m_dots_return = np.zeros(self._n_node_header, dtype=float)

    def _norm_feature(self, value, low_row_1b, high_row_1b, state_1b):
        low = self._minmax2[low_row_1b - 1, state_1b - 1]
        high = self._minmax2[high_row_1b - 1, state_1b - 1]
        denom = high - low
        if abs(denom) <= 1.0e-12:
            return 0.1
        return (value - low) / denom + 0.1

    def _get_loop_state_indices(self, loop_idx, state_1b):
        state_idx = state_1b - 1
        if state_idx < 0 or state_idx >= 4:
            return np.array([], dtype=int)
        if loop_idx < 0 or loop_idx >= self._n_loop_i:
            return np.array([], dtype=int)
        n_state = int(self._ncv_state[state_idx, loop_idx]) if self._ncv_state.shape[1] > loop_idx else 0
        if n_state <= 0:
            return np.array([], dtype=int)
        inds_1b = self._state_indices[state_idx][loop_idx]
        inds_1b = np.array(inds_1b[:n_state], dtype=int)
        return np.clip(inds_1b - 1, 0, self._n_nodes_i - 2)

    def _build_features_for_loop(self, loop_idx, t_bar, dni_array, m_dot_loop):
        n_cv = self._n_nodes_i - 1
        features = np.zeros((n_cv, 7), dtype=float)

        m_dot_eff = m_dot_loop * self._m_dot_var[loop_idx]
        t_amb = self._safe(self.t_amb.v, 300.0)
        wind = self._safe(self.wind.v, 0.0)
        t_sky = self._safe(self.t_sky.v, 300.0)
        w_ap = self._safe(self.w_ap.v, 0.0)

        for state in range(1, 5):
            inds = self._get_loop_state_indices(loop_idx, state)
            if inds.size == 0:
                continue

            features[inds, 0] = self._norm_feature(t_bar[inds], 3, 4, state)
            features[inds, 1] = self._norm_feature(m_dot_eff, 1, 2, state)
            features[inds, 2] = self._norm_feature(t_amb, 7, 8, state)
            features[inds, 3] = self._norm_feature(wind, 11, 12, state)
            features[inds, 4] = self._norm_feature(dni_array[inds] * w_ap, 13, 14, state)

            if state == 4:
                h2_vals = self._h2_pressure[loop_idx]
                for mm, cv_idx in enumerate(inds):
                    h2_p = float(h2_vals[mm]) if mm < len(h2_vals) else 0.0
                    features[cv_idx, 5] = self._norm_feature(h2_p, 5, 6, state)
            else:
                features[inds, 5] = 0.1

            features[inds, 6] = self._norm_feature(t_sky, 9, 10, state)

        return features

    def _refresh_feature_temperature_column(self, loop_idx, features, t_bar_hat):
        for state in range(1, 5):
            inds = self._get_loop_state_indices(loop_idx, state)
            if inds.size == 0:
                continue
            features[inds, 0] = self._norm_feature(t_bar_hat[inds], 3, 4, state)

    def _dt_dtime_nn_surrogate(self, t_nodes, t_bar, features, m_dot, mc_sf, l_segment, vol):
        # TODO(Type4034-NN): Revisit this NN surrogate and replace related
        # forward/derivative subroutines with library-based methods (NumPy/PyTorch/JAX)
        # once NN assets and interfaces are finalized.
        n_nodes = len(t_nodes)
        n_cv = n_nodes - 1

        q_in = np.zeros(n_cv, dtype=float)
        for n in range(n_cv):
            q_wap = (features[n, 4] - 0.1) * (self._minmax2[13, 0] - self._minmax2[12, 0]) + self._minmax2[12, 0]
            q_in[n] = max(q_wap, 0.0)

        rho = np.zeros(n_cv, dtype=float)
        c_bar = np.zeros(n_cv, dtype=float)
        for n in range(n_cv):
            rho[n] = self._fluid_density_ees(self.fluid_id.v, t_bar[n], 0.0)
            c_bar[n] = self._fluid_cp(self.fluid_id.v, t_bar[n])

        dt_bar = (q_in + m_dot * c_bar * (t_nodes[:-1] - t_nodes[1:])) / (mc_sf * rho * vol * c_bar)

        dt_nodes = np.zeros(n_nodes, dtype=float)
        dt_nodes[0] = 0.0
        if n_nodes > 2:
            dt_nodes[1:-1] = 0.5 * (dt_bar[:-1] + dt_bar[1:])
        dt_nodes[-1] = dt_bar[-1]
        return dt_nodes

    def _dt_dt_inlet(self, t_nodes):
        n_nodes = len(t_nodes)
        n_cv = n_nodes - 1
        d_t = np.zeros(n_nodes, dtype=float)
        d_t_bar = np.zeros(n_cv, dtype=float)
        mc_mult = max(self._safe(self.mc_header_mult.v, 1.0), 1.0e-6)
        heat_loss = self._safe(self.inlet_header_heat_loss.v, 0.0)
        p_ref = 0.0

        for n in range(n_cv):
            m_dot = max(self._m_dots_in[n], 0.0)
            vol = max(self._vol_inlet[n], 1.0e-9)
            t_bar = (t_nodes[n] + t_nodes[n + 1]) / 2.0
            c = self._fluid_cp(self.fluid_id.v, t_bar)
            rho = self._fluid_density_ees(self.fluid_id.v, t_bar, p_ref)
            h1 = self._h_dowtherm_a(t_nodes[n])
            h2 = self._h_dowtherm_a(t_nodes[n + 1])
            d_t_bar[n] = (m_dot * (h1 - h2) - heat_loss * self._l_cv_inlet[n]) / max(mc_mult * vol * rho * c, 1.0e-9)

        d_t[0] = 0.0
        for n in range(1, n_nodes - 1):
            d_t[n] = 0.5 * (d_t_bar[n - 1] + d_t_bar[n])
        d_t[n_nodes - 1] = d_t_bar[n_cv - 1]
        return d_t

    def _dt_dt_return(self, t_nodes, t_hold_l, t_hold_r, m_left, m_right):
        n_nodes = len(t_nodes)
        n_cv = n_nodes - 1
        d_t = np.zeros(n_nodes, dtype=float)
        d_t_bar = np.zeros(n_cv, dtype=float)
        mc_mult = max(self._safe(self.mc_header_mult.v, 1.0), 1.0e-6)
        heat_loss = self._safe(self.return_header_heat_loss.v, 0.0)
        p_ref = 0.0

        # Fortran jj = 2 (1-based) -> zero-based index 1
        jj = 1
        for n in range(n_cv):
            rho = self._fluid_density_ees(self.fluid_id.v, (t_nodes[n] + t_nodes[n + 1]) / 2.0, p_ref)
            c = self._fluid_cp(self.fluid_id.v, (t_nodes[n] + t_nodes[n + 1]) / 2.0)
            h1 = self._h_dowtherm_a(t_nodes[n])
            h2 = self._h_dowtherm_a(t_nodes[n + 1])

            if n == 0:
                d_t_bar[n] = (self._m_dots_return[n] * (h1 - h2)) / max(mc_mult * self._vol_return[n] * rho * c, 1.0e-9)
            else:
                has_injection = jj < len(self._inds_header_in) and (n + 1) == int(self._inds_header_in[jj])
                if has_injection and jj < len(m_left) and jj < len(m_right):
                    h_r = self._h_dowtherm_a(t_hold_r[jj])
                    h_l = self._h_dowtherm_a(t_hold_l[jj])
                    d_t_bar[n] = (
                        m_left[jj] * h_l
                        + m_right[jj] * h_r
                        + self._m_dots_return[n - 1] * h1
                        - self._m_dots_return[n] * h2
                        - heat_loss * self._l_cv_return[n]
                    ) / max(mc_mult * self._vol_return[n] * rho * c, 1.0e-9)
                    jj += 1
                else:
                    d_t_bar[n] = (
                        self._m_dots_return[n] * (h1 - h2) - heat_loss * self._l_cv_return[n]
                    ) / max(self._vol_return[n] * rho * c, 1.0e-9)

        d_t[0] = 0.0
        for n in range(1, n_nodes - 1):
            d_t[n] = 0.5 * (d_t_bar[n - 1] + d_t_bar[n])
        d_t[n_nodes - 1] = d_t_bar[n_cv - 1]
        return d_t

    @staticmethod
    def _rk4_step(func, y, dt):
        k1 = func(y)
        k2 = func(y + 0.5 * dt * k1)
        k3 = func(y + 0.5 * dt * k2)
        k4 = func(y + dt * k3)
        return y + (k1 / 6.0 + k2 / 3.0 + k3 / 3.0 + k4 / 6.0) * dt

    def _update_headers_with_rk4(self, t_inlet, m_dot_loop):
        # Update inlet header temperatures and specify inlet loop temperatures
        self._t_header_inlet[0] = t_inlet
        half_loop = max(self._n_loop_i // 2, 1)
        cc = 0
        for n in range(half_loop):
            idx = int(self._inds_header_in[min(n, len(self._inds_header_in) - 1)]) - 1
            idx = max(min(idx, len(self._t_header_inlet) - 1), 0)
            if cc < self._n_loop_i:
                self._loop_t[0, cc] = self._t_header_inlet[idx]
            if cc + 1 < self._n_loop_i:
                self._loop_t[0, cc + 1] = self._t_header_inlet[idx]
            cc += 2

        # Update Return Header Node 1 Temperature
        if self._n_loop_i >= 2:
            num = self._loop_t[-1, self._n_loop_i - 1] * self._m_dot_var[self._n_loop_i - 1] + self._loop_t[-1, self._n_loop_i - 2] * self._m_dot_var[self._n_loop_i - 2]
            den = self._m_dot_var[self._n_loop_i - 1] + self._m_dot_var[self._n_loop_i - 2]
            self._t_header_return[0] = num / max(den, 1.0e-9)

        # Specify mass flow rates out of each return header control volume
        cc = 0
        jj = 1
        for n in range(self._n_cv_header + 1):
            if n == 0:
                ir = int(self._inds_right[min(cc, len(self._inds_right) - 1)]) - 1
                il = int(self._inds_left[min(cc, len(self._inds_left) - 1)]) - 1
                self._m_dots_return[n] = m_dot_loop * (self._m_dot_var[ir] + self._m_dot_var[il])
                cc += 1
            else:
                target = int(self._inds_header_in[min(jj, len(self._inds_header_in) - 1)]) if len(self._inds_header_in) > 0 else -999
                if (n + 1) == target and cc < len(self._inds_right) and cc < len(self._inds_left):
                    ir = int(self._inds_right[cc]) - 1
                    il = int(self._inds_left[cc]) - 1
                    self._m_dots_return[n] = self._m_dots_return[n - 1] + m_dot_loop * (self._m_dot_var[ir] + self._m_dot_var[il])
                    jj += 1
                    cc += 1
                else:
                    self._m_dots_return[n] = self._m_dots_return[n - 1]

        # Specify mass flow rates through each inlet header control volume
        jj = 1
        for n in range(self._n_cv_header):
            target = int(self._inds_header_in[min(jj, len(self._inds_header_in) - 1)]) if len(self._inds_header_in) > 0 else -999
            if (n + 1) == target:
                self._m_dots_in[n] = m_dot_loop * self._n_loop_i * (half_loop - jj) / max(half_loop, 1)
                jj += 1
            else:
                self._m_dots_in[n] = self._m_dots_in[n - 1] if n > 0 else self._m_dots_in[n]

        # Step Return Header Temperatures Through Time (RK4)
        m_left = np.zeros(half_loop, dtype=float)
        m_right = np.zeros(half_loop, dtype=float)
        t_hold_l = np.zeros(half_loop, dtype=float)
        t_hold_r = np.zeros(half_loop, dtype=float)
        for n in range(half_loop):
            il = int(self._inds_left[min(n, len(self._inds_left) - 1)]) - 1
            ir = int(self._inds_right[min(n, len(self._inds_right) - 1)]) - 1
            m_left[n] = m_dot_loop * self._m_dot_var[il]
            m_right[n] = m_dot_loop * self._m_dot_var[ir]
            t_hold_l[n] = self._loop_t[-1, il]
            t_hold_r[n] = self._loop_t[-1, ir]

        dt_seconds = max(self._safe(getattr(self.model.settings, "timestep", 0.0), 0.0) * 3600.0, 0.0)
        self._t_header_return = self._rk4_step(
            lambda y: self._dt_dt_return(y, t_hold_l, t_hold_r, m_left, m_right),
            self._t_header_return,
            dt_seconds,
        )

        # Step Inlet Header Temperatures Through Time (RK4)
        self._t_header_inlet = self._rk4_step(self._dt_dt_inlet, self._t_header_inlet, dt_seconds)
        self._t_header_inlet[0] = t_inlet

    def _initialize_state(self):
        self._n_loop_i = max(int(round(self._safe(self.n_loop.v, 1.0))), 1)
        self._n_nodes_i = max(int(round(self._safe(self.n_nodes_per_loop.v, 2.0))), 2)

        l_tot = max(self._safe(self.l_tot.v, 1.0), 1.0e-6)
        d_receiver = max(self._safe(self.d_receiver.v, 0.07), 1.0e-4)
        l_segment = l_tot / max(self._n_nodes_i - 1, 1)
        self._vol_loop_cv = l_segment * math.pi * (d_receiver / 2.0) ** 2

        t_init_sf = self._safe(self.t_init_sf.v, self._safe(self.temperature.v, 573.15))
        t_init_in = self._safe(self.t_init_in_header.v, t_init_sf)
        t_init_ret = self._safe(self.t_init_return_header.v, t_init_sf)

        self._configure_header_topology()

        self._loop_t = np.full((self._n_nodes_i, self._n_loop_i), t_init_sf, dtype=float)
        self._t_header_inlet = np.full(self._n_node_header, t_init_in, dtype=float)
        self._t_header_return = np.full(self._n_node_header + 1, t_init_ret, dtype=float)
        self._load_loop_file_data()
        self._load_nn_models_from_base()
        self._defocus_mode = np.ones(self._n_loop_i, dtype=float)
        self._time_df = np.zeros(self._n_loop_i, dtype=float)
        self._defocus_groups = np.zeros(8, dtype=float)
        self._temperature_groups = np.zeros(8, dtype=float)
        self._t_bar_sf = np.full(self._n_nodes_i - 1, t_init_sf, dtype=float)
        self._features = np.zeros((self._n_nodes_i - 1, 7), dtype=float)

        fluid = self.fluid_id.v
        rho_init = self._fluid_density(fluid, t_init_sf)
        self._mass_htf_hold = self._n_loop_i * (self._n_nodes_i - 1) * self._vol_loop_cv * rho_init
        self._t4_ave_hold = t_init_sf
        self._is_initialized = True

        # Set the Initial Values of the Outputs (#,Value)
        self.mass_flow_out.v = self._safe(self.mass_flow.v, 0.0)
        self.pressure_out.v = self._safe(self.pressure.v, 0.0)
        self.temperature_out.v = t_init_ret
        self.t4_ave.v = t_init_sf
        self.mass_counter_out.v = self._mass_htf_hold
        self.defocusing.v = 0.0
        self._set_group_outputs()

    def _update_thermal_state_end_of_timestep(self):
        #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        #!! Perform Thermal Computations at the End of Each Timestep !!
        #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        mass_flow = max(self._safe(self.mass_flow.v, 0.0), 0.0)
        ani = max(self._safe(self.ani.v, 0.0), 0.0)
        theta = self._safe(self.theta.v, 0.0)
        phi = self._safe(self.phi.v, 0.0)
        t_inlet = self._safe(self.temperature.v, self._loop_t[0, 0])
        fluid = self.fluid_id.v

        eta_defocus = np.array(
            [
                self._clamp(self._safe(self.eta_defocus_1.v, 1.0)),
                self._clamp(self._safe(self.eta_defocus_2.v, 1.0)),
                self._clamp(self._safe(self.eta_defocus_3.v, 1.0)),
                self._clamp(self._safe(self.eta_defocus_4.v, 1.0)),
            ],
            dtype=float,
        )
        eta_tracking = self._clamp(self._safe(self.eta_tracking.v, 1.0))
        eta_soil = self._clamp(self._safe(self.eta_soil.v, 1.0))
        eta_reflect = self._clamp(self._safe(self.eta_reflect.v, 1.0))
        sf_avail = self._clamp(self._safe(self.sf_avail.v, 1.0))
        eta_tot_base = eta_tracking * eta_soil * eta_reflect

        eta_iam = self._safe(self.iam_a0.v, 1.0) + self._safe(self.iam_a1.v, 0.0) * theta + self._safe(self.iam_a2.v, 0.0) * theta**2
        eta_iam = max(eta_iam, 0.0)

        w_ap = max(self._safe(self.w_ap.v, 0.0), 0.0)
        n_sca = max(int(round(self._safe(self.n_sca.v, 1.0))), 1)
        m_dot_loop = mass_flow / max(self._n_loop_i, 1)

        self._update_headers_with_rk4(t_inlet, m_dot_loop)
        self._defocus_groups[:] = 0.0

        # Define Relevant Tracking Variables
        t_41 = self._safe(self.t_tracking.v, t_inlet) + 2.0
        t_42 = self._safe(self.t_tracking.v, t_inlet) + 4.0
        t_4d = self._safe(self.t_tracking.v, t_inlet) + 6.0

        timestep_s = max(self._safe(getattr(self.model.settings, "timestep", 0.0), 0.0) * 3600.0, 0.0)
        t_bar_sf_accum = np.zeros(self._n_nodes_i - 1, dtype=float)

        #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        #!! Solar Field Loop Calculations              !!
        #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        for loop in range(self._n_loop_i):
            group_index = self._group_index(loop + 1)

            t_sf_hold = self._loop_t[:, loop]
            t_node_out = float(t_sf_hold[-1])
            t_bar = (t_sf_hold[1:] + t_sf_hold[:-1]) / 2.0

            # Defocusing Scheme
            self._defocus_groups[group_index] = (10.0 * self._defocus_groups[group_index]) + 1.0

            mode = self._defocus_mode[loop]
            if mode == 1.0:
                if t_node_out > t_41:
                    mode = 2.0
            elif mode == 2.0:
                if t_node_out > t_42:
                    mode = 3.0
                elif t_node_out < t_41:
                    mode = 1.0
                    self._defocus_groups[group_index] = self._defocus_groups[group_index] - (self._defocus_groups[group_index] % 10.0) + 1.0
                else:
                    self._defocus_groups[group_index] = self._defocus_groups[group_index] - (self._defocus_groups[group_index] % 10.0) + 2.0
            elif mode == 3.0:
                if t_node_out > t_4d:
                    mode = 4.0
                elif t_node_out < t_42:
                    mode = 2.0
                    self._defocus_groups[group_index] = self._defocus_groups[group_index] - (self._defocus_groups[group_index] % 10.0) + 2.0
                else:
                    self._defocus_groups[group_index] = self._defocus_groups[group_index] - (self._defocus_groups[group_index] % 10.0) + 3.0
            elif mode == 4.0:
                if t_node_out < t_4d:
                    mode = 3.0
                    self._defocus_groups[group_index] = self._defocus_groups[group_index] - (self._defocus_groups[group_index] % 10.0) + 3.0
                    self._time_df[loop] = 0.0
                else:
                    self._time_df[loop] = self._time_df[loop] + timestep_s
                    self._defocus_groups[group_index] = self._defocus_groups[group_index] - (self._defocus_groups[group_index] % 10.0) + 4.0

                if self._time_df[loop] > max(self._safe(self.time_lim_df.v, 0.0), 0.0):
                    mode = 5.0
                    self._defocus_groups[group_index] = self._defocus_groups[group_index] - (self._defocus_groups[group_index] % 10.0) + 5.0
            else:
                if t_node_out < t_4d:
                    mode = 3.0
                    self._defocus_groups[group_index] = self._defocus_groups[group_index] - (self._defocus_groups[group_index] % 10.0) + 3.0
                    self._time_df[loop] = 0.0

            self._defocus_mode[loop] = mode

            n_cv = self._n_nodes_i - 1
            dni_array = np.full(n_cv, ani, dtype=float)

            l_sca = max(self._safe(self.l_tot.v, 1.0) / max(n_sca, 1), 1.0e-9)
            end_gain_pred = max(self._safe(self.ave_focal_length.v, 0.0) * math.tan(theta) - self._safe(self.distance_sca.v, 0.0), 0.0) / l_sca
            end_loss_pred = 1.0 - self._safe(self.ave_focal_length.v, 0.0) * math.tan(theta) / l_sca
            if ((loop + 1) % 2) == 0:
                eta_endloss = np.array(
                    [
                        end_loss_pred,
                        end_loss_pred + end_gain_pred,
                        end_loss_pred + end_gain_pred,
                        end_loss_pred,
                    ],
                    dtype=float,
                )
            else:
                eta_endloss = np.array(
                    [
                        end_loss_pred + end_gain_pred,
                        end_loss_pred,
                        end_loss_pred,
                        end_loss_pred + end_gain_pred,
                    ],
                    dtype=float,
                )

            if mode == 5.0:
                dni_array[:] = 0.0
            else:
                n_sca_i = max(n_sca, 1)
                eta_row = self._row_shadow(phi, self._safe(self.row_distance.v, 0.0), w_ap)
                cos_theta = max(abs(math.cos(theta)), 1.0e-12)
                for cc in range(1, n_sca_i + 1):
                    eta_iam = self._safe(self.iam_a0.v, 1.0) + self._safe(self.iam_a1.v, 0.0) * theta / cos_theta + self._safe(self.iam_a2.v, 0.0) * theta**2 / cos_theta
                    eta_iam = max(eta_iam, 0.0)
                    idx0 = int((cc - 1) * n_cv / n_sca_i)
                    idx1 = int(cc * n_cv / n_sca_i)
                    eta_e = eta_endloss[min(cc - 1, len(eta_endloss) - 1)] if cc <= 4 else eta_endloss[-1]
                    eta_df = eta_defocus[int(max(min(mode, 4.0), 1.0)) - 1] if mode > 1.0 and (cc == 3 or cc == 4) else 1.0
                    dni_array[idx0:idx1] = ani * eta_iam * eta_tot_base * eta_row * sf_avail * eta_df * eta_e

            features = self._build_features_for_loop(loop, t_bar, dni_array, m_dot_loop)

            # Step through time with RK-4 (Fortran call sequence)
            l_segment = max(self._safe(self.l_tot.v, 1.0), 1.0e-9) / max(self._n_nodes_i - 1, 1)
            vol = max(self._vol_loop_cv, 1.0e-12)
            mc_sf = max(self._safe(self.mc_receiver_mult.v, 1.0), 1.0e-9)
            m_dot_eff = m_dot_loop * self._m_dot_var[loop]

            if self._nn_loaded:
                k1 = self._dt_dtime_nn(loop, t_sf_hold, t_bar, features, m_dot_eff, mc_sf, l_segment, vol)
            else:
                k1 = self._dt_dtime_nn_surrogate(t_sf_hold, t_bar, features, m_dot_eff, mc_sf, l_segment, vol)
            t_hat = t_sf_hold + k1 * timestep_s / 2.0
            t_bar_hat = (t_hat[1:] + t_hat[:-1]) / 2.0
            self._refresh_feature_temperature_column(loop, features, t_bar_hat)

            if self._nn_loaded:
                k2 = self._dt_dtime_nn(loop, t_hat, t_bar_hat, features, m_dot_eff, mc_sf, l_segment, vol)
            else:
                k2 = self._dt_dtime_nn_surrogate(t_hat, t_bar_hat, features, m_dot_eff, mc_sf, l_segment, vol)
            t_hat = t_sf_hold + k2 * timestep_s / 2.0
            t_bar_hat = (t_hat[1:] + t_hat[:-1]) / 2.0
            self._refresh_feature_temperature_column(loop, features, t_bar_hat)

            if self._nn_loaded:
                k3 = self._dt_dtime_nn(loop, t_hat, t_bar_hat, features, m_dot_eff, mc_sf, l_segment, vol)
            else:
                k3 = self._dt_dtime_nn_surrogate(t_hat, t_bar_hat, features, m_dot_eff, mc_sf, l_segment, vol)
            t_hat = t_sf_hold + k3 * timestep_s
            t_bar_hat = (t_hat[1:] + t_hat[:-1]) / 2.0
            self._refresh_feature_temperature_column(loop, features, t_bar_hat)

            if self._nn_loaded:
                k4 = self._dt_dtime_nn(loop, t_hat, t_bar_hat, features, m_dot_eff, mc_sf, l_segment, vol)
            else:
                k4 = self._dt_dtime_nn_surrogate(t_hat, t_bar_hat, features, m_dot_eff, mc_sf, l_segment, vol)

            self._loop_t[1:, loop] = t_sf_hold[1:] + (
                k1[1:] / 6.0 + k2[1:] / 3.0 + k3[1:] / 3.0 + k4[1:] / 6.0
            ) * timestep_s

            self._features[:, :] = features

            t_bar_sf_accum += (self._loop_t[1:, loop] + self._loop_t[:-1, loop]) / 2.0

        self._t_bar_sf = t_bar_sf_accum / max(self._n_loop_i, 1)

        ind_4 = self._n_nodes_i - int(math.ceil((self._n_nodes_i - 1) / max(n_sca, 1) / 2.0)) - 1
        ind_4 = max(min(ind_4, self._n_nodes_i - 1), 0)
        self._t4_ave_hold = float(np.mean(self._loop_t[ind_4, :]))

        rho_cv = self._fluid_density(fluid, float(np.mean(self._loop_t)))
        self._mass_htf_hold = self._n_loop_i * (self._n_nodes_i - 1) * self._vol_loop_cv * rho_cv

        # TODO(Type4034): Verify header derivative parity numerically against
        # Header_functions equations (dT_dt_inlet and dT_dt_return).
        # TODO(Type4034): Validate Python NN parser against full Fortran NN_data/NN_functions_static
        # for all state files and edge-case formats.
        # NOTE: Loop-file state data are loaded into `_ncv_state`, `_state_indices`, and
        # `_h2_pressure` from `LoopDef_Sector{sf_label}/loops.txt` in Fortran-compatible order.
        # TODO(Type4034): Verify row-shadow/end-loss/per-SCA optics numerically against TRNSYS
        # reference cases from Solar_Position and SF_piping_functions.
        # TODO(Type4034): Verify whether Header_functions should always use H_Dowtherm_A
        # (as in Fortran) or fluid-dependent enthalpy correlations in Python.

    def _update_mass_and_temperature_groups(self):
        # HTF Mass in SF Computation (Only do once because temperatures aren't changing within timestep)
        n_cv_header = max(self._n_cv_header, 1)
        self._mass_htf_hold = 0.0
        self._temperature_groups[:] = 0.0

        # Inlet Header
        vol_header = max(self._vol_loop_cv, 1.0e-9)
        for n in range(n_cv_header):
            t_cv = (self._t_header_inlet[n] + self._t_header_inlet[n + 1]) / 2.0
            self._mass_htf_hold += vol_header * self._fluid_density(self.fluid_id.v, t_cv)

        # Solar Field
        ind_4 = self._n_nodes_i - int(math.ceil((self._n_nodes_i - 1) / max(int(round(self._safe(self.n_sca.v, 1.0))), 1) / 2.0)) - 1
        ind_4 = max(min(ind_4, self._n_nodes_i - 2), 0)

        t4_accum = 0.0
        for loop in range(self._n_loop_i):
            group_index = self._group_index(loop + 1)
            t_ind_4 = float(self._loop_t[ind_4, loop])

            # Temperature levels
            if t_ind_4 <= 370.0:
                d = 1.0
            elif t_ind_4 <= 427.0:
                d = 2.0
            elif t_ind_4 <= 450.0:
                d = 3.0
            elif t_ind_4 <= 483.0:
                d = 4.0
            elif t_ind_4 <= 538.0:
                d = 5.0
            elif t_ind_4 <= 566.0:
                d = 6.0
            elif t_ind_4 <= 594.0:
                d = 7.0
            elif t_ind_4 <= 616.0:
                d = 8.0
            else:
                d = 9.0
            self._temperature_groups[group_index] = 10.0 * self._temperature_groups[group_index] + d

            t_bar = (self._loop_t[1:, loop] + self._loop_t[:-1, loop]) / 2.0
            for n in range(self._n_nodes_i - 1):
                self._mass_htf_hold += self._vol_loop_cv * self._fluid_density(self.fluid_id.v, float(t_bar[n]))

            t4_accum += t_ind_4

        self._t4_ave_hold = t4_accum / max(self._n_loop_i, 1)

        # Return Header
        for n in range(n_cv_header + 1):
            n2 = min(n + 1, len(self._t_header_return) - 1)
            t_cv = (self._t_header_return[n] + self._t_header_return[n2]) / 2.0
            self._mass_htf_hold += vol_header * self._fluid_density(self.fluid_id.v, t_cv)

    def _pressure_drop_surrogate(self):
        #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        #!! Perform all hydraulic calculations here (allow new solution for every iteration) !!
        #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        mass_flow = max(self._safe(self.mass_flow.v, 0.0), 0.0)
        if mass_flow <= 0.0:
            return 0.0

        d_receiver = max(self._safe(self.d_receiver.v, 0.07), 1.0e-4)
        l_tot = max(self._safe(self.l_tot.v, 1.0), 1.0e-6)
        row_distance = max(self._safe(self.row_distance.v, 0.0), 0.0)
        pressure = self._safe(self.pressure.v, 0.0)
        n_loop = max(self._n_loop_i, 1)
        n_cv_header = max(self._n_cv_header, 1)

        t_ref = self._safe(self.temperature_out.v, self._safe(self.temperature.v, self._t4_ave_hold))
        m_dot_htf = mass_flow / n_loop

        # Pressure Drop Across Solar Field Loop
        n_sca = max(int(round(self._safe(self.n_sca.v, 1.0))), 1)
        n_nodes = max(self._n_nodes_i - 1, 1)
        sca_ind_orig = max(n_nodes / n_sca, 1.0)
        sca_ind = sca_ind_orig
        l_segment = l_tot / n_nodes
        dp_loop = 0.0
        for n in range(1, n_nodes + 1):
            if n == 1:
                x1 = 10.0
                x2 = 3.0
            elif n >= sca_ind:
                x1 = 0.0
                x2 = 2.0
                sca_ind += sca_ind_orig
            else:
                x1 = 0.0
                x2 = 0.0
            t_loop_ref = float(self._t_bar_sf[min(n - 1, len(self._t_bar_sf) - 1)]) if len(self._t_bar_sf) > 0 else t_ref
            dp_loop += self._dp_segment(
                m_dot_htf,
                d_receiver,
                l_segment,
                t_loop_ref,
                pressure,
                n_exp=0.0,
                n_con=0.0,
                n_els=x1,
                n_elm=0.0,
                n_ell=0.0,
                n_gav=0.0,
                n_glv=0.0,
                n_chv=0.0,
                n_lw=0.0,
                n_lcv=0.0,
                n_bja=x2,
            )

        # Pressure Drop accounting for inlet/outlet/cross-over piping of loop
        t_iocop_ref = (
            float(self._t_bar_sf[0] + self._t_bar_sf[-1]) / 2.0 if len(self._t_bar_sf) > 1 else t_ref
        )
        dp_iocop = self._dp_segment(
            m_dot_htf,
            d_receiver,
            40.0 + row_distance,
            t_iocop_ref,
            pressure,
            n_exp=0.0,
            n_con=0.0,
            n_els=2.0,
            n_elm=0.0,
            n_ell=0.0,
            n_gav=2.0,
            n_glv=0.0,
            n_chv=0.0,
            n_lw=2.0,
            n_lcv=1.0,
            n_bja=0.0,
        )

        # Pressure Drop Across Inlet Header
        m_dot_header = mass_flow
        dp_inlet_header = 0.0
        cc = 1
        for n in range(n_cv_header):
            x1 = 0.0
            if n > 0:
                if self._d_inlet[n] < self._d_inlet[n - 1]:
                    x1 = 1.0

            if cc > 1:
                if cc == 2:
                    m_dot_header = m_dot_header - 2.0 * m_dot_htf
                cc += 1
                x2 = 2.0
                if cc == 4:
                    cc = 1
            else:
                x2 = 0.0
                cc += 1
                m_dot_header = m_dot_header - 2.0 * m_dot_htf

            t_header_ref = (
                (self._t_header_inlet[min(n, len(self._t_header_inlet) - 1)] + self._t_header_inlet[min(n + 1, len(self._t_header_inlet) - 1)]) / 2.0
            )
            dp_inlet_header += self._dp_segment(
                max(m_dot_header, 0.0),
                max(self._d_inlet[n], 1.0e-4),
                max(self._l_cv_inlet[n], 0.0),
                t_header_ref,
                pressure,
                n_exp=0.0,
                n_con=x1,
                n_els=0.0,
                n_elm=0.0,
                n_ell=x2,
                n_gav=0.0,
                n_glv=0.0,
                n_chv=0.0,
                n_lw=0.0,
                n_lcv=0.0,
                n_bja=0.0,
            )

        # Pressure Drop Across Return Header
        m_dot_header = 0.0
        dp_return_header = 0.0
        cc = 1
        for n in range(n_cv_header):
            x1 = 0.0
            if n > 0:
                if self._d_inlet[n] > self._d_inlet[n - 1]:
                    x1 = 1.0

            if cc > 1:
                if cc == 2:
                    m_dot_header = m_dot_header + 2.0 * m_dot_htf
                cc += 1
                x2 = 2.0
                if cc == 4:
                    cc = 1
            else:
                x2 = 0.0
                cc += 1
                m_dot_header = m_dot_header + 2.0 * m_dot_htf

            t_header_ref = (
                (self._t_header_inlet[min(n, len(self._t_header_inlet) - 1)] + self._t_header_inlet[min(n + 1, len(self._t_header_inlet) - 1)]) / 2.0
            )
            dp_return_header += self._dp_segment(
                min(m_dot_header, mass_flow),
                max(self._d_return[n], 1.0e-4),
                max(self._l_cv_return[n], 0.0),
                t_header_ref,
                pressure,
                n_exp=x1,
                n_con=0.0,
                n_els=0.0,
                n_elm=0.0,
                n_ell=x2,
                n_gav=0.0,
                n_glv=0.0,
                n_chv=0.0,
                n_lw=0.0,
                n_lcv=0.0,
                n_bja=0.0,
            )

        return max(dp_loop + dp_iocop + dp_inlet_header + dp_return_header, 0.0)

        # TODO(Type4034): Confirm fitting-count mapping (elbow vs long-elbow categories)
        # against the original PressureDrop implementation from SF_piping_functions.

    def calculate(self):
        model = getattr(self, "model", None)
        is_first_step = bool(getattr(model, "is_first_step", False))
        is_converged = bool(getattr(model, "is_converged", False))
        timestep_iteration = int(getattr(model, "timestep_iteration", 0)) if model is not None else 0

        if is_first_step or not self._is_initialized:
            self._initialize_state()

        if is_converged:
            self._update_thermal_state_end_of_timestep()

        # HTF Mass in SF Computation (Only do once because temperatures aren't changing within timestep)
        if timestep_iteration == 0:
            self._update_mass_and_temperature_groups()

        p_in = self._safe(self.pressure.v, 0.0)
        dp_tot = self._pressure_drop_surrogate()
        defocusing = float(np.sum(self._defocus_mode) - self._n_loop_i)

        self.mass_flow_out.v = self._safe(self.mass_flow.v, 0.0)
        self.pressure_out.v = p_in - dp_tot
        self.temperature_out.v = float(self._t_header_return[-1])
        self.t4_ave.v = self._t4_ave_hold
        self.mass_counter_out.v = self._mass_htf_hold
        self.defocusing.v = defocusing
        self._set_group_outputs()