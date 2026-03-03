"""
Neural-network helper functions for solar field receiver heat gain.

Converted from solar_field_modules.f90 (NN_data and NN_functions_static modules).

The module provides:
  - ``load_NN``         — load 4 per-state neural networks and their min/max
                          scaling arrays from text files on disk.
  - ``forward_prop``    — 7-ReLU + 1-linear forward propagation through one
                          collector-state network.
  - ``dt_dtime_NN``     — compute nodal HTF temperature time-derivatives for a
                          solar-field loop using NN-predicted heat gain.

All temperatures are in [K], lengths in [m], mass flows in [kg/s].
"""

import os
import numpy as np
from esclab.components.esol_properties import Incompressible as Inc
from esclab.components.flownetwork.sf_piping_helpers import H_Dowtherm_A


# ---------------------------------------------------------------------------
# Module-level NN weight/bias storage (mirrors Fortran module NN_data)
# ---------------------------------------------------------------------------

# w[k] has shape (n_in, n_out, 4) — one slice per collector state (0-based)
# b[k] has shape (n_out, 4)
# minMax has shape (16, 4) — rows are [min_m, max_m, min_T, max_T, …, min_out, max_out]
#   The output scaling sits at rows 14 and 15 (0-based), i.e. minMax[14, state]
#   and minMax[15, state] (matches Fortran minMax(15,state) and minMax(16,state)).

_weights: list[np.ndarray] = []   # list of 8 arrays, each (n_in, n_out, 4)
_biases:  list[np.ndarray] = []   # list of 8 arrays, each (n_out, 4)
_minMax: np.ndarray = np.zeros((16, 4))   # scaling/normalization bounds
_nn_loaded: bool = False
_nn_base: str = ""


# ---------------------------------------------------------------------------
# load_NN — read 4 state-specific NN files and populate module-level arrays
# ---------------------------------------------------------------------------

def load_NN(base: str) -> np.ndarray:
    """Load all four per-state neural networks from text files.

    The Fortran ``load_NN`` subroutine reads four files whose names are
    constructed as ``<state_name><base>.txt`` (e.g. ``pristine_NNv2.txt``).
    Each file has the format::

        <16 minMax values>          ← one line: min/max for all 8 feature pairs
        <row> <col>                 ← dimensions of layer k weight matrix
        <col values>×row            ← weight rows
        <col values>                ← bias row
        ... (repeat for next layer)

    The Fortran code reads weights row-by-row and stores them transposed,
    i.e. ``w1(col, row, 4)`` so that the matrix is (n_in, n_out, 4).  We
    follow the same convention so that ``np.dot(feat, w[:, :, s])`` is correct.

    Parameters
    ----------
    base : str
        Version suffix appended to each state name, e.g. ``"v2"`` produces
        filenames like ``"pristine_NNv2.txt"``.

    Returns
    -------
    np.ndarray, shape (16, 4)
        The ``minMax`` scaling matrix loaded from all four files.
    """
    global _weights, _biases, _minMax, _nn_loaded, _nn_base

    state_names = ["pristine_NN", "vacuumLost_NN", "brokenGlass_NN", "H2_NN"]
    keys = [f"{name}{base}.txt" for name in state_names]

    # Temporary storage per layer per state before we know all shapes
    w_lists: list[list[np.ndarray]] = []   # w_lists[state][layer] = array (n_out, n_in)
    b_lists: list[list[np.ndarray]] = []   # b_lists[state][layer] = array (n_out,)

    for state_idx, fname in enumerate(keys):
        if not os.path.isfile(fname):
            raise FileNotFoundError(
                f"Neural network file not found: '{fname}'. "
                f"Ensure the file exists relative to the current working directory."
            )
        w_state: list[np.ndarray] = []
        b_state: list[np.ndarray] = []

        with open(fname, "r") as fh:
            lines = [ln.strip() for ln in fh if ln.strip()]

        line_idx = 0

        # First line: 16 minMax values for this state
        vals = list(map(float, lines[line_idx].split()))
        _minMax[:, state_idx] = vals[:16]
        line_idx += 1

        # Remaining lines: alternating dimension headers and weight/bias rows
        while line_idx < len(lines):
            dim_parts = lines[line_idx].split()
            if len(dim_parts) < 2:
                line_idx += 1
                continue
            n_row, n_col = int(dim_parts[0]), int(dim_parts[1])
            line_idx += 1

            # Read n_row weight rows (each has n_col values)
            # Fortran stores w(col, row, state) so we build (n_col, n_row) = (n_in, n_out)
            w_layer = np.zeros((n_col, n_row))
            for r in range(n_row):
                row_vals = list(map(float, lines[line_idx].split()))
                w_layer[:, r] = row_vals[:n_col]
                line_idx += 1

            # Read bias row (n_row values)
            bias_vals = list(map(float, lines[line_idx].split()))
            b_layer = np.array(bias_vals[:n_row])
            line_idx += 1

            w_state.append(w_layer)
            b_state.append(b_layer)

        w_lists.append(w_state)
        b_lists.append(b_state)

    # Assemble into per-layer arrays with a trailing state axis
    n_layers = len(w_lists[0])
    _weights = []
    _biases = []
    for k in range(n_layers):
        shapes_w = [w_lists[s][k].shape for s in range(4)]
        n_in_max = max(sh[0] for sh in shapes_w)
        n_out_max = max(sh[1] for sh in shapes_w)
        w_arr = np.zeros((n_in_max, n_out_max, 4))
        b_arr = np.zeros((n_out_max, 4))
        for s in range(4):
            ni, no = shapes_w[s]
            w_arr[:ni, :no, s] = w_lists[s][k]
            b_arr[:no, s] = b_lists[s][k]
        _weights.append(w_arr)
        _biases.append(b_arr)

    _nn_loaded = True
    _nn_base = base
    return _minMax.copy()


# ---------------------------------------------------------------------------
# relu / linear — single forward-propagation layers
# ---------------------------------------------------------------------------

def _relu(curr_feat: np.ndarray, w: np.ndarray, b: np.ndarray) -> np.ndarray:
    """ReLU layer: val = max(curr_feat @ w + b, 0).

    Parameters
    ----------
    curr_feat : np.ndarray, shape (n_curr, n_in)
    w : np.ndarray, shape (n_in, n_out)
    b : np.ndarray, shape (n_out,)

    Returns
    -------
    np.ndarray, shape (n_curr, n_out)
    """
    val = curr_feat @ w + b
    np.maximum(val, 0.0, out=val)
    return val


def _linear(curr_feat: np.ndarray, w: np.ndarray, b: np.ndarray) -> np.ndarray:
    """Linear layer: val = curr_feat @ w + b.

    Parameters
    ----------
    curr_feat : np.ndarray, shape (n_curr, n_in)
    w : np.ndarray, shape (n_in, n_out)
    b : np.ndarray, shape (n_out,)

    Returns
    -------
    np.ndarray, shape (n_curr, n_out)
    """
    return curr_feat @ w + b


# ---------------------------------------------------------------------------
# forward_prop — full 8-layer forward propagation for one state
# ---------------------------------------------------------------------------

def forward_prop(feat_in: np.ndarray, state: int) -> np.ndarray:
    """Forward-propagate features through the neural network for one state.

    Mirrors Fortran ``forward_Prop(feat_in, n_curr, state)``.  The Fortran
    function receives ``feat_in`` as (7, n_curr) (column-major feature matrix)
    and transposes it internally.  We accept the same (7, n_curr) layout.

    Layers 1–7 use ReLU activation; layer 8 uses a linear activation.
    After forward propagation the output is un-normalised::

        heat = (raw_output - 0.1) * (minMax[14,s] - minMax[15,s]) + minMax[15,s]

    which matches Fortran rows 15 and 16 (1-based) → rows 14 and 15 (0-based).

    Parameters
    ----------
    feat_in : np.ndarray, shape (7, n_curr)
        Normalised input feature matrix — one column per collector segment.
    state : int
        Collector state index, 0-based (0=pristine, 1=lostVacuum,
        2=brokenGlass, 3=H2).

    Returns
    -------
    np.ndarray, shape (n_curr,)
        NN heat-gain output in [W/m] (un-normalised).
    """
    if not _nn_loaded:
        raise RuntimeError(
            "Neural network weights not loaded. Call load_NN(base) first."
        )

    s = state
    # Transpose to (n_curr, 7) for matrix multiplication
    feat = feat_in.T  # shape (n_curr, 7)

    # Layers 1–7: ReLU
    n_layers = len(_weights)
    curr = feat
    for k in range(n_layers - 1):
        w = _weights[k][:, :, s]
        b = _biases[k][:, s]
        # Trim to the actual shape for this state
        n_in = curr.shape[1]
        n_out = b.shape[0]
        curr = _relu(curr, w[:n_in, :n_out], b[:n_out])

    # Final layer: linear
    k = n_layers - 1
    w = _weights[k][:, :, s]
    b = _biases[k][:, s]
    n_in = curr.shape[1]
    n_out = b.shape[0]
    curr = _linear(curr, w[:n_in, :n_out], b[:n_out])

    heat = curr[:, 0]

    # Undo normalisation → [W/m]
    # Fortran: heat = (heat - 0.1) * (minMax(15,state) - minMax(16,state)) + minMax(16,state)
    # 0-based: minMax[14, s] and minMax[15, s]
    heat = (heat - 0.1) * (_minMax[14, s] - _minMax[15, s]) + _minMax[15, s]
    return heat


# ---------------------------------------------------------------------------
# dt_dtime_NN — temperature time-derivative for a solar-field loop
# ---------------------------------------------------------------------------

def dt_dtime_NN(
    t: np.ndarray,
    t_bar: np.ndarray,
    features: np.ndarray,
    m_dot: float,
    mc_sf: float,
    L_segment: float,
    Vol: float,
    n_nodes: int,
    nCV_state: np.ndarray,
    inds_pristine: np.ndarray,
    inds_lVacuum: np.ndarray,
    inds_bGlass: np.ndarray,
    inds_H2: np.ndarray,
    fluid: float,
) -> np.ndarray:
    """Compute the HTF nodal temperature time-derivatives for a loop.

    Mirrors Fortran ``dt_dtime_NN`` in ``NN_functions_static``.

    The NN predicts the heat gain [W/m] per control volume for each collector
    state.  HTF properties (density, specific heat) use ``Inc`` (Incompressible)
    in place of the Fortran ``den_NN`` / ``spec_NN`` surrogate functions.

    Parameters
    ----------
    t : np.ndarray, shape (n_nodes,)
        Nodal temperatures [K].
    t_bar : np.ndarray, shape (n_nodes-1,)
        Control-volume average temperatures [K].
    features : np.ndarray, shape (n_nodes-1, 7)
        Normalised NN input features (one row per CV).
    m_dot : float
        Loop mass flow rate [kg/s].
    mc_sf : float
        Receiver thermal-mass multiplier [-].
    L_segment : float
        Length of one control volume [m].
    Vol : float
        Volume of one control volume [m^3].
    n_nodes : int
        Number of nodes in the loop (one more than the number of CVs).
    nCV_state : np.ndarray, shape (4,), dtype int
        Number of CVs in each collector state (pristine, lostVacuum,
        brokenGlass, H2).
    inds_pristine : np.ndarray, shape (n_nodes-1,), dtype int
        0-based CV indices in pristine state.
    inds_lVacuum : np.ndarray, shape (n_nodes-1,), dtype int
        0-based CV indices in lost-vacuum state.
    inds_bGlass : np.ndarray, shape (n_nodes-1,), dtype int
        0-based CV indices in broken-glass state.
    inds_H2 : np.ndarray, shape (n_nodes-1,), dtype int
        0-based CV indices in hydrogen state.
    fluid : float
        HTF fluid identifier passed to Incompressible property functions.
        (Fortran ``dt_dtime_NN`` hard-codes ``fnumd=40`` for Dowtherm A;
        the Python version accepts the fluid ID from the component parameter.)

    Returns
    -------
    np.ndarray, shape (n_nodes,)
        Nodal temperature time-derivatives [K/s].
    """
    n_cv = n_nodes - 1
    state_inds = [inds_pristine, inds_lVacuum, inds_bGlass, inds_H2]

    # ------------------------------------------------------------------
    # Forward-propagate NN for each collector state to get q_in [W] per CV
    # ------------------------------------------------------------------
    q_in = np.zeros(n_cv)
    for state in range(4):
        n = int(nCV_state[state])
        if n > 0:
            curr_inds = state_inds[state][:n]
            # features shape is (n_cv, 7); Fortran passes feat_in as (7, n_curr)
            feat_temp = features[curr_inds, :].T  # shape (7, n)
            # forward_prop expects feat_in (7, n_curr) and returns (n_curr,) in [W/m]
            heat_per_m = forward_prop(feat_temp, state)
            q_in[curr_inds] = heat_per_m * L_segment  # [W]

    # ------------------------------------------------------------------
    # HTF properties at each control volume
    # ------------------------------------------------------------------
    rho = np.array([Inc.density(fluid=fluid, T=t_bar[n], P=0.0) for n in range(n_cv)])
    # AUTO UNITS CONVERSION IMPLEMENTED: spec_NN (Fortran surrogate) returns kJ/kg-K;
    c_bar = np.array([ Inc.specheat(fluid=fluid, T=t_bar[n], P=0.0) for n in range(n_cv)])

    h = np.array([H_Dowtherm_A(t[n]) for n in range(n_nodes)])

    # ------------------------------------------------------------------
    # Compute control-volume temperature time-derivatives
    # Note: Fortran uses the enthalpy-based formulation:
    #   dt_bar = 1/(mc_sf * rho * Vol * c_bar) * (q_in + m_dot * c_bar * (t[n] - t[n+1]))
    # This is equivalent to using H_Dowtherm_A(t[n]) - H_Dowtherm_A(t[n+1]) ≈ c_bar * ΔT
    # for a single-phase fluid; the Fortran code uses c_bar*(t[n]-t[n+1]) directly (see
    # comment in source: "using c_bar*(t[n]-t[n+1])").
    # ------------------------------------------------------------------
    dt_dtheta_bar = (
        1.0 / mc_sf / (rho * Vol * c_bar)
        * (q_in + m_dot * c_bar * (t[:n_cv] - t[1:n_nodes]))
    )

    # ------------------------------------------------------------------
    # Map CV derivatives to nodal derivatives (same stencil as other RK4
    # functions: node 0 = 0, interior = average of adjacent CVs, last = last CV)
    # ------------------------------------------------------------------
    dt_dtheta = np.zeros(n_nodes)
    dt_dtheta[0] = 0.0
    dt_dtheta[1:n_nodes - 1] = 0.5 * (dt_dtheta_bar[:n_cv - 1] + dt_dtheta_bar[1:n_cv])
    dt_dtheta[n_nodes - 1] = dt_dtheta_bar[n_cv - 1]

    return dt_dtheta
