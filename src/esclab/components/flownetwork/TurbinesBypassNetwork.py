"""Type 6028 turbines and bypass network converted from Fortran."""

import math
from functools import lru_cache

import numpy as np
from eeslib import fluid_properties as fp
from scipy.interpolate import RegularGridInterpolator

from esclab.simulate import Component
from esclab.components.esol_properties import Incompressible


def _clip(value, low, high):
    try:
        value = float(value)
    except Exception:
        return low
    if not math.isfinite(value):
        return low
    return max(low, min(high, value))


def _interp1(x0, x1, y0, y1, x):
    x0 = float(x0)
    x1 = float(x1)
    y0 = float(y0)
    y1 = float(y1)
    x = float(x)
    if x1 == x0:
        return y0
    return float(np.interp(x, np.array([x0, x1], dtype=float), np.array([y0, y1], dtype=float)))


@lru_cache(maxsize=16)
def _cached_bilinear_interp(x_axis_key, y_axis_key, table_key):
    x_axis = np.asarray(x_axis_key, dtype=float)
    y_axis = np.asarray(y_axis_key, dtype=float)
    z_table = np.asarray(table_key, dtype=float)
    return RegularGridInterpolator((x_axis, y_axis), z_table, method="linear", bounds_error=False, fill_value=None)


def _bilinear(x, x_axis, y, y_axis, table):
    x_axis_np = np.asarray(x_axis, dtype=float)
    y_axis_np = np.asarray(y_axis, dtype=float)

    x_eval = float(np.clip(float(x), x_axis_np[0], x_axis_np[-1]))
    y_eval = float(np.clip(float(y), y_axis_np[0], y_axis_np[-1]))

    x_key = tuple(float(v) for v in x_axis_np.tolist())
    y_key = tuple(float(v) for v in y_axis_np.tolist())
    table_key = tuple(tuple(float(v) for v in row) for row in table)

    interp = _cached_bilinear_interp(x_key, y_key, table_key)
    return float(interp((x_eval, y_eval)))


def pb_cv_data(valve_type, diameter_m, valve_position):
    """Return C_v from the ESOL PB_CV_data correlation."""
    position = _clip(float(valve_position), 0.0, 1.0)
    if position == 0.0:
        return 1.0e-9

    try:
        valve_type_i = int(round(float(valve_type)))
    except Exception:
        valve_type_i = 1

    try:
        diameter_in = float(diameter_m) * 39.3701
    except Exception:
        diameter_in = 0.0
    if not math.isfinite(diameter_in) or diameter_in < 0.0:
        diameter_in = 0.0

    pos_axis = [0.0, 1.0 / 9.0, 2.0 / 9.0, 3.0 / 9.0, 4.0 / 9.0, 5.0 / 9.0, 6.0 / 9.0, 7.0 / 9.0, 8.0 / 9.0, 1.0]

    if valve_type_i == 1:
        d_axis = [2.0, 3.0, 4.0, 5.0, 6.0, 8.0, 12.0, 16.0, 20.0, 24.0, 28.0, 32.0, 36.0, 40.0]
        cv_table = [
            [0.0, 0.1, 5.0, 12.0, 24.0, 45.0, 64.0, 90.0, 125.0, 135.0],
            [0.0, 0.3, 12.0, 22.0, 39.0, 70.0, 116.0, 183.0, 275.0, 302.0],
            [0.0, 0.5, 17.0, 36.0, 78.0, 139.0, 230.0, 364.0, 546.0, 600.0],
            [0.0, 0.8, 29.0, 61.0, 133.0, 237.0, 392.0, 620.0, 930.0, 1022.0],
            [0.0, 2.0, 45.0, 95.0, 205.0, 366.0, 605.0, 958.0, 1437.0, 1579.0],
            [0.0, 3.0, 89.0, 188.0, 408.0, 727.0, 1202.0, 1903.0, 2854.0, 3136.0],
            [0.0, 5.0, 234.0, 495.0, 1072.0, 1911.0, 3162.0, 5005.0, 7507.0, 8250.0],
            [0.0, 8.0, 464.0, 983.0, 2130.0, 3797.0, 6282.0, 9942.0, 14913.0, 16388.0],
            [0.0, 14.0, 791.0, 1674.0, 3628.0, 6465.0, 10698.0, 16931.0, 25396.0, 27908.0],
            [0.0, 22.0, 1222.0, 2587.0, 5605.0, 9989.0, 16528.0, 26157.0, 39236.0, 43116.0],
            [0.0, 30.0, 1663.0, 3522.0, 7630.0, 12599.0, 20036.0, 30482.0, 46899.0, 58696.0],
            [0.0, 45.0, 2387.0, 4791.0, 8736.0, 13788.0, 20613.0, 31395.0, 48117.0, 68250.0],
            [0.0, 60.0, 3021.0, 6063.0, 11055.0, 17499.0, 26086.0, 39731.0, 60895.0, 86375.0],
            [0.0, 84.0, 4183.0, 8395.0, 15307.0, 24159.0, 36166.0, 55084.0, 84425.0, 119750.0],
        ]
        return max(_bilinear(diameter_in, d_axis, position, pos_axis, cv_table), 1.0e-9)

    if valve_type_i == 2:
        d_axis = [3.0, 4.0, 5.0, 6.0, 8.0, 12.0, 16.0, 20.0, 24.0]
        cv_table = [
            [0.0, 4.7, 16.6, 33.1, 52.1, 78.3, 110.0, 190.0, 210.0, 240.0],
            [0.0, 8.4, 29.3, 58.6, 92.0, 140.0, 200.0, 330.0, 370.0, 420.0],
            [0.0, 13.8, 48.4, 96.5, 160.0, 230.0, 330.0, 550.0, 610.0, 700.0],
            [0.0, 20.9, 73.7, 150.0, 230.0, 350.0, 500.0, 820.0, 930.0, 1050.0],
            [0.0, 38.2, 140.0, 270.0, 420.0, 640.0, 900.0, 1500.0, 1690.0, 1920.0],
            [0.0, 88.4, 310.0, 620.0, 980.0, 1460.0, 2080.0, 3450.0, 3890.0, 4420.0],
            [0.0, 150.0, 530.0, 1060.0, 1660.0, 2490.0, 3540.0, 5870.0, 6620.0, 7520.0],
            [0.0, 270.0, 930.0, 1850.0, 2900.0, 4350.0, 6190.0, 10300.0, 11600.0, 13200.0],
            [0.0, 420.0, 1450.0, 2890.0, 4530.0, 6800.0, 9680.0, 16100.0, 18100.0, 20600.0],
        ]
        return max(_bilinear(diameter_in, d_axis, position, pos_axis, cv_table), 1.0e-9)

    if valve_type_i == 3:
        d_axis = [2.0, 3.0, 4.0, 5.0, 6.0, 8.0, 12.0, 16.0, 20.0, 24.0, 28.0, 32.0, 36.0, 40.0]
        cv_max = [135.0, 302.0, 600.0, 1022.0, 1579.0, 3136.0, 8250.0, 16388.0, 27908.0, 43116.0, 58696.0, 68250.0, 86375.0, 119750.0]
        if diameter_in <= d_axis[0]:
            cv_max_d = cv_max[0] * diameter_in / d_axis[0]
        elif diameter_in >= d_axis[-1]:
            cv_max_d = cv_max[-1] * diameter_in / d_axis[-1]
        else:
            i1 = next(idx for idx in range(1, len(d_axis)) if d_axis[idx] >= diameter_in)
            i0 = i1 - 1
            cv_max_d = _interp1(d_axis[i0], d_axis[i1], cv_max[i0], cv_max[i1], diameter_in)
        return max((cv_max_d / 2.0) * position, 1.0e-9)

    return 1.0e-9


def eta_scc_hpt(m_dot_in, m_dot_rated, n_parallel, pitch_diameter_m, p_in_pa, v_in_m3kg, v_design_in_m3kg, design_ep_pa, n_cv, n_row):
    """HPT efficiency from ESOL6015 eta_SCC_hpt correlation."""
    m_dot_in_eng = m_dot_in * 7936.64
    m_dot_rated_eng = m_dot_rated * 7936.64
    pd_eng = pitch_diameter_m * 39.3701
    p_in_eng = p_in_pa / 6894.76
    design_ep_eng = design_ep_pa / 6894.76
    v_in_eng = v_in_m3kg * 16.0185
    v_design_in_eng = v_design_in_m3kg * 16.0185
    vol_dot_rated_eng = max(m_dot_rated_eng * v_design_in_eng, 1.0e-9)

    if n_row == 1.0:
        tfr = m_dot_in_eng / max(m_dot_rated_eng, 1.0e-9)
        eta_base = 87.0

        delta_eta = (1005200.0 * n_parallel / vol_dot_rated_eng) / 100.0
        eta_base = eta_base - delta_eta * eta_base

        delta_eta = (-0.115 * pd_eng + 4.37) / 100.0
        eta_base = eta_base + delta_eta * eta_base

        x = design_ep_eng / max(p_in_eng, 1.0e-9)
        y = math.log(max(vol_dot_rated_eng, 1.0e-9))
        delta_eta = (11.151 - 63.0 * x - 0.50091 * y + 2.83 * x * y) / 100.0
        eta_base = eta_base + delta_eta * eta_base

        delta_eta = (-21.8085 + 21.8085 * tfr + 0.573908 * pd_eng - 0.573908 * tfr * pd_eng) / 100.0
        eta_base = eta_base + delta_eta * eta_base

        y = math.log(max(p_in_eng / max(design_ep_eng, 1.0e-9), 1.0e-9))
        delta_eta = (
            -60.75
            + 66.85 * tfr
            + 29.75 * tfr ** 2
            - 35.85 * tfr ** 3
            + 17.50 * y
            - 20.02 * y * tfr
            - 0.525 * y * tfr ** 2
            + 3.045 * y * tfr ** 3
        ) / 100.0
        eta_base = eta_base + delta_eta * eta_base

        x = tfr
        y = n_cv
        delta_eta = (-5.4 + 4.395 * x + 0.45 * n_cv - 0.36625 * x * y) / 100.0
        eta_base = eta_base + delta_eta * eta_base
    else:
        tfr = m_dot_in_eng / max(m_dot_rated_eng, 1.0e-9)
        eta_base = 84.0

        delta_eta = (1350000.0 * n_parallel / max(m_dot_in_eng * v_in_eng, 1.0e-9)) / 100.0
        eta_base = eta_base - delta_eta * eta_base

        x = design_ep_eng / max(p_in_eng, 1.0e-9)
        y = math.log(max(vol_dot_rated_eng, 1.0e-9))
        delta_eta = (25.665 - 145.0 * x - 1.33281 * y + 7.53 * x * y) / 100.0
        eta_base = eta_base + delta_eta * eta_base

        x = 1.0 - tfr
        y = p_in_eng / max(design_ep_eng, 1.0e-9)
        delta_eta = (
            42.676909 * x
            - 89.391147 * x ** 2
            + 9.0376638 * x ** 3
            - 26.221836 * x * y
            + 25.549385 * y * x ** 2
            + 8.8283868 * y * x ** 3
            + 4.0479550 * y ** 2 * x
            - 1.4725197 * x ** 2 * y ** 2
            - 4.0183332 * y ** 2 * x ** 3
            - 0.14502211 * y ** 3 * x
            - 0.18580363 * y ** 3 * x ** 2
            + 0.42657518 * x ** 3 * y ** 3
        ) / 100.0
        eta_base = eta_base + delta_eta * eta_base

        x = tfr
        y = n_cv
        delta_eta = (-5.4 + 4.395 * x + 0.45 * n_cv - 0.36625 * x * y) / 100.0
        eta_base = eta_base + delta_eta * eta_base

    return _clip(eta_base / 100.0, 0.05, 1.0)


def eta_scc_lpt(m_dot_in, m_dot_rated, n_parallel, p_in_pa, t_in_k, v_in_m3kg, v_design_in_m3kg, design_ep_pa):
    """LPT efficiency from ESOL6015 eta_SCC_lpt correlation."""
    m_dot_in_eng = m_dot_in * 7936.64
    m_dot_rated_eng = m_dot_rated * 7936.64
    p_in_eng = p_in_pa / 6894.76
    v_in_eng = v_in_m3kg * 16.0185
    v_design_in_eng = v_design_in_m3kg * 16.0185
    vol_dot_rated_eng = max(m_dot_rated_eng * v_design_in_eng, 1.0e-9)

    eta_base = 91.93

    delta_eta = (1270000.0 * n_parallel / vol_dot_rated_eng) / 100.0
    eta_base = eta_base - delta_eta * eta_base

    s_in = float(fp.entropy("water", P=max(p_in_pa, 1.0), T=max(t_in_k, 273.15)))
    h_in = float(fp.enthalpy("water", P=max(p_in_pa, 1.0), T=max(t_in_k, 273.15)))
    s_in_eng = 0.0002388459 * s_in
    h_in_eng = 0.0004299226 * h_in

    fig14 = [
        [28.232252, -92.390491, -625.79590, 207.2301, 70.251642, -22.516388],
        [-0.047796308, 1.2844571, 0.38556961, -0.039652999, -0.27180357, 0.064869467],
        [-0.69791427e-3, -0.17037268e-2, 0.86563845e-3, -0.59510660e-3, 0.39705804e-3, -0.73533255e-4],
        [0.12050837e-5, 0.26826382e-6, -0.67887771e-6, 0.52886157e-6, -0.24106229e-6, 0.37881801e-7],
        [-0.50719109e-9, 0.26393497e-9, 0.38021911e-10, -0.10149993e-9, 0.47757232e-10, -0.70989561e-11],
    ]

    x = math.log10(max(p_in_eng, 1.0e-9))
    if s_in_eng > 2.0041:
        y = min(h_in_eng, 1154.0 + 80.0 * x + 88.0 * x ** 2)
    else:
        y = h_in_eng

    delta_eta = 0.0
    for j in range(5):
        for i in range(6):
            delta_eta += fig14[j][i] * (x ** i) * (y ** j)
    delta_eta /= 100.0
    eta_base = eta_base + delta_eta * eta_base

    return _clip(eta_base / 100.0, 0.05, 1.0)


def h_lpt_stage(h_guess, h_min, h_max, p_stage_pa, delta_s, a0, a1, a2, tol):
    """Solve stage enthalpy using ESOL h_lpt_stage root iteration."""
    h_guess = float(h_guess)
    h_min = float(h_min)
    h_max = float(h_max)
    if h_guess < h_min or h_guess > h_max:
        h_guess = 0.5 * (h_min + h_max)

    error = tol + 10.0
    iterations = 0
    alpha = 0.5
    max_iterations = 50
    h_prev = h_guess
    error_prev = error

    while abs(error) > tol:
        iterations += 1
        s_elep = delta_s + a0 + a1 * h_guess + a2 * (h_guess ** 2)
        s_fit = float(fp.entropy("water", P=max(p_stage_pa, 1.0), h=max(h_guess, 1.0)))
        error = s_elep - s_fit

        if abs(error) < tol:
            break

        if iterations > 1:
            if h_guess != h_prev:
                slope = (error - error_prev) / (h_guess - h_prev)
                intercept = error - slope * h_guess
                if slope != 0.0:
                    h_new = _clip(-intercept / slope, h_min, h_max)
                else:
                    h_new = min(h_guess + 10.0, h_max) if error > 0.0 else max(h_guess - 10.0, h_min)
            else:
                h_new = min(h_guess + 10.0, h_max) if error > 0.0 else max(h_guess - 10.0, h_min)
        else:
            h_new = 0.5 * (h_guess + h_max) if error > 0.0 else 0.5 * (h_guess + h_min)

        error_prev = error
        h_prev = h_guess
        h_guess = h_guess + (h_new - h_guess) * alpha
        if iterations >= max_iterations:
            break

    return _clip(h_guess, h_min, h_max)


def cp_water(pressure_pa, temperature_k):
    """Type6015 f_cp_water(P,T) approximation in J/kg-K using enthalpy differences."""
    p = min(max(float(pressure_pa), 1.0), 2.2e7)
    t = float(temperature_k)
    t_sat = float(fp.temperature("water", P=p, Q=0.0))

    def _h(pt):
        return float(fp.enthalpy("water", P=p, T=max(pt, 273.15)))

    if abs(t - t_sat) > 1.0:
        t_high = t + 1.0
        t_low = t - 1.0
    elif t <= t_sat:
        t_high = t_sat
        t_low = t - 1.0
    else:
        t_high = t + 1.0
        t_low = t_sat

    h_high = _h(t_high)
    h_low = _h(t_low)
    return max((h_high - h_low) / max(t_high - t_low, 1.0e-6), 1.0)


def _u_ph(pressure_pa, enthalpy_jkg):
    p = max(float(pressure_pa), 1.0)
    h = max(float(enthalpy_jkg), 1.0)
    rho = max(float(fp.density("water", P=p, h=h)), 1.0e-9)
    return h - p / rho


def drhodhcp(pressure_pa, enthalpy_jkg, dh):
    """Type6015 drhodhcp(P,h,dh): (drho/dh)_P in kg^2/(m^3·J)."""
    p = min(max(float(pressure_pa), 1.0), 2.0e8)
    h = max(float(enthalpy_jkg), 1.0)
    dh = max(abs(float(dh)), 1.0)
    h_low = max(h - dh, 1.0)
    h_high = h + dh
    rho_low = float(fp.density("water", P=p, h=h_low))
    rho_high = float(fp.density("water", P=p, h=h_high))
    return (rho_high - rho_low) / max(h_high - h_low, 1.0e-9)


def drhodpch(pressure_pa, enthalpy_jkg, dp):
    """Type6015 drhodpch(P,h,dP): (drho/dP)_h in kg/(m^3·Pa)."""
    p = min(max(float(pressure_pa), 1.0), 2.0e8)
    h = max(float(enthalpy_jkg), 1.0)
    dp = max(abs(float(dp)), 1.0)
    p_low = max(p - dp, 1.0)
    p_high = min(p + dp, 2.0e8)
    rho_low = float(fp.density("water", P=p_low, h=h))
    rho_high = float(fp.density("water", P=p_high, h=h))
    return (rho_high - rho_low) / max(p_high - p_low, 1.0e-9)


def dudhcp(pressure_pa, enthalpy_jkg, dh):
    """Type6015 dudhcp(P,h,dh): (du/dh)_P as dimensionless slope."""
    p = min(max(float(pressure_pa), 1.0), 2.0e8)
    h = max(float(enthalpy_jkg), 1.0)
    dh = max(abs(float(dh)), 1.0)
    h_low = max(h - dh, 1.0)
    h_high = h + dh
    u_low = _u_ph(p, h_low)
    u_high = _u_ph(p, h_high)
    return (u_high - u_low) / max(h_high - h_low, 1.0e-9)


def dudpch(pressure_pa, enthalpy_jkg, dp):
    """Type6015 dudpch(P,h,dP): (du/dP)_h in J/(kg·Pa)."""
    p = min(max(float(pressure_pa), 1.0), 2.0e8)
    h = max(float(enthalpy_jkg), 1.0)
    dp = max(abs(float(dp)), 1.0)
    p_low = max(p - dp, 1.0)
    p_high = min(p + dp, 2.0e8)
    u_low = _u_ph(p_low, h)
    u_high = _u_ph(p_high, h)
    return (u_high - u_low) / max(p_high - p_low, 1.0e-9)


def _fric_factor_ic(rough, reynold, guess):
    """Type6015 FricFactor_IC Colebrook iteration."""
    reynold = max(float(reynold), 1.0)
    rough = max(float(rough), 0.0)
    guess = max(float(guess), 1.0e-5)
    if reynold < 2750.0:
        return 64.0 / reynold

    x = 1.0 / math.sqrt(guess)
    test_old = x + 2.0 * math.log10(rough / 3.7 + 2.51 * x / reynold)
    x_old = x
    x *= 0.7
    for _ in range(20):
        test = x + 2.0 * math.log10(rough / 3.7 + 2.51 * x / reynold)
        if abs(test - test_old) <= 1.0e-5:
            return max(1.0 / max(x * x, 1.0e-12), 1.0e-6)
        slope = (test - test_old) / max(x - x_old, 1.0e-12)
        x_old = x
        test_old = test
        if slope != 0.0:
            x = max((slope * x - test) / slope, 1.0e-5)
        else:
            break
    return max(guess, 1.0e-6)


def _viscosity_steam(temperature_k):
    t = max(float(temperature_k), 1.0)
    a = 0.00000000001856
    b = 4209.0
    c = 0.04527
    d = -0.00003376
    return max(a * math.exp(b / t + c * t + d * t ** 2), 1.0e-8)


def convection_dynamicpipe(pressure_pa, enthalpy_jkg, diameter_m, m_dot, ff_guess, t_metal_k):
    """Type6015 convection_dynamicpipe(P,h,D,m_dot,ff_guess,T_metal) in W/m^2-K."""
    p = min(max(float(pressure_pa), 1.0), 2.0e8)
    h = max(float(enthalpy_jkg), 1.0)
    d_pipe = max(float(diameter_m), 1.0e-4)
    m_dot = float(m_dot)
    ff_guess = max(float(ff_guess), 1.0e-5)
    t_metal = float(t_metal_k)

    area = math.pi * d_pipe ** 2 / 4.0

    t_steam = float(fp.temperature("water", P=p, h=h))
    rho_steam = max(float(fp.density("water", P=p, h=h)), 1.0e-9)
    x = min(max(float(fp.quality("water", P=p, h=h)), 0.0), 1.0)
    mu_steam = max(float(fp.viscosity("water", P=p, h=h)), 1.0e-9)

    p_sat = min(max(p, 1.0), 2.2e7)
    t_sat = float(fp.temperature("water", P=p_sat, Q=1.0))
    h_sat_g = float(fp.enthalpy("water", P=p_sat, Q=1.0))
    rho_sat_g = max(float(fp.density("water", P=p_sat, Q=1.0)), 1.0e-9)
    mu_g = max(float(fp.viscosity("water", P=p_sat, Q=1.0)), 1.0e-9)

    h_sat_f = float(fp.enthalpy("water", P=p_sat, Q=0.0))
    rho_sat_f = max(float(fp.density("water", P=p_sat, Q=0.0)), 1.0e-9)
    mu_l = max(float(fp.viscosity("water", P=p_sat, Q=0.0)), 1.0e-9)

    k_steam = max(float(fp.conductivity("water", P=p, h=h)), 1.0e-6)
    k_l = max(float(fp.conductivity("water", P=p_sat, Q=0.0)), 1.0e-6)

    vel = m_dot / max(rho_steam * area, 1.0e-9)
    cp_l = 4200.0
    p_crit = 22064000.0

    if h >= h_sat_g:
        re = rho_steam * vel * d_pipe / max(mu_steam, 1.0e-9)
        if re > 2300.0:
            cp = cp_water(p, t_steam)
            pr = max(mu_steam * cp / k_steam, 1.0e-9)
            ff = _fric_factor_ic(0.0, re, ff_guess)
            nu = ((ff / 8.0) * (re - 1000.0) * pr) / max(1.0 + 12.7 * (ff / 8.0) ** 0.5 * (pr ** (2.0 / 3.0) - 1.0), 1.0e-9)
            return max(nu * k_steam / d_pipe, 0.0)
        return max(3.66 * k_steam / d_pipe, 0.0)

    vel_g = m_dot / max(rho_sat_g * area, 1.0e-9)
    re_g = rho_sat_g * vel_g * d_pipe / max(mu_g, 1.0e-9)
    if re_g >= 35000.0 and x > 0.0:
        z = (1.0 / x - 1.0) ** 0.8 * (p / p_crit) ** 0.4
        g_tot = m_dot / max(area, 1.0e-9)
        j_g = x * g_tot / max(math.sqrt(9.81 * d_pipe * rho_sat_g * max(rho_sat_f - rho_sat_g, 1.0e-9)), 1.0e-9)
        j_g_i = 0.98 * (z + 0.263) ** (-0.62)
        re_l = rho_sat_f * vel_g * d_pipe / max(mu_l, 1.0e-9)
        if j_g >= j_g_i:
            pr_l = mu_l * cp_l / max(k_l, 1.0e-9)
            h_l = 0.23 * re_l ** 0.8 * pr_l ** 0.4 * k_l / d_pipe
            return max(h_l * (1.0 + 3.8 / max(z ** 0.95, 1.0e-9)) * (mu_l / max(14.0 * mu_g, 1.0e-9)) ** (0.0058 + 0.557 * p / p_crit), 0.0)

        j_g_iii = 0.95 * (1.254 + 2.27 * z ** 1.249) ** (-1.0)
        h_nu = 1.32 * re_l ** (1.0 / 3.0) * ((rho_sat_f * max(rho_sat_f - rho_sat_g, 1.0e-9) * 9.81 * k_l ** 3) / max(mu_l ** 2, 1.0e-18)) ** (1.0 / 3.0)
        if j_g <= j_g_iii:
            return max(h_nu, 0.0)
        pr_l = mu_l * cp_l / max(k_l, 1.0e-9)
        h_l = 0.23 * re_l ** 0.8 * pr_l ** 0.4 * k_l / d_pipe
        h_i = h_l * (1.0 + 3.8 / max(z ** 0.95, 1.0e-9)) * (mu_l / max(14.0 * _viscosity_steam(t_steam), 1.0e-9)) ** (0.0058 + 0.557 * p / p_crit)
        return max(h_i + h_nu, 0.0)

    if abs(t_sat - t_metal) > 0.01:
        h_fg = (h_sat_g - h_sat_f) + 3.0 / 8.0 * cp_l * (t_sat - t_metal)
        h_bar = 0.555 * ((9.81 * rho_sat_f * max(rho_sat_f - rho_sat_g, 1.0e-9) * k_l ** 3 * h_fg) / (max(mu_l, 1.0e-9) * abs(t_sat - t_metal) * d_pipe)) ** 0.25
        return max(h_bar, 0.0)
    return 0.0


def stodola_stage(p_in_d, p_out_d, p_in, h_in_d, h_in, m_dot_in_d, m_dot_in):
    """Pressure drop through a turbine stage using Stodola's ellipse relation."""
    p_in_d = max(float(p_in_d), 1.0)
    p_in = max(float(p_in), 1.0)
    h_in_d = max(float(h_in_d), 1.0)
    h_in = max(float(h_in), 1.0)

    v_in_d = max(float(fp.volume("water", P=p_in_d, h=h_in_d)), 1.0e-9)
    v_in = max(float(fp.volume("water", P=p_in, h=h_in)), 1.0e-9)

    pr_d = max(min(float(p_out_d) / p_in_d, 0.999999), 0.001)
    phi_d = max(float(m_dot_in_d) / math.sqrt(p_in_d / v_in_d), 1.0e-12)
    phi_max = math.sqrt(0.999999 / max(1.0 - pr_d ** 2, 1.0e-12)) * phi_d
    phi_a = min(float(m_dot_in) / math.sqrt(p_in / v_in), phi_max)
    term = max(1.0 - (phi_a / phi_d) ** 2 * (1.0 - pr_d ** 2), 0.0)
    pr_a = math.sqrt(term)
    return pr_a * p_in


def valve_massflow(cv, p_in, h_in, p_out):
    """Water valve mass flow approximation converted from Type6015 `valve_massflow`."""
    cv = max(float(cv), 0.0)
    p_in = max(float(p_in), 1.0)
    p_out = max(float(p_out), 1.0)
    h_in = max(float(h_in), 1.0)
    delta_p = p_in - p_out
    if cv <= 0.0 or delta_p <= 0.0:
        return 0.0

    kv = cv / 1.156
    p_sat = min(max(p_in, 1.0), 2.2e7)
    t_in = float(fp.temperature("water", P=p_in, h=h_in))
    rho_in = max(float(fp.density("water", P=p_in, h=h_in)), 1.0e-9)
    h_sat_f = float(fp.enthalpy("water", P=p_sat, Q=0.0))
    h_sat_g = float(fp.enthalpy("water", P=p_sat, Q=1.0))
    t_sat = float(fp.temperature("water", P=p_sat, Q=0.0))

    if h_in <= h_sat_f:
        p_crit = 22064000.0
        p_ref = 87726.1
        t_ref = 373.0
        ln_p1p2 = 8.314 * (1.0 / t_ref - 1.0 / max(t_in, 1.0))
        p_v = p_ref * math.exp(ln_p1p2)
        f_l = 0.62
        f_f = 0.96 - 0.28 * math.sqrt(max(p_v / p_crit, 0.0))
        n1 = 0.1
        rho_ref = 1000.0
        if delta_p < (f_l ** 2 * (p_in - f_f * p_v)):
            q_dot = kv * n1 * math.sqrt(max(delta_p / 1000.0, 0.0) / max(rho_in / rho_ref, 1.0e-12))
            return max(q_dot * rho_in / 3600.0, 0.0)
        q_dot = kv * n1 * f_l * math.sqrt(max((p_in / 1000.0 - f_f * p_v / 1000.0), 0.0) / max(rho_in / rho_ref, 1.0e-12))
        return max(q_dot * rho_in / 3600.0, 0.0)

    if h_in > h_sat_g:
        x_t = 0.35
        n6 = 3.16
        x = delta_p / p_in
        cp = cp_water(p_in, max(t_in, t_sat + 2.0))
        cv_w = max(cp * 0.72, 1.0)
        f_y = cp / cv_w / 1.4
        if x < f_y * x_t:
            y = 1.0 - x / max(3.0 * f_y * x_t, 1.0e-9)
            m_dot = kv * n6 * y * math.sqrt(max(x * p_in / 1000.0 * rho_in, 0.0))
        else:
            y = 0.667
            m_dot = kv * n6 * y * math.sqrt(max(f_y * x_t * p_in / 1000.0 * rho_in, 0.0))
        return max(m_dot / 3600.0, 0.0)

    qual = min(max((h_in - h_sat_f) / max(h_sat_g - h_sat_f, 1.0e-9), 0.0), 1.0)
    rho_sat_g = max(float(fp.density("water", P=p_sat, Q=1.0)), 1.0e-9)

    x_t = 0.35
    n6 = 3.16
    x = delta_p / p_in
    cp = cp_water(p_in, t_sat + 2.0)
    cv_w = max(cp * 0.72, 1.0)
    f_y = cp / cv_w / 1.4
    if x < f_y * x_t:
        y = 1.0 - x / max(3.0 * f_y * x_t, 1.0e-9)
        m_dot_c = kv * n6 * y * math.sqrt(max(x * p_in / 1000.0 * rho_sat_g, 0.0)) / 3600.0
    else:
        y = 0.667
        m_dot_c = kv * n6 * y * math.sqrt(max(f_y * x_t * p_in / 1000.0 * rho_sat_g, 0.0)) / 3600.0

    p_crit = 22064000.0
    p_ref = 87726.1
    t_ref = 373.0
    ln_p1p2 = 8.314 * (1.0 / t_ref - 1.0 / max(t_in, 1.0))
    p_v = p_ref * math.exp(ln_p1p2)
    f_l = 0.62
    f_f = 0.96 - 0.28 * math.sqrt(max(p_v / p_crit, 0.0))
    n1 = 0.1
    rho_liq = 1000.0
    rho_ref = 1000.0
    if delta_p < (f_l ** 2 * (p_in - f_f * p_v)):
        q_dot = kv * n1 * math.sqrt(max(delta_p / 1000.0, 0.0) / max(rho_liq / rho_ref, 1.0e-12))
        m_dot_inc = q_dot * rho_liq / 3600.0
    else:
        q_dot = kv * n1 * f_l * math.sqrt(max((p_in / 1000.0 - f_f * p_v / 1000.0), 0.0) / max(rho_liq / rho_ref, 1.0e-12))
        m_dot_inc = q_dot * rho_liq / 3600.0

    return max(m_dot_c * qual + m_dot_inc * (1.0 - qual), 0.0)


class TurbinesBypassNetwork(Component):
    """
    TRNSYS Type 6028 (ESOL6028-Turbines&BypassNetwork).

    Parameters
    ----------
    p_hpmain_ini ... hx_steamhr_trip
        Type6028 parameters 1..116 exposed as named ``Component.Parameter``
        members.

    Inputs
    ------
    turbine_on_signal ... pid_hptcv
        Type6028 inputs 1..36 exposed as named ``Component.Input`` members.

    Outputs
    -------
    turbine_on_state ... trip_hx_steam_hr
        Type6028 outputs 1..123 exposed as named ``Component.Output`` members.
    """

    p_hpmain_ini = Component.Parameter()
    t_hpmain_ini = Component.Parameter()
    length_hpmain = Component.Parameter()
    d_hpmain = Component.Parameter()
    p_lpmain_ini = Component.Parameter()
    t_lpmain_ini = Component.Parameter()
    length_lpmain = Component.Parameter()
    d_lpmain = Component.Parameter()
    p_aux_ini = Component.Parameter()
    t_aux_ini = Component.Parameter()
    length_auxline = Component.Parameter()
    d_auxline = Component.Parameter()
    mc_hpmain_pipe = Component.Parameter()
    mc_lpmain_pipe = Component.Parameter()
    mc_aux_pipe = Component.Parameter()
    hp_bypass_d = Component.Parameter()
    hp_bypass_vs = Component.Parameter()
    hp_bypass_vt = Component.Parameter()
    hp_aux_d = Component.Parameter()
    hp_aux_vs = Component.Parameter()
    hp_aux_vt = Component.Parameter()
    hp_warmup_d = Component.Parameter()
    hp_warmup_vs = Component.Parameter()
    hp_warmup_vt = Component.Parameter()
    hp_drain_d = Component.Parameter()
    hp_drain_vs = Component.Parameter()
    hp_drain_vt = Component.Parameter()
    aux_da_d = Component.Parameter()
    aux_da_vs = Component.Parameter()
    aux_da_vt = Component.Parameter()
    lp_bypass_d = Component.Parameter()
    lp_bypass_vs = Component.Parameter()
    lp_bypass_vt = Component.Parameter()
    lp_aux_d = Component.Parameter()
    lp_aux_vs = Component.Parameter()
    lp_aux_vt = Component.Parameter()
    lp_warmup_d = Component.Parameter()
    lp_warmup_vs = Component.Parameter()
    lp_warmup_vt = Component.Parameter()
    lp_drain_d = Component.Parameter()
    lp_drain_vs = Component.Parameter()
    lp_drain_vt = Component.Parameter()
    m_dot_turbine_seals = Component.Parameter()
    p_ts_req = Component.Parameter()
    hpt_no_gs = Component.Parameter()
    hpt_parallel_sects = Component.Parameter()
    hpt_cv_number = Component.Parameter()
    hpt_cv_d = Component.Parameter()
    hpt_cv_vs = Component.Parameter()
    hpt_cv_vpd = Component.Parameter()
    gs_diameter = Component.Parameter()
    m_dot_hpt_d = Component.Parameter()
    p_hpt_in_d = Component.Parameter()
    t_hpt_in_d = Component.Parameter()
    p_hpt1_d = Component.Parameter()
    p_hpt_exh_d = Component.Parameter()
    m_dot_hpt1_d = Component.Parameter()
    m_dot_hpt2_d = Component.Parameter()
    hx_ua_d = Component.Parameter()
    m_dot_fw_hx_d = Component.Parameter()
    m_dot_htf_hx_d = Component.Parameter()
    hx_exp = Component.Parameter()
    hx_no_shell = Component.Parameter()
    hx_length = Component.Parameter()
    hx_tube_od = Component.Parameter()
    hx_tube_th = Component.Parameter()
    hx_no_tubes = Component.Parameter()
    fluid_id = Component.Parameter()
    lpt_parallel_sects = Component.Parameter()
    lpt_exp_a0 = Component.Parameter()
    lpt_exp_a1 = Component.Parameter()
    lpt_exp_a2 = Component.Parameter()
    m_dot_lpt_d = Component.Parameter()
    p_lpt_in_d = Component.Parameter()
    t_lpt_in_d = Component.Parameter()
    p_lpt1_d = Component.Parameter()
    p_lpt2_d = Component.Parameter()
    p_lpt3_d = Component.Parameter()
    p_lpt4_d = Component.Parameter()
    m_dot_lpt1_d = Component.Parameter()
    m_dot_lpt2_d = Component.Parameter()
    m_dot_lpt3_d = Component.Parameter()
    m_dot_lpt4_d = Component.Parameter()
    p_cond_d = Component.Parameter()
    hpt_sh_alarm_fl = Component.Parameter()
    hpt_sh_trip_fl = Component.Parameter()
    hpt_sh_alarm_pl = Component.Parameter()
    hpt_sh_trip_pl = Component.Parameter()
    partial_load = Component.Parameter()
    hpt_hightemp_alarm = Component.Parameter()
    hpt_hightemp_timedtrip = Component.Parameter()
    hpt_hightemp_trip = Component.Parameter()
    hpt_exhpres_alarm = Component.Parameter()
    hpt_exhpres_trip = Component.Parameter()
    lpt_sh_alarm = Component.Parameter()
    lpt_sh_trip = Component.Parameter()
    lpt_hightemp_alarm = Component.Parameter()
    lpt_hightemp_timedtrip = Component.Parameter()
    lpt_hightemp_trip = Component.Parameter()
    t_hppipe_ini = Component.Parameter()
    t_lppipe_ini = Component.Parameter()
    t_auxpipe_ini = Component.Parameter()
    hx_hightemp_alarm = Component.Parameter()
    hx_hightemp_trip = Component.Parameter()
    hx_lowtemp_alarm = Component.Parameter()
    hx_lowtemp_trip = Component.Parameter()
    hx_highflow_alarm = Component.Parameter()
    hx_highflow_trip = Component.Parameter()
    hx_highpressure_alarm = Component.Parameter()
    hx_highpressure_trip = Component.Parameter()
    hx_highdt_alarm = Component.Parameter()
    hx_highdt_trip = Component.Parameter()
    hx_htfhr_alarm = Component.Parameter()
    hx_htfhr_trip = Component.Parameter()
    hx_steamhr_alarm = Component.Parameter()
    hx_steamhr_trip = Component.Parameter()

    turbine_on_signal = Component.Input()
    if_auto_hpbypass = Component.Input()
    if_auto_lpbypass = Component.Input()
    if_auto_hpaux = Component.Input()
    if_auto_lpaux = Component.Input()
    if_auto_daaux = Component.Input()
    if_auto_hptcv = Component.Input()
    p_da = Component.Input()
    hpt_cv_vpi = Component.Input()
    hp_bypass_vpi = Component.Input()
    hp_aux_vpi = Component.Input()
    hp_warmup_vpi = Component.Input()
    hp_drain_vpi = Component.Input()
    aux_da_vpi = Component.Input()
    lp_bypass_vpi = Component.Input()
    lp_aux_vpi = Component.Input()
    lp_warmup_vpi = Component.Input()
    lp_drain_vpi = Component.Input()
    p_sgt_in = Component.Input()
    h_sgt_in = Component.Input()
    m_dot_htf = Component.Input()
    htf_p_in = Component.Input()
    htf_t_in = Component.Input()
    p_cond = Component.Input()
    m_dot_hpt1 = Component.Input()
    m_dot_hpt2 = Component.Input()
    m_dot_lpt1 = Component.Input()
    m_dot_lpt2 = Component.Input()
    m_dot_lpt3 = Component.Input()
    m_dot_lpt4 = Component.Input()
    pid_hpbypass = Component.Input()
    pid_lpbypass = Component.Input()
    pid_hpaux = Component.Input()
    pid_lpaux = Component.Input()
    pid_daaux = Component.Input()
    pid_hptcv = Component.Input()

    turbine_on_state = Component.Output()
    hpt_cv_position = Component.Output()
    hp_bypass_position = Component.Output()
    hp_aux_position = Component.Output()
    hp_warmup_position = Component.Output()
    hp_drain_position = Component.Output()
    aux_da_position = Component.Output()
    lp_bypass_position = Component.Output()
    lp_aux_position = Component.Output()
    lp_warmup_position = Component.Output()
    lp_drain_position = Component.Output()
    m_dot_sgt_req = Component.Output()
    m_dot_cond = Component.Output()
    h_cond = Component.Output()
    m_dot_da = Component.Output()
    h_da = Component.Output()
    w_dot_total = Component.Output()
    m_dot_htf_out = Component.Output()
    vol_dot_htf_out = Component.Output()
    htf_p_in_out = Component.Output()
    htf_t_out = Component.Output()
    m_dot_hpt_in = Component.Output()
    p_gs_out = Component.Output()
    h_gs_out = Component.Output()
    m_dot_hpts1 = Component.Output()
    p_hpt1 = Component.Output()
    h_hpt1 = Component.Output()
    m_dot_hpts2 = Component.Output()
    p_hpt2 = Component.Output()
    h_hpt2 = Component.Output()
    m_dot_hpt_exh = Component.Output()
    p_hpt_exh = Component.Output()
    h_hpt_exh = Component.Output()
    t_hpt2 = Component.Output()
    m_dot_ss_drain = Component.Output()
    vol_dot_ss_drain = Component.Output()
    p_drain = Component.Output()
    h_ss_drain = Component.Output()
    t_ss_drain = Component.Output()
    m_dot_ss_steam = Component.Output()
    vol_dot_ss_steam = Component.Output()
    p_steam = Component.Output()
    h_ss_steam = Component.Output()
    t_ss_steam = Component.Output()
    m_dot_steam_out = Component.Output()
    vol_dot_steam_out = Component.Output()
    t_steam_out = Component.Output()
    p_steam_out = Component.Output()
    q_dot_hx = Component.Output()
    eta_od = Component.Output()
    m_dot_lpt_stage1 = Component.Output()
    p_lpt1 = Component.Output()
    t_lpt1 = Component.Output()
    h_lpt1 = Component.Output()
    m_dot_lpt_stage2 = Component.Output()
    p_lpt2 = Component.Output()
    t_lpt2 = Component.Output()
    h_lpt2 = Component.Output()
    m_dot_lpt_stage3 = Component.Output()
    p_lpt3 = Component.Output()
    t_lpt3 = Component.Output()
    h_lpt3 = Component.Output()
    m_dot_lpt_stage4 = Component.Output()
    t_lpt4 = Component.Output()
    p_lpt4 = Component.Output()
    h_lpt4 = Component.Output()
    m_dot_lpt_exh = Component.Output()
    vol_dot_lpt_exh = Component.Output()
    t_lpt_exh = Component.Output()
    h_lpt_exh = Component.Output()
    p_hpmain = Component.Output()
    t_hpmain = Component.Output()
    x_hpmain = Component.Output()
    t_hppipe = Component.Output()
    p_lpmain = Component.Output()
    t_lpmain = Component.Output()
    x_lpmain = Component.Output()
    t_lppipe = Component.Output()
    p_aux = Component.Output()
    t_aux = Component.Output()
    x_aux = Component.Output()
    t_auxpipe = Component.Output()
    m_dot_hp_bypass = Component.Output()
    m_dot_hp_aux = Component.Output()
    m_dot_hp_drain = Component.Output()
    m_dot_hp_warmup = Component.Output()
    m_dot_lp_aux = Component.Output()
    m_dot_lp_bypass = Component.Output()
    m_dot_lp_drain = Component.Output()
    m_dot_lp_warmup = Component.Output()
    m_dot_aux_da = Component.Output()
    friction_factor_hp = Component.Output()
    friction_factor_lp = Component.Output()
    friction_factor_aux = Component.Output()
    trip_turbine_seals = Component.Output()
    alarm_hpt_superheat = Component.Output()
    trip_hpt_superheat = Component.Output()
    alarm_hpt_hightemp = Component.Output()
    timedtrip_hpt_hightemp = Component.Output()
    trip_hpt_hightemp = Component.Output()
    alarm_hpt_exhpres = Component.Output()
    trip_hpt_exhpres = Component.Output()
    alarm_lpt_superheat = Component.Output()
    trip_lpt_superheat = Component.Output()
    alarm_lpt_hightemp = Component.Output()
    timedtrip_lpt_hightemp = Component.Output()
    trip_lpt_hightemp = Component.Output()
    hr_htf_in = Component.Output()
    hr_steam = Component.Output()
    alarm_hx_hightemp_in = Component.Output()
    trip_hx_hightemp_in = Component.Output()
    alarm_hx_lowtemp_out = Component.Output()
    trip_hx_lowtemp_out = Component.Output()
    alarm_hx_highflow = Component.Output()
    trip_hx_highflow = Component.Output()
    alarm_hx_highpressure = Component.Output()
    trip_hx_highpressure = Component.Output()
    alarm_hx_highdt = Component.Output()
    trip_hx_highdt = Component.Output()
    alarm_hx_htf_hr = Component.Output()
    trip_hx_htf_hr = Component.Output()
    alarm_hx_steam_hr = Component.Output()
    trip_hx_steam_hr = Component.Output()

    _all_output_names = (
        "turbine_on_state", "hpt_cv_position", "hp_bypass_position", "hp_aux_position", "hp_warmup_position",
        "hp_drain_position", "aux_da_position", "lp_bypass_position", "lp_aux_position", "lp_warmup_position",
        "lp_drain_position", "m_dot_sgt_req", "m_dot_cond", "h_cond", "m_dot_da", "h_da", "w_dot_total",
        "m_dot_htf_out", "vol_dot_htf_out", "htf_p_in_out", "htf_t_out", "m_dot_hpt_in", "p_gs_out", "h_gs_out",
        "m_dot_hpts1", "p_hpt1", "h_hpt1", "m_dot_hpts2", "p_hpt2", "h_hpt2", "m_dot_hpt_exh", "p_hpt_exh",
        "h_hpt_exh", "t_hpt2", "m_dot_ss_drain", "vol_dot_ss_drain", "p_drain", "h_ss_drain", "t_ss_drain",
        "m_dot_ss_steam", "vol_dot_ss_steam", "p_steam", "h_ss_steam", "t_ss_steam", "m_dot_steam_out",
        "vol_dot_steam_out", "t_steam_out", "p_steam_out", "q_dot_hx", "eta_od", "m_dot_lpt_stage1", "p_lpt1", "t_lpt1",
        "h_lpt1", "m_dot_lpt_stage2", "p_lpt2", "t_lpt2", "h_lpt2", "m_dot_lpt_stage3", "p_lpt3", "t_lpt3",
        "h_lpt3", "m_dot_lpt_stage4", "t_lpt4", "p_lpt4", "h_lpt4", "m_dot_lpt_exh", "vol_dot_lpt_exh",
        "t_lpt_exh", "h_lpt_exh", "p_hpmain", "t_hpmain", "x_hpmain", "t_hppipe", "p_lpmain", "t_lpmain",
        "x_lpmain", "t_lppipe", "p_aux", "t_aux", "x_aux", "t_auxpipe", "m_dot_hp_bypass", "m_dot_hp_aux",
        "m_dot_hp_drain", "m_dot_hp_warmup", "m_dot_lp_aux", "m_dot_lp_bypass", "m_dot_lp_drain", "m_dot_lp_warmup",
        "m_dot_aux_da", "friction_factor_hp", "friction_factor_lp", "friction_factor_aux", "trip_turbine_seals",
        "alarm_hpt_superheat", "trip_hpt_superheat", "alarm_hpt_hightemp", "timedtrip_hpt_hightemp", "trip_hpt_hightemp",
        "alarm_hpt_exhpres", "trip_hpt_exhpres", "alarm_lpt_superheat", "trip_lpt_superheat", "alarm_lpt_hightemp",
        "timedtrip_lpt_hightemp", "trip_lpt_hightemp", "hr_htf_in", "hr_steam", "alarm_hx_hightemp_in", "trip_hx_hightemp_in",
        "alarm_hx_lowtemp_out", "trip_hx_lowtemp_out", "alarm_hx_highflow", "trip_hx_highflow", "alarm_hx_highpressure",
        "trip_hx_highpressure", "alarm_hx_highdt", "trip_hx_highdt", "alarm_hx_htf_hr", "trip_hx_htf_hr",
        "alarm_hx_steam_hr", "trip_hx_steam_hr",
    )
    _inc_props = Incompressible()

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
        for output_name in component._all_output_names:
            getattr(component, output_name).v = 0.0

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

    @staticmethod
    def _vp_new(vp_current, vp_request, valve_speed_dps, ts_sec):
        vp_current = max(min(vp_current, 1.0), 0.0)
        vp_request = max(min(vp_request, 1.0), 0.0)
        valve_speed_dps = min(max(valve_speed_dps, 1.0), 90.0)

        vp_current_deg = vp_current * 90.0
        vp_request_deg = vp_request * 90.0
        if vp_current_deg == vp_request_deg:
            return vp_current
        if vp_request_deg > vp_current_deg:
            return min(vp_current_deg + valve_speed_dps * ts_sec, vp_request_deg) / 90.0
        return max(vp_current_deg - valve_speed_dps * ts_sec, vp_request_deg) / 90.0

    def calculate(self):
        ts_sec = max(self.model.settings.timestep * 3600.0, 1.0e-9)
        n_int = max(int(math.ceil(60.0 / ts_sec)), 1)
        is_first_iteration = bool(getattr(self.model, "is_first_iteration", True))
        is_converged = bool(getattr(self.model, "is_converged", False))

        turbine_on = 1.0 if self._safe(self.turbine_on_signal.v, 0.0) == 1.0 else 0.0
        p_sgt_in = max(self._safe(self.p_sgt_in.v, 8.0e6), 1.0)
        h_sgt_in = self._safe(self.h_sgt_in.v, 2.8e6)
        m_dot_sgt_in = max(self._safe(self.m_dot_sgt_req.v, self._safe(self.m_dot_hpt_d.v, 0.0)), 0.0)

        # Previous valve positions are used for valve-speed-limited position updates.
        hpt_cv_prev = max(min(self._safe(self.hpt_cv_position.v, 0.0), 1.0), 0.0)
        hp_bypass_prev = max(min(self._safe(self.hp_bypass_position.v, 0.0), 1.0), 0.0)
        hp_aux_prev = max(min(self._safe(self.hp_aux_position.v, 0.0), 1.0), 0.0)
        hp_warmup_prev = max(min(self._safe(self.hp_warmup_position.v, 0.0), 1.0), 0.0)
        hp_drain_prev = max(min(self._safe(self.hp_drain_position.v, 0.0), 1.0), 0.0)
        aux_da_prev = max(min(self._safe(self.aux_da_position.v, 0.0), 1.0), 0.0)
        lp_bypass_prev = max(min(self._safe(self.lp_bypass_position.v, 0.0), 1.0), 0.0)
        lp_aux_prev = max(min(self._safe(self.lp_aux_position.v, 0.0), 1.0), 0.0)
        lp_warmup_prev = max(min(self._safe(self.lp_warmup_position.v, 0.0), 1.0), 0.0)
        lp_drain_prev = max(min(self._safe(self.lp_drain_position.v, 0.0), 1.0), 0.0)

        # AUTO/MANUAL valve command selection from control inputs.
        hpt_cv_req = self._safe(self.pid_hptcv.v, self._safe(self.hpt_cv_vpi.v, 1.0)) if self._safe(self.if_auto_hptcv.v, 0.0) == 1.0 else self._safe(self.hpt_cv_vpi.v, 1.0)
        hp_bypass_req = self._safe(self.pid_hpbypass.v, self._safe(self.hp_bypass_vpi.v, 0.0)) if self._safe(self.if_auto_hpbypass.v, 0.0) == 1.0 else self._safe(self.hp_bypass_vpi.v, 0.0)
        hp_aux_req = self._safe(self.pid_hpaux.v, self._safe(self.hp_aux_vpi.v, 0.0)) if self._safe(self.if_auto_hpaux.v, 0.0) == 1.0 else self._safe(self.hp_aux_vpi.v, 0.0)
        aux_da_req = self._safe(self.pid_daaux.v, self._safe(self.aux_da_vpi.v, 0.0)) if self._safe(self.if_auto_daaux.v, 0.0) == 1.0 else self._safe(self.aux_da_vpi.v, 0.0)
        lp_bypass_req = self._safe(self.pid_lpbypass.v, self._safe(self.lp_bypass_vpi.v, 0.0)) if self._safe(self.if_auto_lpbypass.v, 0.0) == 1.0 else self._safe(self.lp_bypass_vpi.v, 0.0)
        lp_aux_req = self._safe(self.pid_lpaux.v, self._safe(self.lp_aux_vpi.v, 0.0)) if self._safe(self.if_auto_lpaux.v, 0.0) == 1.0 else self._safe(self.lp_aux_vpi.v, 0.0)
        hp_warmup_req = self._safe(self.hp_warmup_vpi.v, 0.0)
        hp_drain_req = self._safe(self.hp_drain_vpi.v, 0.0)
        lp_warmup_req = self._safe(self.lp_warmup_vpi.v, 0.0)
        lp_drain_req = self._safe(self.lp_drain_vpi.v, 0.0)

        # Valve speed model from Type6015 VP_new (deg/s and timestep-limited movement).
        # Apply valve movement only on first iteration of a timestep.
        if is_first_iteration:
            hpt_cv_vpi = self._vp_new(hpt_cv_prev, hpt_cv_req, self._safe(self.hpt_cv_vs.v, 90.0), ts_sec)
            hp_bypass_vpi = self._vp_new(hp_bypass_prev, hp_bypass_req, self._safe(self.hp_bypass_vs.v, 90.0), ts_sec)
            hp_aux_vpi = self._vp_new(hp_aux_prev, hp_aux_req, self._safe(self.hp_aux_vs.v, 90.0), ts_sec)
            hp_warmup_vpi = self._vp_new(hp_warmup_prev, hp_warmup_req, self._safe(self.hp_warmup_vs.v, 90.0), ts_sec)
            hp_drain_vpi = self._vp_new(hp_drain_prev, hp_drain_req, self._safe(self.hp_drain_vs.v, 90.0), ts_sec)
            aux_da_vpi = self._vp_new(aux_da_prev, aux_da_req, self._safe(self.aux_da_vs.v, 90.0), ts_sec)
            lp_bypass_vpi = self._vp_new(lp_bypass_prev, lp_bypass_req, self._safe(self.lp_bypass_vs.v, 90.0), ts_sec)
            lp_aux_vpi = self._vp_new(lp_aux_prev, lp_aux_req, self._safe(self.lp_aux_vs.v, 90.0), ts_sec)
            lp_warmup_vpi = self._vp_new(lp_warmup_prev, lp_warmup_req, self._safe(self.lp_warmup_vs.v, 90.0), ts_sec)
            lp_drain_vpi = self._vp_new(lp_drain_prev, lp_drain_req, self._safe(self.lp_drain_vs.v, 90.0), ts_sec)
        else:
            hpt_cv_vpi = hpt_cv_prev
            hp_bypass_vpi = hp_bypass_prev
            hp_aux_vpi = hp_aux_prev
            hp_warmup_vpi = hp_warmup_prev
            hp_drain_vpi = hp_drain_prev
            aux_da_vpi = aux_da_prev
            lp_bypass_vpi = lp_bypass_prev
            lp_aux_vpi = lp_aux_prev
            lp_warmup_vpi = lp_warmup_prev
            lp_drain_vpi = lp_drain_prev

        # Reheater and condenser side inputs
        m_dot_htf = self._safe(self.m_dot_htf.v, 0.0)
        p_htf_in = self._safe(self.htf_p_in.v, 2.0e6)
        t_htf_in = self._safe(self.htf_t_in.v, 560.0)
        fluid_name = str(self.fluid_id.v)
        p_cond = self._safe(self.p_cond.v, 1.0e5)

        self._set_all_outputs_zero(self)

        # Outputs 1-11 are control/status passthroughs
        self.turbine_on_state.v = turbine_on
        self.hpt_cv_position.v = hpt_cv_vpi
        self.hp_bypass_position.v = hp_bypass_vpi
        self.hp_aux_position.v = hp_aux_vpi
        self.hp_warmup_position.v = hp_warmup_vpi
        self.hp_drain_position.v = hp_drain_vpi
        self.aux_da_position.v = aux_da_vpi
        self.lp_bypass_position.v = lp_bypass_vpi
        self.lp_aux_position.v = lp_aux_vpi
        self.lp_warmup_position.v = lp_warmup_vpi
        self.lp_drain_position.v = lp_drain_vpi

        # Previous-timestep pipe states used by turbine, bypass, and separator blocks
        p_hpmain_prev = self._safe(self.p_hpmain.v, self._safe(self.p_hpmain_ini.v, p_sgt_in))
        if p_hpmain_prev <= 0.0:
            p_hpmain_prev = self._safe(self.p_hpmain_ini.v, p_sgt_in)
        p_hpmain_prev = max(p_hpmain_prev, 700.0)
        t_hpmain_prev = self._safe(self.t_hpmain.v, self._safe(self.t_hpmain_ini.v, 520.0))
        if t_hpmain_prev <= 0.0:
            t_hpmain_prev = self._safe(self.t_hpmain_ini.v, 520.0)
        p_lpmain_prev = self._safe(self.p_lpmain.v, self._safe(self.p_lpmain_ini.v, p_cond))
        if p_lpmain_prev <= 0.0:
            p_lpmain_prev = self._safe(self.p_lpmain_ini.v, p_cond)
        p_lpmain_prev = max(p_lpmain_prev, 700.0)
        t_lpmain_prev = self._safe(self.t_lpmain.v, self._safe(self.t_lpmain_ini.v, 380.0))
        if t_lpmain_prev <= 0.0:
            t_lpmain_prev = self._safe(self.t_lpmain_ini.v, 380.0)
        p_aux_prev = self._safe(self.p_aux.v, self._safe(self.p_aux_ini.v, p_cond))
        if p_aux_prev <= 0.0:
            p_aux_prev = self._safe(self.p_aux_ini.v, p_cond)
        p_aux_prev = max(p_aux_prev, 700.0)
        t_aux_prev = self._safe(self.t_aux.v, self._safe(self.t_aux_ini.v, 360.0))
        if t_aux_prev <= 0.0:
            t_aux_prev = self._safe(self.t_aux_ini.v, 360.0)
        t_hppipe_prev = self._safe(self.t_hppipe.v, self._safe(self.t_hppipe_ini.v, t_hpmain_prev))
        t_lppipe_prev = self._safe(self.t_lppipe.v, self._safe(self.t_lppipe_ini.v, t_lpmain_prev))
        t_auxpipe_prev = self._safe(self.t_auxpipe.v, self._safe(self.t_auxpipe_ini.v, t_aux_prev))

        h_hpmain_prev = float(fp.enthalpy("water", P=max(p_hpmain_prev, 1.0), T=max(t_hpmain_prev, 273.15)))
        h_lpmain_prev = float(fp.enthalpy("water", P=max(p_lpmain_prev, 1.0), T=max(t_lpmain_prev, 273.15)))
        h_aux_prev = float(fp.enthalpy("water", P=max(p_aux_prev, 1.0), T=max(t_aux_prev, 273.15)))

        # HIGH PRESSURE TURBINE MASS FLOW
        # Turbine block using ESOL6015 helper correlations
        if turbine_on == 1.0:
            m_dot_hpt_d = max(self._safe(self.m_dot_hpt_d.v, 1.0), 1.0e-6)
            m_dot_hpt1_extr = max(self._safe(self.m_dot_hpt1.v, 0.0), 0.0)
            m_dot_hpt2_extr = max(self._safe(self.m_dot_hpt2.v, 0.0), 0.0)

            cv_design = pb_cv_data(1, self.hpt_cv_d.v, max(min(self.hpt_cv_vpd.v, 1.0), 1.0e-6))
            cv_current = pb_cv_data(1, self.hpt_cv_d.v, hpt_cv_vpi)
            hpt_flow_fraction = max(cv_current / max(cv_design, 1.0e-9), 0.0)
            hpt_cv_limit = max(m_dot_sgt_in * hpt_flow_fraction, 1.0e-6)
            m_dot_hpt_in = min(max(self._safe(self.m_dot_hpt_in.v, m_dot_hpt_d), 1.0e-6), hpt_cv_limit)
            m_dot_hpt_max = min(max(1.25 * m_dot_hpt_d, 1.0e-6), max(m_dot_sgt_in, 1.0e-6))

            p_hpt1 = max(self._safe(self.p_hpt1.v, 0.65 * p_hpmain_prev), p_cond)
            p_hpt2 = max(self._safe(self.p_hpt2.v, 0.35 * p_hpmain_prev), p_cond)
            p_hpt_exh = p_hpt2
            h_hpt1 = self._safe(self.h_hpt1.v, h_hpmain_prev)
            h_hpt2 = self._safe(self.h_hpt2.v, h_hpmain_prev)

            # Iterate HPT mass flow and stage pressure until HPT exhaust pressure matches LP main pressure.
            for _ in range(35):
                p_hpt_in = p_hpmain_prev
                h_hpt_in = h_hpmain_prev

                rho_hpt_in = max(float(fp.density("water", P=max(p_hpt_in, 1.0), h=max(h_hpt_in, 1.0))), 1.0e-6)
                rho_hpt_design = max(float(fp.density("water", P=max(self._safe(self.p_hpt_in_d.v, p_hpt_in), 1.0), T=max(self._safe(self.t_hpt_in_d.v, 773.15), 273.15))), 1.0e-6)
                v_hpt_in = 1.0 / rho_hpt_in
                v_hpt_design = 1.0 / rho_hpt_design

                # Find pressure loss from passing through the governing stage based on control valves.
                cv_control = pb_cv_data(1, self.hpt_cv_d.v, hpt_cv_vpi)
                n_cv = max(self._safe(self.hpt_cv_number.v, 1.0), 1.0)
                m_dot_hpt_cv = m_dot_hpt_in / n_cv
                hpt_vol_in_gpm = m_dot_hpt_cv * v_hpt_in * 15850.3
                if cv_control > 1.0e-9 and hpt_vol_in_gpm > 1.0e-9:
                    spec_grav = 0.001 / max(v_hpt_in, 1.0e-9)
                    delta_p_psi = spec_grav / max((cv_control / hpt_vol_in_gpm) ** 2, 1.0e-12)
                    delta_p_cv = delta_p_psi * 6894.76
                else:
                    delta_p_cv = 0.0

                p_gs_out = max(p_hpt_in - delta_p_cv, p_cond + 1000.0)
                h_gs_out = h_hpt_in

                eta_hpt = eta_scc_hpt(
                    m_dot_hpt_in,
                    m_dot_hpt_d,
                    max(self._safe(self.hpt_parallel_sects.v, 1.0), 1.0),
                    max(self._safe(self.gs_diameter.v, 0.5), 0.01),
                    p_hpt_in,
                    v_hpt_in,
                    v_hpt_design,
                    max(self._safe(self.p_hpt_exh_d.v, p_hpt_exh), p_cond),
                    n_cv,
                    1.0 if self._safe(self.hpt_no_gs.v, 1.0) < 1.5 else 2.0,
                )

                m_dot_hpts2_d = max(m_dot_hpt_d - max(self._safe(self.m_dot_hpt1_d.v, 0.0), 0.0), 1.0e-6)
                m_dot_hpt_s2 = max(m_dot_hpt_in - m_dot_hpt1_extr, 1.0e-6)

                p_hpt1 = stodola_stage(
                    max(self._safe(self.p_hpt_in_d.v, p_hpt_in), 1.0),
                    max(self._safe(self.p_hpt1_d.v, p_hpt_in * 0.7), p_cond),
                    p_gs_out,
                    max(self._safe(self.h_gs_out.v, h_gs_out), 1.0),
                    h_gs_out,
                    m_dot_hpt_d,
                    m_dot_hpt_in,
                )
                p_hpt1 = max(p_hpt1, p_cond + 1000.0)

                p_hpt2 = stodola_stage(
                    max(self._safe(self.p_hpt1_d.v, p_hpt1), p_cond + 1000.0),
                    max(self._safe(self.p_hpt_exh_d.v, p_hpt1 * 0.6), p_cond),
                    p_hpt1,
                    max(self._safe(self.h_hpt1.v, h_hpt1), 1.0),
                    h_hpt1,
                    m_dot_hpts2_d,
                    m_dot_hpt_s2,
                )
                p_hpt2 = max(p_hpt2, p_cond + 1000.0)
                p_hpt_exh = p_hpt2

                s_gs_out = float(fp.entropy("water", P=max(p_gs_out, 1.0), h=max(h_gs_out, 1.0)))
                h_hpt2_s = float(fp.enthalpy("water", P=max(p_hpt2, 1.0), s=s_gs_out))
                h_hpt2 = h_gs_out - eta_hpt * max(h_gs_out - h_hpt2_s, 0.0)
                h_hpt1_s = float(fp.enthalpy("water", P=max(p_hpt1, 1.0), s=s_gs_out))
                h_hpt1 = h_gs_out - eta_hpt * max(h_gs_out - h_hpt1_s, 0.0)

                pressure_error_hpt = p_hpt2 - p_lpmain_prev
                if abs(pressure_error_hpt) <= 1000.0 and pressure_error_hpt >= 0.0:
                    break
                step = 0.02 * m_dot_hpt_d * pressure_error_hpt / max(p_hpt_in, 1.0)
                m_dot_hpt_in = max(min(m_dot_hpt_in + step, m_dot_hpt_max), 1.0e-6)
                m_dot_hpt_in = min(m_dot_hpt_in, hpt_cv_limit)

            m_dot_hpt_s1 = m_dot_hpt_in
            m_dot_hpt_s2 = max(m_dot_hpt_s1 - m_dot_hpt1_extr, 0.0)
            m_dot_hpt_exh = max(m_dot_hpt_s2 - m_dot_hpt2_extr, 0.0)
            t_hpt2 = float(fp.temperature("water", P=max(p_hpt_exh, 1.0), h=max(h_hpt2, 1.0)))
            h_hpt_exh = h_hpt2

            # LOW PRESSURE TURBINE MASS FLOW
            # LPT mass-flow/pressure iteration to match condenser pressure.
            p_lpt_in = max(p_lpmain_prev, p_cond + 1000.0)
            h_lpt_in = h_lpmain_prev
            m_dot_lpt_d = max(self._safe(self.m_dot_lpt_d.v, m_dot_hpt_exh), 1.0e-6)
            m_dot_lpt_in = max(min(m_dot_hpt_exh, 1.25 * m_dot_lpt_d), 1.0e-6)

            p_lpt1 = max(self._safe(self.p_lpt1.v, 0.75 * p_lpt_in), p_cond)
            p_lpt2 = max(self._safe(self.p_lpt2.v, 0.55 * p_lpt_in), p_cond)
            p_lpt3 = max(self._safe(self.p_lpt3.v, 0.35 * p_lpt_in), p_cond)
            p_lpt4 = max(self._safe(self.p_lpt4.v, 0.20 * p_lpt_in), p_cond)
            h_lpt1 = self._safe(self.h_lpt1.v, h_lpt_in)
            h_lpt2 = self._safe(self.h_lpt2.v, h_lpt_in)
            h_lpt3 = self._safe(self.h_lpt3.v, h_lpt_in)
            h_lpt4 = self._safe(self.h_lpt4.v, h_lpt_in)

            m_dot_lpt1_extr = max(self._safe(self.m_dot_lpt1.v, 0.0), 0.0)
            m_dot_lpt2_extr = max(self._safe(self.m_dot_lpt2.v, 0.0), 0.0)
            m_dot_lpt3_extr = max(self._safe(self.m_dot_lpt3.v, 0.0), 0.0)
            m_dot_lpt4_extr = max(self._safe(self.m_dot_lpt4.v, 0.0), 0.0)

            m_dot_lpts2_d = max(m_dot_lpt_d - max(self._safe(self.m_dot_lpt1_d.v, 0.0), 0.0), 1.0e-6)
            m_dot_lpts3_d = max(m_dot_lpts2_d - max(self._safe(self.m_dot_lpt2_d.v, 0.0), 0.0), 1.0e-6)
            m_dot_lpts4_d = max(m_dot_lpts3_d - max(self._safe(self.m_dot_lpt3_d.v, 0.0), 0.0), 1.0e-6)
            m_dot_lpt_exh_d = max(m_dot_lpts4_d - max(self._safe(self.m_dot_lpt4_d.v, 0.0), 0.0), 1.0e-6)

            # Iterate with same mass flow rate until LPT stage pressures and enthalpies converge.
            for _ in range(35):
                m_dot_lpt_s1 = m_dot_lpt_in
                m_dot_lpt_s2 = max(m_dot_lpt_s1 - m_dot_lpt1_extr, 1.0e-6)
                m_dot_lpt_s3 = max(m_dot_lpt_s2 - m_dot_lpt2_extr, 1.0e-6)
                m_dot_lpt_s4 = max(m_dot_lpt_s3 - m_dot_lpt3_extr, 1.0e-6)
                m_dot_lpt_exh = max(m_dot_lpt_s4 - m_dot_lpt4_extr, 1.0e-6)

                rho_lpt_in = max(float(fp.density("water", P=max(p_lpt_in, 1.0), h=max(h_lpt_in, 1.0))), 1.0e-6)
                rho_lpt_design = max(float(fp.density("water", P=max(self._safe(self.p_lpt_in_d.v, p_lpt_in), 1.0), T=max(self._safe(self.t_lpt_in_d.v, 673.15), 273.15))), 1.0e-6)
                v_lpt_in = 1.0 / rho_lpt_in
                v_lpt_design = 1.0 / rho_lpt_design
                t_lpt_in = float(fp.temperature("water", P=max(p_lpt_in, 1.0), h=max(h_lpt_in, 1.0)))

                eta_lpt = eta_scc_lpt(
                    m_dot_lpt_in,
                    m_dot_lpt_d,
                    max(self._safe(self.lpt_parallel_sects.v, 1.0), 1.0),
                    p_lpt_in,
                    t_lpt_in,
                    v_lpt_in,
                    v_lpt_design,
                    max(self._safe(self.p_cond_d.v, p_cond), p_cond),
                )

                p_lpt1 = stodola_stage(
                    max(self._safe(self.p_lpt_in_d.v, p_lpt_in), p_cond + 1000.0),
                    max(self._safe(self.p_lpt1_d.v, p_lpt_in * 0.8), p_cond),
                    p_lpt_in,
                    max(self._safe(self.h_lpt_exh.v, h_lpt_in), 1.0),
                    h_lpt_in,
                    m_dot_lpt_d,
                    m_dot_lpt_in,
                )
                p_lpt2 = stodola_stage(max(self._safe(self.p_lpt1_d.v, p_lpt1), p_cond + 1000.0), max(self._safe(self.p_lpt2_d.v, p_lpt1 * 0.75), p_cond), p_lpt1, max(self._safe(self.h_lpt1.v, h_lpt1), 1.0), h_lpt1, m_dot_lpts2_d, m_dot_lpt_s2)
                p_lpt3 = stodola_stage(max(self._safe(self.p_lpt2_d.v, p_lpt2), p_cond + 1000.0), max(self._safe(self.p_lpt3_d.v, p_lpt2 * 0.70), p_cond), p_lpt2, max(self._safe(self.h_lpt2.v, h_lpt2), 1.0), h_lpt2, m_dot_lpts3_d, m_dot_lpt_s3)
                p_lpt4 = stodola_stage(max(self._safe(self.p_lpt3_d.v, p_lpt3), p_cond + 1000.0), max(self._safe(self.p_lpt4_d.v, p_lpt3 * 0.65), p_cond), p_lpt3, max(self._safe(self.h_lpt3.v, h_lpt3), 1.0), h_lpt3, m_dot_lpts4_d, m_dot_lpt_s4)
                p_lpt_exh = stodola_stage(max(self._safe(self.p_lpt4_d.v, p_lpt4), p_cond + 1000.0), max(self._safe(self.p_cond_d.v, p_cond), p_cond), p_lpt4, max(self._safe(self.h_lpt4.v, h_lpt4), 1.0), h_lpt4, m_dot_lpt_exh_d, m_dot_lpt_exh)

                p_lpt1 = max(p_lpt1, p_cond + 1000.0)
                p_lpt2 = max(p_lpt2, p_cond + 1000.0)
                p_lpt3 = max(p_lpt3, p_cond + 1000.0)
                p_lpt4 = max(p_lpt4, p_cond + 1000.0)
                p_lpt_exh = max(p_lpt_exh, p_cond + 1000.0)

                s_lpt_in = float(fp.entropy("water", P=max(p_lpt_in, 1.0), h=max(h_lpt_in, 1.0)))
                h_lpt_exh_s = float(fp.enthalpy("water", P=max(p_lpt_exh, 1.0), s=s_lpt_in))
                h_lpt_exh = h_lpt_in - eta_lpt * max(h_lpt_in - h_lpt_exh_s, 0.0)
                h_lpt_exh = min(max(h_lpt_exh, 1.0e5), h_lpt_in)

                s_lpt_exh = float(fp.entropy("water", P=max(p_lpt_exh, 1.0), h=max(h_lpt_exh, 1.0)))
                a0 = self._safe(self.lpt_exp_a0.v, 0.0)
                a1 = self._safe(self.lpt_exp_a1.v, 0.0)
                a2 = self._safe(self.lpt_exp_a2.v, 0.0)
                s_prime = a0 + a1 * h_lpt_exh + a2 * h_lpt_exh ** 2
                delta_s = s_lpt_exh - s_prime

                h_lpt1 = h_lpt_stage(h_lpt1, h_lpt_exh, h_lpt_in, p_lpt1, delta_s, a0, a1, a2, 1.0e-3)
                h_lpt2 = h_lpt_stage(h_lpt2, h_lpt_exh, h_lpt_in, p_lpt2, delta_s, a0, a1, a2, 1.0e-3)
                h_lpt3 = h_lpt_stage(h_lpt3, h_lpt_exh, h_lpt_in, p_lpt3, delta_s, a0, a1, a2, 1.0e-3)
                h_lpt4 = h_lpt_stage(h_lpt4, h_lpt_exh, h_lpt_in, p_lpt4, delta_s, a0, a1, a2, 1.0e-3)

                pressure_error_lpt = p_lpt_exh - p_cond
                if abs(pressure_error_lpt) <= 2000.0 and pressure_error_lpt >= 0.0:
                    break
                step = 0.02 * m_dot_lpt_d * pressure_error_lpt / max(p_lpt_in, 1.0)
                m_dot_lpt_in = max(min(m_dot_lpt_in + step, 1.25 * m_dot_lpt_d), 1.0e-6)

            m_dot_lpt_s1 = m_dot_lpt_in
            m_dot_lpt_s2 = max(m_dot_lpt_s1 - m_dot_lpt1_extr, 0.0)
            m_dot_lpt_s3 = max(m_dot_lpt_s2 - m_dot_lpt2_extr, 0.0)
            m_dot_lpt_s4 = max(m_dot_lpt_s3 - m_dot_lpt3_extr, 0.0)
            m_dot_lpt_exh = max(m_dot_lpt_s4 - m_dot_lpt4_extr, 0.0)
            t_lpt_exh = float(fp.temperature("water", P=max(p_lpt_exh, 1.0), h=max(h_lpt_exh, 1.0)))

            t_lpt1 = float(fp.temperature("water", P=max(p_lpt1, 1.0), h=max(h_lpt1, 1.0)))
            t_lpt2 = float(fp.temperature("water", P=max(p_lpt2, 1.0), h=max(h_lpt2, 1.0)))
            t_lpt3 = float(fp.temperature("water", P=max(p_lpt3, 1.0), h=max(h_lpt3, 1.0)))
            t_lpt4 = float(fp.temperature("water", P=max(p_lpt4, 1.0), h=max(h_lpt4, 1.0)))

            # Solve for work produced by turbine and temperatures of each turbine stage.
            w_hpt = max(m_dot_hpt_s1 * (h_hpmain_prev - h_hpt1) + m_dot_hpt_s2 * (h_hpt1 - h_hpt2), 0.0)
            w_lpt = max(
                m_dot_lpt_s1 * (h_lpt_in - h_lpt1)
                + m_dot_lpt_s2 * (h_lpt1 - h_lpt2)
                + m_dot_lpt_s3 * (h_lpt2 - h_lpt3)
                + m_dot_lpt_s4 * (h_lpt3 - h_lpt4)
                + m_dot_lpt_exh * (h_lpt4 - h_lpt_exh),
                0.0,
            )
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
            t_lpt1 = 0.0
            t_lpt2 = 0.0
            t_lpt3 = 0.0
            t_lpt4 = 0.0
            w_dot_total = 0.0

        # BYPASS & DRAIN VALVE FLOWS + ENTHALPIES
        # Solve for the mass flow rates and enthalpies through all bypass and drain valves.
        p_atm = 101325.0
        p_da = max(self._safe(self.p_da.v, p_aux_prev), 1.0)

        hp_bypass_cv = pb_cv_data(self.hp_bypass_vt.v, self.hp_bypass_d.v, hp_bypass_vpi)
        hp_aux_cv = pb_cv_data(self.hp_aux_vt.v, self.hp_aux_d.v, hp_aux_vpi)
        hp_warmup_cv = pb_cv_data(self.hp_warmup_vt.v, self.hp_warmup_d.v, hp_warmup_vpi)
        hp_drain_cv = pb_cv_data(self.hp_drain_vt.v, self.hp_drain_d.v, hp_drain_vpi)
        aux_da_cv = pb_cv_data(self.aux_da_vt.v, self.aux_da_d.v, aux_da_vpi)
        lp_bypass_cv = pb_cv_data(self.lp_bypass_vt.v, self.lp_bypass_d.v, lp_bypass_vpi)
        lp_aux_cv = pb_cv_data(self.lp_aux_vt.v, self.lp_aux_d.v, lp_aux_vpi)
        lp_warmup_cv = pb_cv_data(self.lp_warmup_vt.v, self.lp_warmup_d.v, lp_warmup_vpi)
        lp_drain_cv = pb_cv_data(self.lp_drain_vt.v, self.lp_drain_d.v, lp_drain_vpi)

        m_dot_hp_bypass = valve_massflow(hp_bypass_cv, p_sgt_in, h_sgt_in, p_lpmain_prev) if hp_bypass_vpi > 0.0 else 0.0
        h_hp_bypass = h_sgt_in
        m_dot_hp_aux = valve_massflow(hp_aux_cv, p_sgt_in, h_sgt_in, p_aux_prev) if hp_aux_vpi > 0.0 else 0.0
        h_hp_aux = h_sgt_in
        m_dot_hp_warmup = valve_massflow(hp_warmup_cv, p_hpmain_prev, h_hpmain_prev, p_cond) if hp_warmup_vpi > 0.0 else 0.0
        h_hp_warmup = h_hpmain_prev

        area_hpmain = math.pi * max(self.d_hpmain.v, 1.0e-3) ** 2 / 4.0
        vol_hpmain = area_hpmain * max(self.length_hpmain.v, 1.0)
        m_hpmain_prev = vol_hpmain * max(float(fp.density("water", P=max(p_hpmain_prev, 1.0), h=max(h_hpmain_prev, 1.0))), 1.0e-6)
        x_hp = min(max(float(fp.quality("water", P=max(p_hpmain_prev, 1.0), h=max(h_hpmain_prev, 1.0))), 0.0), 1.0)
        if hp_drain_vpi > 0.0 and x_hp < 1.0:
            h_sat_f_hp = float(fp.enthalpy("water", P=max(p_hpmain_prev, 1.0), Q=0.0))
            m_dot_hp_drain_cap = max(m_hpmain_prev * (1.0 - x_hp) / ts_sec, 0.0)
            m_dot_hp_drain = min(valve_massflow(hp_drain_cv, p_hpmain_prev, h_sat_f_hp, p_atm), m_dot_hp_drain_cap)
            h_hp_drain = h_sat_f_hp
        else:
            m_dot_hp_drain = 0.0
            h_hp_drain = 0.0

        m_dot_aux_da = valve_massflow(aux_da_cv, p_aux_prev, h_aux_prev, p_da) if aux_da_vpi > 0.0 else 0.0
        h_aux_da = h_aux_prev if m_dot_aux_da > 0.0 else 0.0

        m_dot_lp_bypass = valve_massflow(lp_bypass_cv, p_lpmain_prev, h_lpmain_prev, p_cond) if lp_bypass_vpi > 0.0 else 0.0
        h_lp_bypass = h_lpmain_prev
        m_dot_lp_aux = valve_massflow(lp_aux_cv, p_lpmain_prev, h_lpmain_prev, p_aux_prev) if lp_aux_vpi > 0.0 else 0.0
        h_lp_aux = h_lpmain_prev
        m_dot_lp_warmup = valve_massflow(lp_warmup_cv, p_lpmain_prev, h_lpmain_prev, p_cond) if lp_warmup_vpi > 0.0 else 0.0
        h_lp_warmup = h_lpmain_prev

        area_lpmain = math.pi * max(self.d_lpmain.v, 1.0e-3) ** 2 / 4.0
        vol_lpmain = area_lpmain * max(self.length_lpmain.v, 1.0)
        m_lpmain_prev = vol_lpmain * max(float(fp.density("water", P=max(p_lpmain_prev, 1.0), h=max(h_lpmain_prev, 1.0))), 1.0e-6)
        x_lp = min(max(float(fp.quality("water", P=max(p_lpmain_prev, 1.0), h=max(h_lpmain_prev, 1.0))), 0.0), 1.0)
        if lp_drain_vpi > 0.0 and x_lp < 1.0:
            h_sat_f_lp = float(fp.enthalpy("water", P=max(p_lpmain_prev, 1.0), Q=0.0))
            m_dot_lp_drain_cap = max(m_lpmain_prev * (1.0 - x_lp) / ts_sec, 0.0)
            m_dot_lp_drain = min(valve_massflow(lp_drain_cv, p_lpmain_prev, h_sat_f_lp, p_atm), m_dot_lp_drain_cap)
            h_lp_drain = h_sat_f_lp
        else:
            m_dot_lp_drain = 0.0
            h_lp_drain = 0.0

        # STEAM SEPARATOR CODE
        # Combine high pressure turbine exhaust with high pressure bypass mass flow.
        m_dot_ss_in = m_dot_hpt_exh + m_dot_hp_bypass
        if m_dot_ss_in > 0.0:
            h_ss_in = (m_dot_hpt_exh * h_hpt_exh + m_dot_hp_bypass * h_sgt_in) / m_dot_ss_in
            p_ss_sat = min(max(p_lpmain_prev, 1.0), 2.2e7)
            h_sat_g = float(fp.enthalpy("water", P=p_ss_sat, Q=1.0))
            h_sat_f = float(fp.enthalpy("water", P=p_ss_sat, Q=0.0))
            t_sat = float(fp.temperature("water", P=p_ss_sat, Q=0.0))
            rho_sat_g = max(float(fp.density("water", P=p_ss_sat, Q=1.0)), 1.0e-6)
            rho_sat_f = max(float(fp.density("water", P=p_ss_sat, Q=0.0)), 1.0e-6)

            if h_ss_in >= h_sat_g:
                m_dot_ss_drain = 0.0
                vol_dot_ss_drain = 0.0
                h_ss_drain = 0.0
                t_ss_drain = 0.0

                m_dot_ss_steam = m_dot_ss_in
                rho_ss_in = max(float(fp.density("water", P=max(p_lpmain_prev, 1.0), h=max(h_ss_in, 1.0))), 1.0e-6)
                vol_dot_ss_steam = m_dot_ss_steam / rho_ss_in
                h_ss_steam = h_ss_in
                t_ss_steam = float(fp.temperature("water", P=max(p_lpmain_prev, 1.0), h=max(h_ss_steam, 1.0)))
            else:
                x_ss_in = min(max((h_ss_in - h_sat_f) / max(h_sat_g - h_sat_f, 1.0e-6), 0.0), 1.0)
                x_ss_in = round(x_ss_in, 3)

                m_dot_ss_drain = m_dot_ss_in * (1.0 - x_ss_in)
                vol_dot_ss_drain = m_dot_ss_drain / rho_sat_f
                h_ss_drain = h_sat_f
                t_ss_drain = t_sat

                m_dot_ss_steam = m_dot_ss_in - m_dot_ss_drain
                vol_dot_ss_steam = m_dot_ss_steam / rho_sat_g
                h_ss_steam = h_sat_g
                t_ss_steam = t_sat
        else:
            m_dot_ss_drain = 0.0
            vol_dot_ss_drain = 0.0
            m_dot_ss_steam = 0.0
            vol_dot_ss_steam = 0.0
            h_ss_drain = 0.0
            h_ss_steam = 0.0
            t_ss_drain = 0.0
            t_ss_steam = 0.0

        # REHEATER CODE
        if m_dot_ss_steam > 0.0 and m_dot_htf > 0.0:
            cp_htf = max(float(self._inc_props.specheat(fluid_name, float(t_htf_in), float(p_htf_in))), 1.0)
            cp_steam = cp_water(p_lpmain_prev, 0.5 * (t_htf_in + t_ss_steam))
            hx_tube_od = self._safe(self.hx_tube_od.v, 0.05)
            hx_tube_th = self._safe(self.hx_tube_th.v, 0.002)
            hx_length = self._safe(self.hx_length.v, 1.0)
            hx_no_shell = max(self._safe(self.hx_no_shell.v, 1.0), 1.0)
            hx_no_tubes = max(self._safe(self.hx_no_tubes.v, 1.0), 1.0)
            hx_ua_d = max(self._safe(self.hx_ua_d.v, 0.0), 0.0)
            m_dot_htf_hx_d = max(self._safe(self.m_dot_htf_hx_d.v, m_dot_htf), 1.0e-9)
            hx_exp = max(self._safe(self.hx_exp.v, 0.0), 0.0)

            area_surface = (
                math.pi
                * max(hx_tube_od - 2.0 * hx_tube_th, 1.0e-6)
                * max(hx_length, 1.0e-6)
                * hx_no_shell
                * hx_no_tubes
            )
            ua_rated = hx_ua_d * area_surface
            ua_od = ua_rated * (m_dot_htf / m_dot_htf_hx_d) ** hx_exp

            cap_steam = m_dot_ss_steam * cp_steam
            cap_htf = m_dot_htf * cp_htf
            cap_min = max(min(cap_steam, cap_htf), 1.0e-9)
            cap_max = max(cap_steam, cap_htf)
            cr = max(cap_min / max(cap_max, 1.0e-9), 0.001)
            ntu_od = ua_od / cap_min

            exp_arg = math.exp(-ntu_od * math.sqrt(1.0 + cr ** 2))
            eta_1pass_den = 1.0 + cr + math.sqrt(1.0 + cr ** 2) * ((1.0 + exp_arg) / max(1.0 - exp_arg, 1.0e-9))
            eta_1pass = 2.0 / max(eta_1pass_den, 1.0e-9)

            ratio = (1.0 - eta_1pass * cr) / max(1.0 - eta_1pass, 1.0e-9)
            ratio_n = ratio ** hx_no_shell
            eta_od = (ratio_n - 1.0) / max(ratio_n - cr, 1.0e-9)
            eta_od = min(max(eta_od, 0.0), 1.0)

            # Check that HTF temperature is higher than steam temperature entering reheater.
            if t_ss_steam < t_htf_in:
                h_hx_out_s = float(fp.enthalpy("water", P=max(p_lpmain_prev, 1.0), T=max(t_htf_in, 273.15)))
                q_dot_hx = min(
                    eta_od * m_dot_ss_steam * max(h_hx_out_s - h_ss_steam, 0.0),
                    eta_od * m_dot_htf * cp_htf * max(t_htf_in - t_ss_steam, 0.0),
                )
                h_hx_out = (m_dot_ss_steam * h_ss_steam + q_dot_hx) / max(m_dot_ss_steam, 1.0e-9)
                t_htf_out = (m_dot_htf * cp_htf * t_htf_in - q_dot_hx) / max(m_dot_htf * cp_htf, 1.0e-9)
            else:
                h_hx_out_s = float(fp.enthalpy("water", P=max(p_lpmain_prev, 1.0), T=max(t_htf_in, 273.15)))
                q_dot_hx = min(
                    eta_od * m_dot_ss_steam * max(h_ss_steam - h_hx_out_s, 0.0),
                    eta_od * m_dot_htf * cp_htf * max(t_ss_steam - t_htf_in, 0.0),
                )
                h_hx_out = (m_dot_ss_steam * h_ss_steam - q_dot_hx) / max(m_dot_ss_steam, 1.0e-9)
                t_htf_out = (m_dot_htf * cp_htf * t_htf_in + q_dot_hx) / max(m_dot_htf * cp_htf, 1.0e-9)

            rho_htf = max(float(self._inc_props.density(fluid_name, float(t_htf_out), float(p_htf_in))), 1.0e-6)
            vol_dot_htf = m_dot_htf / max(rho_htf, 1.0e-9)

            if not math.isfinite(h_hx_out):
                h_hx_out = h_ss_steam
            if not math.isfinite(t_htf_out):
                t_htf_out = t_htf_in
            if not math.isfinite(q_dot_hx):
                q_dot_hx = 0.0

            t_hx_out = float(fp.temperature("water", P=max(p_lpmain_prev, 1.0), h=max(h_hx_out, 1.0)))
            rho_steam = max(float(fp.density("water", P=max(p_lpmain_prev, 1.0), h=max(h_hx_out, 1.0))), 0.6)
            vol_dot_ss_steam = m_dot_ss_steam / rho_steam
        else:
            q_dot_hx = 0.0
            h_hx_out = h_ss_steam
            t_hx_out = t_ss_steam
            t_htf_out = t_htf_in
            vol_dot_htf = m_dot_htf / max(float(self._inc_props.density(fluid_name, float(t_htf_out), float(p_htf_in))), 1.0e-9)
            eta_od = 0.0

        # HP main / LP main / AUX piping code
        # Dynamic pipe-state updates with mass-energy balance and pipe thermal inertia.
        def _pipe_update(p_prev, t_prev, t_pipe_prev, diameter, length, mc_pipe, m_dot_in, h_in, m_dot_out, h_out, ff_guess):
            diameter = max(diameter, 1.0e-4)
            length = max(length, 1.0)
            area = math.pi * diameter * diameter / 4.0
            volume = area * length
            m_dot_ave = 0.5 * (m_dot_in + m_dot_out)
            dh_fd = 25000.0
            dp_fd = 25000.0

            p_state = max(p_prev, 700.0)
            h_state = float(fp.enthalpy("water", P=max(p_state, 1.0), T=max(t_prev, 273.15)))
            t_pipe_state = t_pipe_prev
            ff_state = max(ff_guess, 1.0e-4)

            n_sub = max(int(math.ceil(ts_sec / 0.1)), 1)
            t_crit = ts_sec / n_sub

            def _rhs(p_local, h_local, t_pipe_local, ff_local):
                p_local = max(p_local, 700.0)
                h_local = max(h_local, 1.0)
                p_prop = min(max(p_local, 700.0), 3.0e7)
                h_prop = min(max(h_local, 1.0), 5.0e6)
                rho_local = max(float(fp.density("water", P=max(p_prop, 1.0), h=max(h_prop, 1.0))), 1.0e-6)
                t_fluid = float(fp.temperature("water", P=max(p_prop, 1.0), h=max(h_prop, 1.0)))
                u_local = h_prop - p_prop / max(rho_local, 1.0e-9)

                h_bar = convection_dynamicpipe(p_prop, h_prop, diameter, m_dot_ave, ff_local, t_pipe_local)
                q_dot_max = mc_pipe * (t_fluid - t_pipe_local) / t_crit if mc_pipe > 0.0 else 0.0
                q_dot_raw = h_bar * math.pi * diameter * length * (t_fluid - t_pipe_local)
                if q_dot_max > 0.0:
                    q_dot = min(q_dot_raw, q_dot_max)
                else:
                    q_dot = max(q_dot_raw, q_dot_max)

                drho_dh = drhodhcp(p_prop, h_prop, dh_fd)
                drho_dp = drhodpch(p_prop, h_prop, dp_fd)
                du_dh = dudhcp(p_prop, h_prop, dh_fd)
                du_dp = dudpch(p_prop, h_prop, dp_fd)

                denom = volume * rho_local * (drho_dh * du_dp - drho_dp * du_dh)
                if abs(denom) > 1.0e-16 and abs(drho_dp) > 1.0e-16:
                    dh_dt = -(
                        h_in * drho_dp * m_dot_in
                        - h_out * drho_dp * m_dot_out
                        - u_local * drho_dp * m_dot_in
                        + u_local * drho_dp * m_dot_out
                        - du_dp * m_dot_in * rho_local
                        + du_dp * m_dot_out * rho_local
                        - q_dot * drho_dp
                    ) / denom
                    dp_dt = ((m_dot_in - m_dot_out) / volume - dh_dt * drho_dh) / drho_dp
                else:
                    mass_local = max(rho_local * volume, 1.0e-6)
                    dh_dt = (m_dot_in * (h_in - h_local) - m_dot_out * (h_out - h_local) - q_dot) / mass_local
                    dp_dt = 0.0

                dt_pipe_dt = q_dot / mc_pipe if mc_pipe > 0.0 else 0.0
                return dp_dt, dh_dt, dt_pipe_dt

            for _ in range(n_sub):
                k1_p, k1_h, k1_tp = _rhs(p_state, h_state, t_pipe_state, ff_state)
                k2_p, k2_h, k2_tp = _rhs(
                    p_state + 0.5 * t_crit * k1_p,
                    h_state + 0.5 * t_crit * k1_h,
                    t_pipe_state + 0.5 * t_crit * k1_tp,
                    ff_state,
                )
                k3_p, k3_h, k3_tp = _rhs(
                    p_state + 0.5 * t_crit * k2_p,
                    h_state + 0.5 * t_crit * k2_h,
                    t_pipe_state + 0.5 * t_crit * k2_tp,
                    ff_state,
                )
                k4_p, k4_h, k4_tp = _rhs(
                    p_state + t_crit * k3_p,
                    h_state + t_crit * k3_h,
                    t_pipe_state + t_crit * k3_tp,
                    ff_state,
                )

                p_state = max(p_state + (t_crit / 6.0) * (k1_p + 2.0 * k2_p + 2.0 * k3_p + k4_p), 700.0)
                h_state = max(h_state + (t_crit / 6.0) * (k1_h + 2.0 * k2_h + 2.0 * k3_h + k4_h), 1.0)
                t_pipe_state = t_pipe_state + (t_crit / 6.0) * (k1_tp + 2.0 * k2_tp + 2.0 * k3_tp + k4_tp)

                t_now = float(fp.temperature("water", P=max(p_state, 1.0), h=max(h_state, 1.0)))
                mu_now = max(float(fp.viscosity("water", P=max(p_state, 1.0), h=max(h_state, 1.0))), 1.0e-6)
                re_now = max(4.0 * abs(m_dot_ave) / (math.pi * diameter * mu_now), 1.0)
                if re_now < 2300.0:
                    ff_state = min(max(64.0 / re_now, 1.0e-4), 0.5)
                else:
                    ff_state = min(max(0.3164 / (re_now ** 0.25), 1.0e-4), 0.2)

            p_new = p_state
            h_new = h_state
            t_pipe_new = t_pipe_state
            t_new = float(fp.temperature("water", P=max(p_new, 1.0), h=max(h_new, 1.0)))

            mu = max(float(fp.viscosity("water", P=max(p_new, 1.0), h=max(h_new, 1.0))), 1.0e-6)
            re = max(4.0 * abs(m_dot_ave) / (math.pi * diameter * mu), 1.0)
            if re < 2300.0:
                ff = min(max(64.0 / re, 1.0e-4), 0.5)
            else:
                ff = min(max(0.3164 / (re ** 0.25), 1.0e-4), 0.2)

            x_new = min(max(float(fp.quality("water", P=max(p_new, 1.0), h=max(h_new, 1.0))), 0.0), 1.0)

            return p_new, t_new, x_new * 100.0, t_pipe_new, h_new, ff

        p_hpmain, t_hpmain, x_hpmain_pct, t_hppipe, h_hpmain, ff_hp = _pipe_update(
            p_hpmain_prev,
            t_hpmain_prev,
            t_hppipe_prev,
            self._safe(self.d_hpmain.v, 0.5),
            self._safe(self.length_hpmain.v, 100.0),
            max(self._safe(self.mc_hpmain_pipe.v, 1.0), 1.0),
            max(m_dot_sgt_in - m_dot_hp_bypass - m_dot_hp_aux, 0.0),
            h_sgt_in,
            m_dot_hp_warmup + m_dot_hp_drain + m_dot_hpt_in,
            h_hpmain_prev,
            self._safe(self.friction_factor_hp.v, 0.02),
        )

        p_lpmain, t_lpmain, x_lpmain_pct, t_lppipe, h_lpmain, ff_lp = _pipe_update(
            p_lpmain_prev,
            t_lpmain_prev,
            t_lppipe_prev,
            self._safe(self.d_lpmain.v, 0.5),
            self._safe(self.length_lpmain.v, 120.0),
            max(self._safe(self.mc_lpmain_pipe.v, 1.0), 1.0),
            m_dot_ss_steam,
            h_hx_out,
            m_dot_lp_bypass + m_dot_lp_warmup + m_dot_lp_aux + m_dot_lp_drain + m_dot_lpt_s1,
            h_lpmain_prev,
            self._safe(self.friction_factor_lp.v, 0.02),
        )

        m_dot_ts_req = max(self._safe(self.m_dot_turbine_seals.v, 0.0), 0.0)
        p_ts_req = max(self._safe(self.p_ts_req.v, 1.0e5), 1.0)
        p_ts_min = 0.7 * p_ts_req
        if turbine_on != 1.0 and p_aux_prev <= p_ts_req:
            m_dot_ts = max(m_dot_ts_req * (p_aux_prev - p_ts_min) / max(p_ts_req - p_ts_min, 1.0), 0.0)
        else:
            m_dot_ts = m_dot_ts_req if turbine_on != 1.0 else 0.0

        p_aux, t_aux, x_aux_pct, t_auxpipe, h_aux, ff_aux = _pipe_update(
            p_aux_prev,
            t_aux_prev,
            t_auxpipe_prev,
            self._safe(self.d_auxline.v, 0.3),
            self._safe(self.length_auxline.v, 80.0),
            max(self._safe(self.mc_aux_pipe.v, 1.0), 1.0),
            m_dot_lp_aux + m_dot_hp_aux,
            (m_dot_lp_aux * h_lp_aux + m_dot_hp_aux * h_hp_aux) / max(m_dot_lp_aux + m_dot_hp_aux, 1.0e-9),
            m_dot_ts + m_dot_aux_da,
            h_aux_prev,
            self._safe(self.friction_factor_aux.v, 0.02),
        )

        # Turbine-seal trip logic (output 95)
        m_dot_turbine_seals = max(self._safe(self.m_dot_turbine_seals.v, 0.0), 0.0)
        if turbine_on != 1.0:
            if p_aux > p_ts_req:
                trip_ts = 0.0
            else:
                trip_ts = 1.0
        else:
            trip_ts = 0.0

        # Smoothed steam-drum demand (output 12)
        m_dot_sgt_req_raw = m_dot_hp_bypass + m_dot_hp_aux + m_dot_hpt_in
        m_dot_sgt_prev = self._safe(self.m_dot_sgt_req.v, m_dot_sgt_req_raw)
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
            h_da = (m_dot_aux_da * h_aux_da + m_dot_ss_drain * h_ss_drain) / m_dot_da
        else:
            h_da = 0.0

        # Main outputs 12-95 mapping
        self.m_dot_sgt_req.v = m_dot_sgt_req
        self.m_dot_cond.v = m_dot_cond
        self.h_cond.v = h_cond
        self.m_dot_da.v = m_dot_da
        self.h_da.v = h_da
        self.w_dot_total.v = w_dot_total

        self.m_dot_htf_out.v = m_dot_htf
        self.vol_dot_htf_out.v = vol_dot_htf
        self.htf_p_in_out.v = p_htf_in
        self.htf_t_out.v = t_htf_out

        self.m_dot_hpt_in.v = m_dot_hpt_in
        self.p_gs_out.v = p_hpt1
        self.h_gs_out.v = h_hpt1
        self.m_dot_hpts1.v = m_dot_hpt_in
        self.p_hpt1.v = p_hpt1
        self.h_hpt1.v = h_hpt1
        self.m_dot_hpts2.v = m_dot_hpt_s2
        self.p_hpt2.v = p_hpt2
        self.h_hpt2.v = h_hpt2
        self.m_dot_hpt_exh.v = m_dot_hpt_exh
        self.p_hpt_exh.v = p_hpt_exh
        self.h_hpt_exh.v = h_hpt_exh
        self.t_hpt2.v = t_hpt2

        self.m_dot_ss_drain.v = m_dot_ss_drain
        self.vol_dot_ss_drain.v = vol_dot_ss_drain
        self.p_drain.v = p_lpmain
        self.h_ss_drain.v = h_ss_drain
        self.t_ss_drain.v = t_ss_drain

        self.m_dot_ss_steam.v = m_dot_ss_steam
        self.vol_dot_ss_steam.v = vol_dot_ss_steam
        self.p_steam.v = p_lpmain
        self.h_ss_steam.v = h_ss_steam
        self.t_ss_steam.v = t_ss_steam

        self.m_dot_steam_out.v = m_dot_ss_steam
        self.vol_dot_steam_out.v = vol_dot_ss_steam
        self.t_steam_out.v = t_hx_out
        self.p_steam_out.v = p_lpmain
        self.q_dot_hx.v = q_dot_hx
        self.eta_od.v = eta_od

        self.m_dot_lpt_stage1.v = m_dot_lpt_s1
        self.p_lpt1.v = p_lpt1
        self.t_lpt1.v = t_lpt1
        self.h_lpt1.v = h_lpt1
        self.m_dot_lpt_stage2.v = m_dot_lpt_s2
        self.p_lpt2.v = p_lpt2
        self.t_lpt2.v = t_lpt2
        self.h_lpt2.v = h_lpt2
        self.m_dot_lpt_stage3.v = m_dot_lpt_s3
        self.p_lpt3.v = p_lpt3
        self.t_lpt3.v = t_lpt3
        self.h_lpt3.v = h_lpt3
        self.m_dot_lpt_stage4.v = m_dot_lpt_s4
        self.t_lpt4.v = t_lpt4
        self.p_lpt4.v = p_lpt4
        self.h_lpt4.v = h_lpt4
        self.m_dot_lpt_exh.v = m_dot_lpt_exh
        self.vol_dot_lpt_exh.v = m_dot_lpt_exh / max(float(fp.density("water", P=max(max(p_lpt_exh, p_cond), 1.0), h=max(h_lpt_exh, 1.0))), 1.0e-9)
        self.t_lpt_exh.v = t_lpt_exh
        self.h_lpt_exh.v = h_lpt_exh

        self.p_hpmain.v = p_hpmain
        self.t_hpmain.v = t_hpmain
        self.x_hpmain.v = x_hpmain_pct
        self.t_hppipe.v = t_hppipe

        self.p_lpmain.v = p_lpmain
        self.t_lpmain.v = t_lpmain
        self.x_lpmain.v = x_lpmain_pct
        self.t_lppipe.v = t_lppipe

        self.p_aux.v = p_aux
        self.t_aux.v = t_aux
        self.x_aux.v = x_aux_pct
        self.t_auxpipe.v = t_auxpipe

        self.m_dot_hp_bypass.v = m_dot_hp_bypass
        self.m_dot_hp_aux.v = m_dot_hp_aux
        self.m_dot_hp_drain.v = m_dot_hp_drain
        self.m_dot_hp_warmup.v = m_dot_hp_warmup
        self.m_dot_lp_aux.v = m_dot_lp_aux
        self.m_dot_lp_bypass.v = m_dot_lp_bypass
        self.m_dot_lp_drain.v = m_dot_lp_drain
        self.m_dot_lp_warmup.v = m_dot_lp_warmup
        self.m_dot_aux_da.v = m_dot_aux_da

        self.friction_factor_hp.v = ff_hp
        self.friction_factor_lp.v = ff_lp
        self.friction_factor_aux.v = ff_aux
        self.trip_turbine_seals.v = trip_ts

        # End-of-timestep alarm/trip channels (96-123)
        if turbine_on != 1.0:
            self.alarm_hpt_superheat.v = 0.0
            self.trip_hpt_superheat.v = 0.0
            self.alarm_hpt_hightemp.v = 0.0
            self.timedtrip_hpt_hightemp.v = 0.0
            self.trip_hpt_hightemp.v = 0.0
            self.alarm_hpt_exhpres.v = 0.0
            self.trip_hpt_exhpres.v = 0.0
            self.alarm_lpt_superheat.v = 0.0
            self.trip_lpt_superheat.v = 0.0
            self.alarm_lpt_hightemp.v = 0.0
            self.timedtrip_lpt_hightemp.v = 0.0
            self.trip_lpt_hightemp.v = 0.0
        else:
            hpt_sh_alarm_fl = self._safe(self.hpt_sh_alarm_fl.v, 50.0)
            hpt_sh_trip_fl = self._safe(self.hpt_sh_trip_fl.v, 20.0)
            hpt_sh_alarm_pl = self._safe(self.hpt_sh_alarm_pl.v, 50.0)
            hpt_sh_trip_pl = self._safe(self.hpt_sh_trip_pl.v, 20.0)
            partial_load = self._safe(self.partial_load.v, 1.0e8)

            hpt_high_t_alarm = self._safe(self.hpt_hightemp_alarm.v, 850.0)
            hpt_high_t_timed = self._safe(self.hpt_hightemp_timedtrip.v, 900.0)
            hpt_high_t_trip = self._safe(self.hpt_hightemp_trip.v, 950.0)
            hpt_exh_p_alarm = self._safe(self.hpt_exhpres_alarm.v, 3.0e6)
            hpt_exh_p_trip = self._safe(self.hpt_exhpres_trip.v, 4.0e6)

            lpt_sh_alarm = self._safe(self.lpt_sh_alarm.v, 50.0)
            lpt_sh_trip = self._safe(self.lpt_sh_trip.v, 20.0)
            lpt_high_t_alarm = self._safe(self.lpt_hightemp_alarm.v, 750.0)
            lpt_high_t_timed = self._safe(self.lpt_hightemp_timedtrip.v, 800.0)
            lpt_high_t_trip = self._safe(self.lpt_hightemp_trip.v, 850.0)

            t_sat_hp = float(fp.temperature("water", P=min(max(self.p_hpmain.v, 1.0), 2.2e7), Q=0.0))
            superheat_hp = self.t_hpmain.v - t_sat_hp
            if self.w_dot_total.v > partial_load:
                alarm, trip = self._alarm_trip(superheat_hp, hpt_sh_alarm_fl, hpt_sh_trip_fl, high=False)
            else:
                alarm, trip = self._alarm_trip(superheat_hp, hpt_sh_alarm_pl, hpt_sh_trip_pl, high=False)
            self.alarm_hpt_superheat.v = alarm
            self.trip_hpt_superheat.v = trip

            alarm, timed_trip, trip = self._alarm_timed_trip(self.t_hpmain.v, hpt_high_t_alarm, hpt_high_t_timed, hpt_high_t_trip)
            self.alarm_hpt_hightemp.v = alarm
            self.timedtrip_hpt_hightemp.v = timed_trip
            self.trip_hpt_hightemp.v = trip

            alarm, trip = self._alarm_trip(self.p_lpmain.v, hpt_exh_p_alarm, hpt_exh_p_trip, high=True)
            self.alarm_hpt_exhpres.v = alarm
            self.trip_hpt_exhpres.v = trip

            t_sat_lp = float(fp.temperature("water", P=min(max(self.p_lpmain.v, 1.0), 2.2e7), Q=0.0))
            superheat_lp = self.t_lpmain.v - t_sat_lp
            alarm, trip = self._alarm_trip(superheat_lp, lpt_sh_alarm, lpt_sh_trip, high=False)
            self.alarm_lpt_superheat.v = alarm
            self.trip_lpt_superheat.v = trip

            alarm, timed_trip, trip = self._alarm_timed_trip(self.t_lpmain.v, lpt_high_t_alarm, lpt_high_t_timed, lpt_high_t_trip)
            self.alarm_lpt_hightemp.v = alarm
            self.timedtrip_lpt_hightemp.v = timed_trip
            self.trip_lpt_hightemp.v = trip

        # Reheater alarms/trips
        hi_htf_in_alarm = self._safe(self.hx_hightemp_alarm.v, 900.0)
        hi_htf_in_trip = self._safe(self.hx_hightemp_trip.v, 950.0)
        lo_htf_out_alarm = self._safe(self.hx_lowtemp_alarm.v, 450.0)
        lo_htf_out_trip = self._safe(self.hx_lowtemp_trip.v, 400.0)
        hi_htf_flow_alarm = self._safe(self.hx_highflow_alarm.v, 1.0e9)
        hi_htf_flow_trip = self._safe(self.hx_highflow_trip.v, 1.0e9)
        hi_htf_p_alarm = self._safe(self.hx_highpressure_alarm.v, 1.0e9)
        hi_htf_p_trip = self._safe(self.hx_highpressure_trip.v, 1.0e9)
        hi_dt_alarm = self._safe(self.hx_highdt_alarm.v, 200.0)
        hi_dt_trip = self._safe(self.hx_highdt_trip.v, 250.0)

        alarm, trip = self._alarm_trip(t_htf_in, hi_htf_in_alarm, hi_htf_in_trip, high=True)
        self.alarm_hx_hightemp_in.v = alarm
        self.trip_hx_hightemp_in.v = trip

        alarm, trip = self._alarm_trip(t_htf_out, lo_htf_out_alarm, lo_htf_out_trip, high=False)
        self.alarm_hx_lowtemp_out.v = alarm
        self.trip_hx_lowtemp_out.v = trip

        alarm, trip = self._alarm_trip(m_dot_htf, hi_htf_flow_alarm, hi_htf_flow_trip, high=True)
        self.alarm_hx_highflow.v = alarm
        self.trip_hx_highflow.v = trip

        alarm, trip = self._alarm_trip(p_htf_in, hi_htf_p_alarm, hi_htf_p_trip, high=True)
        self.alarm_hx_highpressure.v = alarm
        self.trip_hx_highpressure.v = trip

        delta_t = t_htf_in - self.t_ss_steam.v
        alarm, trip = self._alarm_trip(delta_t, hi_dt_alarm, hi_dt_trip, high=True)
        self.alarm_hx_highdt.v = alarm
        self.trip_hx_highdt.v = trip

        # Heating-rate channels (108,109,120-123)
        # Use end-of-timestep update semantics to avoid counting inner convergence iterations.
        if is_converged:
            self._ensure_rate_buffers(n_int, t_htf_in, self.t_ss_steam.v)
            hr_htf = self._compute_heating_rate(self._htf_rate_hist, t_htf_in, ts_sec)
            hr_steam = self._compute_heating_rate(self._steam_rate_hist, self.t_ss_steam.v, ts_sec)

            hr_htf_alarm = self._safe(self.hx_htfhr_alarm.v, 999.0)
            hr_htf_trip = self._safe(self.hx_htfhr_trip.v, 9999.0)
            hr_steam_alarm = self._safe(self.hx_steamhr_alarm.v, 999.0)
            hr_steam_trip = self._safe(self.hx_steamhr_trip.v, 9999.0)

            self.hr_htf_in.v = hr_htf
            self.hr_steam.v = hr_steam
            self.alarm_hx_htf_hr.v = 1.0 if abs(hr_htf) >= hr_htf_alarm else 0.0
            self.trip_hx_htf_hr.v = 1.0 if abs(hr_htf) >= hr_htf_trip else 0.0
            self.alarm_hx_steam_hr.v = 1.0 if abs(hr_steam) >= hr_steam_alarm else 0.0
            self.trip_hx_steam_hr.v = 1.0 if abs(hr_steam) >= hr_steam_trip else 0.0

            self._advance_history(self._htf_rate_hist, t_htf_in)
            self._advance_history(self._steam_rate_hist, self.t_ss_steam.v)
