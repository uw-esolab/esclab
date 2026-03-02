"""
esol6015_helpers.py

Python conversion of helper functions from ESOL6015_myfunctions.f90
(TRNSYS Fortran module ESOL6015_myfunctions).

Units follow the SI convention unless otherwise noted:
  Pressures  : Pa
  Enthalpies : J/kg
  Temperatures: K
  Densities  : kg/m^3
  Mass flows : kg/s
  Volumes    : m^3
  Lengths    : m

eeslib property calls require:
  P in kPa  (divide Pa by 1000)
  h in kJ/kg (divide J/kg by 1000)
  Returned h/s in kJ/kg (multiply by 1000 to get J/kg)
"""

import math
import numpy as np
from eeslib import fluid_properties as fp
from esclab.components.esol_properties import Incompressible

# Module-level Incompressible instance for HTF property lookups (specheat, density)
_incompressible = Incompressible()

# ---------------------------------------------------------------------------
# Module-level lookup tables for PB_CV_data
# ---------------------------------------------------------------------------

# Concentric Butterfly Valve diameters [inches] and positions [0-1]
_D_CONCENTRIC = np.array([
    2., 3., 4., 5., 6., 8., 12., 16., 20., 24., 28., 32., 36., 40.
])

_POS_CONCENTRIC = np.array([
    0., 1./9., 2./9., 3./9., 4./9., 5./9., 6./9., 7./9., 8./9., 1.
])

# CV table for Concentric Butterfly: shape (14 diameters, 10 positions)
# Source: https://www.valvesonline.com.au/references/flow-rates/butterfly-valves/
_CV_CONCENTRIC = np.array([
    [0.,   0.1,   5.,   12.,   24.,   45.,   64.,   90.,   125.,   135.],    # 2 inch
    [0.,   0.3,  12.,   22.,   39.,   70.,  116.,  183.,   275.,   302.],    # 3 inch
    [0.,   0.5,  17.,   36.,   78.,  139.,  230.,  364.,   546.,   600.],    # 4 inch
    [0.,   0.8,  29.,   61.,  133.,  237.,  392.,  620.,   930.,  1022.],    # 5 inch
    [0.,   2.,   45.,   95.,  205.,  366.,  605.,  958.,  1437.,  1579.],    # 6 inch
    [0.,   3.,   89.,  188.,  408.,  727., 1202., 1903.,  2854.,  3136.],    # 8 inch
    [0.,   5.,  234.,  495., 1072., 1911., 3162., 5005.,  7507.,  8250.],    # 12 inch
    [0.,   8.,  464.,  983., 2130., 3797., 6282., 9942., 14913., 16388.],    # 16 inch
    [0.,  14.,  791., 1674., 3628., 6465.,10698.,16931., 25396., 27908.],    # 20 inch
    [0.,  22., 1222., 2587., 5605., 9989.,16528.,26157., 39236., 43116.],    # 24 inch
    [0.,  30., 1663., 3522., 7630.,12599.,20036.,30482., 46899., 58696.],    # 28 inch
    [0.,  45., 2387., 4791., 8736.,13788.,20613.,31395., 48117., 68250.],    # 32 inch
    [0.,  60., 3021., 6063.,11055.,17499.,26086.,39731., 60895., 86375.],    # 36 inch
    [0.,  84., 4183., 8395.,15307.,24159.,36166.,55084., 84425.,119750.],    # 40 inch
])

# Maximum CV values for each concentric diameter (used for Linear Butterfly, case 3)
_CV_MAX_CONCENTRIC = np.array([
    135., 302., 600., 1022., 1579., 3136., 8250., 16388.,
    27908., 43116., 58696., 68250., 86375., 119750.
])

# Triple Offset Butterfly Valve diameters [inches] and positions [0-1]
_D_TRIPLE = np.array([3., 4., 5., 6., 8., 12., 16., 20., 24.])

_POS_TRIPLE = np.array([
    0., 1./9., 2./9., 3./9., 4./9., 5./9., 6./9., 7./9., 8./9., 1.
])

# CV table for Triple Offset Butterfly: shape (9 diameters, 10 positions)
# Source: https://www.valvesonline.com.au/references/flow-rates/butterfly-valves/
_CV_TRIPLE = np.array([
    [0.,   4.7,  16.6,   33.1,   52.1,   78.3,  110.,   190.,    210.,    240.],   # 3 inch NPS
    [0.,   8.4,  29.3,   58.6,   92.,   140.,   200.,   330.,    370.,    420.],   # 4 inch NPS
    [0.,  13.8,  48.4,   96.5,  160.,   230.,   330.,   550.,    610.,    700.],   # 5 inch NPS
    [0.,  20.9,  73.7,  150.,   230.,   350.,   500.,   820.,    930.,   1050.],   # 6 inch NPS
    [0.,  38.2, 140.,   270.,   420.,   640.,   900.,  1500.,   1690.,   1920.],   # 8 inch NPS
    [0.,  88.4, 310.,   620.,   980.,  1460.,  2080.,  3450.,   3890.,   4420.],   # 12 inch NPS
    [0., 150.,  530.,  1060.,  1660.,  2490.,  3540.,  5870.,   6620.,   7520.],   # 16 inch NPS
    [0., 270.,  930.,  1850.,  2900.,  4350.,  6190., 10300.,  11600.,  13200.],   # 20 inch NPS
    [0., 420., 1450.,  2890.,  4530.,  6800.,  9680., 16100.,  18100.,  20600.],   # 24 inch NPS
])


# ---------------------------------------------------------------------------
# Numerical derivative helpers
# ---------------------------------------------------------------------------

def drhodhcp(P_tank: float, h_tank: float, dh: float = 1000.0) -> float:
    """
    Numerical derivative of density with respect to enthalpy at constant pressure.

    Inputs
    ------
    P_tank : float
        Pressure [Pa]
    h_tank : float
        Specific enthalpy [J/kg]
    dh : float, optional
        Finite-difference step size in enthalpy [J/kg], default 1000.0

    Returns
    -------
    float
        d(rho)/d(h) at constant P  [kg^2 / (m^3 · J)]
    """
    # Solve with a symmetric finite difference over enthalpy
    h_low  = h_tank - dh
    h_high = h_tank + dh
    rho_low  = fp.density("water", P=P_tank / 1000.0, h=h_low  / 1000.0)
    rho_high = fp.density("water", P=P_tank / 1000.0, h=h_high / 1000.0)
    return (rho_high - rho_low) / (h_high - h_low)


def drhodpch(P_tank: float, h_tank: float, dP: float = 1000.0) -> float:
    """
    Numerical derivative of density with respect to pressure at constant enthalpy.

    Inputs
    ------
    P_tank : float
        Pressure [Pa]
    h_tank : float
        Specific enthalpy [J/kg]
    dP : float, optional
        Finite-difference step size in pressure [Pa], default 1000.0

    Returns
    -------
    float
        d(rho)/d(P) at constant h  [kg / (m^3 · Pa)]
    """
    P_low  = P_tank - dP
    P_high = P_tank + dP
    rho_low  = fp.density("water", P=P_low  / 1000.0, h=h_tank / 1000.0)
    rho_high = fp.density("water", P=P_high / 1000.0, h=h_tank / 1000.0)
    return (rho_high - rho_low) / (P_high - P_low)


def dudhcp(P_tank: float, h_tank: float, dh: float = 1000.0) -> float:
    """
    Numerical derivative of specific internal energy with respect to enthalpy
    at constant pressure.

    Inputs
    ------
    P_tank : float
        Pressure [Pa]
    h_tank : float
        Specific enthalpy [J/kg]
    dh : float, optional
        Finite-difference step size in enthalpy [J/kg], default 1000.0

    Returns
    -------
    float
        d(u)/d(h) at constant P  [-]

    Notes
    -----
    fp.internalenergy returns u in kJ/kg; multiply by 1000 to get J/kg before
    differencing (mirrors the Fortran ``u*1000.d0`` conversion).
    """
    h_low  = h_tank - dh
    h_high = h_tank + dh
    # fp.internalenergy returns kJ/kg → convert to J/kg
    u_low  = fp.internalenergy("water", P=P_tank / 1000.0, h=h_low  / 1000.0) * 1000.0
    u_high = fp.internalenergy("water", P=P_tank / 1000.0, h=h_high / 1000.0) * 1000.0
    return (u_high - u_low) / (h_high - h_low)


def dudpch(P_tank: float, h_tank: float, dP: float = 1000.0) -> float:
    """
    Numerical derivative of specific internal energy with respect to pressure
    at constant enthalpy.

    Inputs
    ------
    P_tank : float
        Pressure [Pa]
    h_tank : float
        Specific enthalpy [J/kg]
    dP : float, optional
        Finite-difference step size in pressure [Pa], default 1000.0

    Returns
    -------
    float
        d(u)/d(P) at constant h  [J / (kg · Pa)]

    Notes
    -----
    fp.internalenergy returns u in kJ/kg; multiply by 1000 to get J/kg.
    """
    P_low  = P_tank - dP
    P_high = P_tank + dP
    # fp.internalenergy returns kJ/kg → convert to J/kg
    u_low  = fp.internalenergy("water", P=P_low  / 1000.0, h=h_tank / 1000.0) * 1000.0
    u_high = fp.internalenergy("water", P=P_high / 1000.0, h=h_tank / 1000.0) * 1000.0
    return (u_high - u_low) / (P_high - P_low)


# ---------------------------------------------------------------------------
# Specific heat of water (numerical derivative near saturation boundary)
# ---------------------------------------------------------------------------

def f_cp_water(P: float, T: float) -> float:
    """
    Specific heat at constant pressure for water/steam.

    Uses a numerical enthalpy derivative with special logic to avoid
    crossing the saturation boundary (vapor dome).

    Inputs
    ------
    P : float
        Pressure [Pa]
    T : float
        Temperature [K]

    Returns
    -------
    float
        Specific heat at constant pressure  [J / (kg · K)]
    """
    P_kPa = P / 1000.0

    # Current enthalpy at (P, T) in J/kg
    h = fp.enthalpy("water", P=P_kPa, T=T) * 1000.0

    # Saturation temperature and saturated-liquid enthalpy at P
    T_sat   = fp.temperature("water", P=P_kPa, Q=0.0)
    h_sat_f = fp.enthalpy("water",    P=P_kPa, Q=0.0) * 1000.0  # J/kg

    if abs(T - T_sat) > 1.0:
        # No risk of crossing the vapor dome — simple central difference
        dT     = 1.0
        T_high = T + dT
        T_low  = T - dT
        h_high = fp.enthalpy("water", T=T_high, P=P_kPa) * 1000.0
        h_low  = fp.enthalpy("water", T=T_low,  P=P_kPa) * 1000.0
        return (h_high - h_low) / (T_high - T_low)

    else:
        # Near or inside the vapor dome — differentiate away from saturation boundary
        h_sat_g = fp.enthalpy("water", P=P_kPa, Q=1.0) * 1000.0  # J/kg

        if h <= h_sat_f:
            # Flow is subcooled but very close to vapor dome
            T_high = T_sat
            T_low  = T - 1.0
            h_high = fp.enthalpy("water", T=T_high, P=P_kPa) * 1000.0
            h_low  = fp.enthalpy("water", T=T_low,  P=P_kPa) * 1000.0
            return (h_high - h_low) / (T_high - T_low)
        else:
            # Flow is superheated but very close to the vapor dome
            T_low  = T_sat
            T_high = T + 1.0
            h_high = fp.enthalpy("water", T=T_high, P=P_kPa) * 1000.0
            h_low  = fp.enthalpy("water", T=T_low,  P=P_kPa) * 1000.0
            return (h_high - h_low) / (T_high - T_low)


# ---------------------------------------------------------------------------
# Liquid level in a horizontal cylindrical tank (with hemispherical ends)
# ---------------------------------------------------------------------------

def tank_level(
    Volume_liquid: float,
    tank_diameter: float,
    tank_length: float,
    previous_level: float,
    Vol_tolerance: float,
) -> float:
    """
    Solve for the liquid level height in a horizontal cylindrical tank
    with hemispherical end caps using a secant method with bisection fallback.

    Inputs
    ------
    Volume_liquid  : float  Liquid volume [m^3]
    tank_diameter  : float  Inner tank diameter [m]
    tank_length    : float  Total tank length (tip-to-tip) [m]
    previous_level : float  Previous timestep liquid level [m] — used as initial guess
    Vol_tolerance  : float  Convergence tolerance on volume [m^3]

    Returns
    -------
    float
        Liquid level height [m]
    """
    D_tank = tank_diameter
    R_tank = D_tank / 2.0
    h_guess = previous_level
    alpha = 0.5           # learning rate for secant step blending
    R_D = 1.0             # ratio used in hemispherical end-cap volume formula
    h_max = D_tank
    h_min = 0.0

    error_new = Vol_tolerance + 1.0   # ensure at least one iteration
    whileiterations = 0
    error_prev = 0.0
    h_guess_prev = h_guess

    while abs(error_new) > Vol_tolerance:
        whileiterations += 1

        if whileiterations > 1:
            # Compute cross-sectional area of liquid in horizontal cylinder
            Area_guess = (
                3.14 * D_tank**2.0 / 4.0
                * (0.5 - (np.arcsin(1.0 - 2.0 * h_guess / D_tank) / 3.14))
                - np.sqrt((D_tank / 2.0)**2.0 - (D_tank / 2.0 - h_guess)**2.0)
                * (D_tank / 2.0 - h_guess)
            )
            # Volume in cylindrical section + hemispherical end caps
            Vol_guess = (
                Area_guess * (tank_length - D_tank / R_D)
                + 3.14 * (h_guess**2.0 * D_tank / 2.0 - h_guess**3.0 / 2.0) / R_D
            )
            error_new = Vol_guess - Volume_liquid

            # Tighten bounds
            if error_new > 0.0:
                h_max = min(h_max, h_guess)
            else:
                h_min = max(h_min, h_guess)

            if abs(error_new) <= Vol_tolerance:
                pass  # converged, keep h_guess
            elif h_guess != h_guess_prev:
                m_line = (error_new - error_prev) / (h_guess - h_guess_prev)
                b_line = error_new - m_line * h_guess
                if m_line != 0.0:
                    h_new = -b_line / m_line
                    h_new = max(min(h_new, h_max), h_min)
                    h_guess_prev = h_guess
                    error_prev   = error_new
                    h_guess      = h_guess + (h_new - h_guess) * alpha
                else:
                    # Slope of the line is zero — step away from bound
                    h_guess_prev = h_guess
                    error_prev   = error_new
                    if error_new > 0.0:
                        h_guess = max(h_guess - 0.01, h_min)
                    else:
                        h_guess = min(h_guess + 0.01, h_max)
            else:
                h_guess_prev = h_guess
                error_prev   = error_new
                if error_new > 0.0:
                    h_guess = max(h_guess - 0.01, h_min)
                else:
                    h_guess = min(h_guess + 0.01, h_max)

            # Bisection fallback when reaching bounds
            if h_guess == h_max:
                h_guess = (h_max + h_min) / 2.0
            if h_guess == h_min:
                h_guess = (h_max + h_min) / 2.0

        else:
            # First iteration: compute volume with initial guess
            Area_guess = (
                3.14 * D_tank**2.0 / 4.0
                * (0.5 - (np.arcsin(1.0 - 2.0 * h_guess / D_tank) / 180.0))
                - np.sqrt((D_tank / 2.0)**2.0 - (D_tank / 2.0 - h_guess)**2.0)
                * (D_tank / 2.0 - h_guess)
            )
            Vol_guess = (
                Area_guess * (tank_length - D_tank / R_D)
                + 3.14 * (h_guess**2.0 * D_tank / 2.0 - h_guess**3.0 / 2.0) / R_D
            )
            error_new = Vol_guess - Volume_liquid
            error_prev = error_new
            h_guess_prev = h_guess

            if abs(error_prev) < Vol_tolerance:
                pass  # initial guess already converged
            elif error_prev > 0.0:
                # Guess area is greater than actual — decrease height
                h_guess = max(h_guess - 0.01, 0.0)
            else:
                # Guess area is smaller than actual — increase height
                h_guess = min(h_guess + 0.01, tank_diameter)

        # Safety limit to prevent infinite loop
        if whileiterations > 100:
            error_new = Vol_tolerance  # force exit

    return h_guess_prev


# ---------------------------------------------------------------------------
# Power-block Cv data (distinct larger tables from the simple Valve Cv tables)
# ---------------------------------------------------------------------------

def PB_CV_data(Valve_Type: int, D_in: float, Valve_position: float) -> float:
    """
    Flow coefficient (Cv) for power-block butterfly valves.

    IMPORTANT: These lookup tables are DIFFERENT and LARGER than the CV_data
    tables used in the basic Valve component.

    Inputs
    ------
    Valve_Type     : int    Valve type selector:
                              1 = Concentric Butterfly (14 diameters, 2–40 in)
                              2 = Triple Offset Butterfly (9 diameters, 3–24 in)
                              3 = Linear Opening Butterfly (uses concentric CV_max)
    D_in           : float  Valve inlet diameter [m]
    Valve_position : float  Valve position [0–1], where 1 = fully open

    Returns
    -------
    float
        Flow coefficient Cv [-]
    """
    # Return a tiny non-zero Cv so the hydraulic model will not crash
    if Valve_position == 0.0:
        return 1.0e-9

    # Convert diameter from metres to inches
    D = D_in * 39.3701

    if Valve_Type == 1:
        # ---- Concentric Butterfly ----
        D_arr   = _D_CONCENTRIC
        Pos_arr = _POS_CONCENTRIC
        CV_arr  = _CV_CONCENTRIC

        # --- Diameter index ---
        match_D = False
        ind_D   = -1
        found   = False
        ind_D_low = ind_D_high = 0
        for n in range(len(D_arr)):
            if abs(D - D_arr[n]) < 0.05:
                match_D = True
                ind_D = n
            elif D_arr[n] > D and not found:
                ind_D_low  = n - 1
                ind_D_high = n
                found = True

        # --- Position index ---
        match_Pos = False
        ind_Pos   = -1
        found_pos = False
        ind_Pos_low = ind_Pos_high = 0
        for n in range(len(Pos_arr)):
            if Pos_arr[n] == Valve_position:
                match_Pos = True
                ind_Pos = n
            elif Pos_arr[n] > Valve_position and not found_pos:
                ind_Pos_low  = n - 1
                ind_Pos_high = n
                found_pos = True

        # --- Bilinear interpolation ---
        if match_Pos and match_D:
            C_v = CV_arr[ind_D, ind_Pos]
        elif match_D:
            # Interpolate in position only
            x0 = Pos_arr[ind_Pos_low]
            x1 = Pos_arr[ind_Pos_high]
            y0 = CV_arr[ind_D, ind_Pos_low]
            y1 = CV_arr[ind_D, ind_Pos_high]
            C_v = y0 + (Valve_position - x0) * (y1 - y0) / (x1 - x0)
        elif match_Pos:
            # Interpolate in diameter only
            x0 = D_arr[ind_D_low]
            x1 = D_arr[ind_D_high]
            y0 = CV_arr[ind_D_low,  ind_Pos]
            y1 = CV_arr[ind_D_high, ind_Pos]
            C_v = y0 + (D - x0) * (y1 - y0) / (x1 - x0)
        else:
            # Full bilinear interpolation
            x0  = Pos_arr[ind_Pos_low]
            x1  = Pos_arr[ind_Pos_high]
            yL0 = CV_arr[ind_D_low,  ind_Pos_low]
            yL1 = CV_arr[ind_D_low,  ind_Pos_high]
            yH0 = CV_arr[ind_D_high, ind_Pos_low]
            yH1 = CV_arr[ind_D_high, ind_Pos_high]
            y0  = yL0 + (Valve_position - x0) * (yL1 - yL0) / (x1 - x0)
            y1  = yH0 + (Valve_position - x0) * (yH1 - yH0) / (x1 - x0)
            x0  = D_arr[ind_D_low]
            x1  = D_arr[ind_D_high]
            C_v = y0 + (D - x0) * (y1 - y0) / (x1 - x0)

    elif Valve_Type == 2:
        # ---- Triple Offset Butterfly ----
        D_arr   = _D_TRIPLE
        Pos_arr = _POS_TRIPLE
        CV_arr  = _CV_TRIPLE

        # --- Diameter index ---
        match_D = False
        ind_D   = -1
        found   = False
        ind_D_low = ind_D_high = 0
        for n in range(len(D_arr)):
            if abs(D - D_arr[n]) < 0.05:
                match_D = True
                ind_D = n
            elif D_arr[n] > D and not found:
                ind_D_low  = n - 1
                ind_D_high = n
                found = True

        # --- Position index ---
        match_Pos = False
        ind_Pos   = -1
        found_pos = False
        ind_Pos_low = ind_Pos_high = 0
        for n in range(len(Pos_arr)):
            if Pos_arr[n] == Valve_position:
                match_Pos = True
                ind_Pos = n
            elif Pos_arr[n] > Valve_position and not found_pos:
                ind_Pos_low  = n - 1
                ind_Pos_high = n
                found_pos = True

        # --- Bilinear interpolation ---
        if match_Pos and match_D:
            C_v = CV_arr[ind_D, ind_Pos]
        elif match_D:
            x0 = Pos_arr[ind_Pos_low]
            x1 = Pos_arr[ind_Pos_high]
            y0 = CV_arr[ind_D, ind_Pos_low]
            y1 = CV_arr[ind_D, ind_Pos_high]
            C_v = y0 + (Valve_position - x0) * (y1 - y0) / (x1 - x0)
        elif match_Pos:
            x0 = D_arr[ind_D_low]
            x1 = D_arr[ind_D_high]
            y0 = CV_arr[ind_D_low,  ind_Pos]
            y1 = CV_arr[ind_D_high, ind_Pos]
            C_v = y0 + (D - x0) * (y1 - y0) / (x1 - x0)
        else:
            x0  = Pos_arr[ind_Pos_low]
            x1  = Pos_arr[ind_Pos_high]
            yL0 = CV_arr[ind_D_low,  ind_Pos_low]
            yL1 = CV_arr[ind_D_low,  ind_Pos_high]
            yH0 = CV_arr[ind_D_high, ind_Pos_low]
            yH1 = CV_arr[ind_D_high, ind_Pos_high]
            y0  = yL0 + (Valve_position - x0) * (yL1 - yL0) / (x1 - x0)
            y1  = yH0 + (Valve_position - x0) * (yH1 - yH0) / (x1 - x0)
            x0  = D_arr[ind_D_low]
            x1  = D_arr[ind_D_high]
            C_v = y0 + (D - x0) * (y1 - y0) / (x1 - x0)

    elif Valve_Type == 3:
        # ---- Linear Opening Butterfly ----
        # Uses maximum Cv from the concentric table scaled linearly by valve position
        D_arr    = _D_CONCENTRIC
        CV_max   = _CV_MAX_CONCENTRIC

        # --- Diameter index ---
        match_D = False
        ind_D   = -1
        found   = False
        ind_D_low = ind_D_high = 0
        for n in range(len(D_arr)):
            if abs(D - D_arr[n]) < 0.05:
                match_D = True
                ind_D = n
            elif D_arr[n] > D and not found:
                ind_D_low  = n - 1
                ind_D_high = n
                found = True

        if match_D:
            # Diameter matched — no interpolation needed
            CV_max_d = CV_max[ind_D]
        else:
            if ind_D_low > 0 and found:
                # Valve diameter is within lookup table range — interpolate
                x0 = D_arr[ind_D_low]
                x1 = D_arr[ind_D_high]
                y0 = CV_max[ind_D_low]
                y1 = CV_max[ind_D_high]
                CV_max_d = y0 + (D - x0) * (y1 - y0) / (x1 - x0)
            elif ind_D_low == 0:
                # Valve diameter is smaller than smallest entry — extrapolate
                CV_max_d = CV_max[0] * D / D_arr[0]
            else:
                # Valve diameter is larger than last entry — extrapolate
                CV_max_d = CV_max[13] * D / D_arr[13]

        # Scale CV_max by valve position (linear relationship)
        CV_max_d = CV_max_d / 2.0
        C_v = CV_max_d * Valve_position

    else:
        raise ValueError(f"Unsupported Valve_Type={Valve_Type}. Must be 1, 2, or 3.")

    return C_v


# ---------------------------------------------------------------------------
# Valve position ramp function
# ---------------------------------------------------------------------------

def VP_new(
    VP_Current: float,
    VP_Request: float,
    valve_speed: float,
    ts: float,
) -> float:
    """
    Ramp valve position toward a requested value at a given speed.

    Inputs
    ------
    VP_Current  : float  Current valve position [0–1 fraction open]
    VP_Request  : float  Requested valve position [0–1 fraction open]
    valve_speed : float  Maximum valve angular speed [deg/s]
    ts          : float  Simulation timestep [s]

    Returns
    -------
    float
        New valve position [0–1]
    """
    # Cap valve speed to physical limits
    if valve_speed > 90.0:
        valve_speed = 90.0
    if valve_speed <= 0.0:
        # Valve position cannot change — set to 1 deg/s minimum
        valve_speed = 1.0

    # Convert fractional openings to angular degrees (0–90°)
    VP_current_d = VP_Current * 90.0
    VP_request_d = VP_Request * 90.0

    if VP_current_d != VP_request_d:
        if VP_request_d > VP_current_d:
            # Open valve
            return min(VP_current_d + valve_speed * ts, VP_request_d) / 90.0
        else:
            # Close valve
            return max(VP_current_d - valve_speed * ts, VP_request_d) / 90.0
    else:
        return VP_Current


# ---------------------------------------------------------------------------
# General secant-method helper
# ---------------------------------------------------------------------------

def error_function(
    error_new: float,
    error_prev: float,
    h_new: float,
    h_prev: float,
    h_min: float,
    h_max: float,
    iterations: float,
    tol: float,
) -> float:
    """
    General secant method helper that returns the next guess for an iterative solver.

    Uses a secant step when two data points are available and falls back to
    a midpoint bisection when the next guess would land on a bound.

    Inputs
    ------
    error_new   : float  Current residual
    error_prev  : float  Previous residual
    h_new       : float  Current guess value
    h_prev      : float  Previous guess value
    h_min       : float  Lower bound for the guess
    h_max       : float  Upper bound for the guess
    iterations  : float  Iteration counter (Fortran stores as double; >1 means ≥2 pts)
    tol         : float  Convergence tolerance

    Returns
    -------
    float
        Next guess value
    """
    if abs(error_new) <= tol:
        return h_new

    if iterations > 1.0:
        if h_new == h_prev:
            # Degenerate case — step away from current value
            if error_new > 0.0:
                result = min(h_new + 100000.0, h_max)
            else:
                result = max(h_new - 100000.0, h_min)
        else:
            m = (error_new - error_prev) / (h_new - h_prev)
            b = error_new - m * h_new
            if m != 0.0:
                result = max(min(-b / m, h_max), h_min)
            else:
                if error_new > 0.0:
                    result = min(h_new + 100000.0, h_max)
                else:
                    result = max(h_new - 100000.0, h_min)
    else:
        # First iteration — large step in appropriate direction
        if error_new > 0.0:
            result = min(h_new + 1000000.0, h_max)
        else:
            result = max(h_new - 1000000.0, h_min)

    # Bisection fallback — keeps from guessing min and max again
    if result == h_max:
        result = (h_max + h_min) / 2.0  # helps converge faster
    elif result == h_min:
        result = (h_max + h_min) / 2.0  # helps converge faster

    return result


# ---------------------------------------------------------------------------
# Stodola's ellipse stage model
# ---------------------------------------------------------------------------

def StodolaStage(
    P_in_d: float,
    P_out_d: float,
    P_in: float,
    h_in_d: float,
    h_in: float,
    m_dot_in_d: float,
    m_dot_in: float,
) -> float:
    """
    Solve for the outlet pressure of a turbine stage using Stodola's ellipse.

    Inputs
    ------
    P_in_d    : float  Design inlet pressure [Pa]
    P_out_d   : float  Design outlet pressure [Pa]
    P_in      : float  Actual inlet pressure [Pa]
    h_in_d    : float  Design inlet specific enthalpy [J/kg]
    h_in      : float  Actual inlet specific enthalpy [J/kg]
    m_dot_in_d: float  Design mass flow rate [kg/s]
    m_dot_in  : float  Actual mass flow rate [kg/s]

    Returns
    -------
    float
        Outlet pressure of the stage [Pa]

    Notes
    -----
    fp.specific_volume returns v in m^3/kg.
    """
    # Specific volumes at design and actual conditions
    v_in_d = fp.specific_volume("water", P=P_in_d / 1000.0, h=h_in_d / 1000.0)
    v_in   = fp.specific_volume("water", P=P_in   / 1000.0, h=h_in   / 1000.0)

    PR_d    = P_out_d / P_in_d
    phi_D   = m_dot_in_d / np.sqrt(P_in_d / v_in_d)
    phi_max = np.sqrt(0.999999 / (1.0 - PR_d**2.0)) * phi_D
    phi_a   = min(m_dot_in / np.sqrt(P_in / v_in), phi_max)
    PR_a    = np.sqrt(1.0 - (phi_a / phi_D)**2.0 * (1.0 - PR_d**2.0))
    return PR_a * P_in


# ---------------------------------------------------------------------------
# LPT stage enthalpy solver (expansion line + steam tables intersection)
# ---------------------------------------------------------------------------

def h_lpt_stage(
    h_guess: float,
    h_min: float,
    h_max: float,
    P_stage: float,
    DELTA_S: float,
    A0: float,
    A1: float,
    A2: float,
    tol: float,
) -> float:
    """
    Find actual specific enthalpy at a LPT stage exit by matching the expansion
    line entropy to the (P, h) entropy from steam tables using a secant solver.

    Inputs (all SI)
    ------
    h_guess  : float  Initial guess enthalpy [J/kg]
    h_min    : float  Minimum bound (exhaust enthalpy) [J/kg]
    h_max    : float  Maximum bound (inlet enthalpy) [J/kg]
    P_stage  : float  Stage outlet pressure from Stodola's ellipse [Pa]
    DELTA_S  : float  Offset applied to expansion line entropy [J/(kg·K)]
    A0, A1, A2 : float  Expansion line polynomial coefficients
    tol      : float  Convergence tolerance on entropy residual [J/(kg·K)]

    Returns
    -------
    float
        Solved stage exit enthalpy [J/kg]
    """
    if h_guess < h_min:
        h_guess = (h_min + h_max) / 2.0  # using enthalpy between max and min
    if h_guess > h_max:
        h_guess = (h_min + h_max) / 2.0  # using enthalpy between max and min

    error = tol + 10.0    # making error greater than tolerance for the first iteration
    whileiterations = 0.0  # reset to zero
    alpha = 0.5            # learning rate used for while iterations
    maxiterations = 50.0
    h_prev = h_guess
    error_prev = 0.0
    h_new = h_guess        # initialise h_new before the loop

    while abs(error) > tol:
        whileiterations += 1.0
        s_elep = DELTA_S + A0 + A1 * h_guess + A2 * h_guess ** 2.0
        # Converting from kJ/kg-K to J/kg-K
        s_fit = fp.entropy("water", P=P_stage / 1000.0, h=h_guess / 1000.0) * 1000.0
        error = s_elep - s_fit

        if abs(error) < tol:
            break

        elif whileiterations > 1.0:
            if h_guess != h_prev:  # check to make sure that the slope is not undefined
                m = (error - error_prev) / (h_guess - h_prev)  # slope of error line
                b = error - m * h_guess                          # y-intercept of error line
                if m != 0.0:   # make sure slope does not equal zero
                    h_new = -b / m   # new guess value
                    h_new = min(h_max, h_new)  # check that it's lower than the maximum guess
                    h_new = max(h_min, h_new)  # check that it's greater than the minimum guess
                else:
                    if error > 0.0:  # left of convergence point, decrease enthalpy guess
                        # TODO-NEEDS CONVERSION REVIEW: Fortran uses max(h_guess+10, h_max) here;
                        # this may step outside bounds — ported faithfully
                        h_new = max(h_guess + 10.0, h_max)
                    else:  # right of convergence point, increase enthalpy guess
                        h_new = min(h_guess - 10.0, h_min)
            else:
                if error > 0.0:  # left of convergence point, decrease enthalpy guess
                    # TODO-NEEDS CONVERSION REVIEW: see note above
                    h_new = max(h_guess + 10.0, h_max)
                else:
                    h_new = min(h_guess - 10.0, h_min)

        else:  # first iteration
            if error > 0.0:  # left of convergence point, decrease enthalpy guess
                h_new = (h_guess + h_max) / 2.0
            else:  # right of the convergence point, increase enthalpy guess
                h_new = (h_guess + h_min) / 2.0

        error_prev = error  # set new error to previous error
        h_prev = h_guess    # set current guess to previous guess
        h_guess = h_guess + (h_new - h_guess) * alpha

        if maxiterations == whileiterations:
            break  # force exit with current alpha-blended h_guess (mirrors Fortran do-while exit)

    return h_guess


# ---------------------------------------------------------------------------
# Spencer-Cotton-Cannon turbine efficiency models
# ---------------------------------------------------------------------------

def eta_SCC_hpt(
    m_dot_in: float,
    m_dot_rated: float,
    N_parallel: float,
    PD: float,
    P_in: float,
    v_in: float,
    v_design_in: float,
    Design_EP: float,
    N_CV: float,
    N_row: float,
) -> float:
    """
    Spencer-Cotton-Cannon efficiency model for the High Pressure Turbine.

    Inputs (SI units)
    -----------------
    m_dot_in     : float  Actual mass flow entering the turbine [kg/s]
    m_dot_rated  : float  Design mass flow for the turbine [kg/s]
    N_parallel   : float  Number of parallel sections in the turbine [-]
    PD           : float  Pitch diameter of the governing stage [m]
    P_in         : float  Inlet pressure (actual conditions) [Pa]
    v_in         : float  Inlet specific volume (actual conditions) [m^3/kg]
    v_design_in  : float  Inlet specific volume (design conditions) [m^3/kg]
    Design_EP    : float  Design exhaust pressure [Pa]
    N_CV         : float  Number of control valves ahead of the turbine [-]
    N_row        : float  Number of governing stages (1 or 2) [-]

    Returns
    -------
    float
        Isentropic efficiency [-]
    """
    # AUTO UNITS CONVERSION IMPLEMENTED: SI → English units for SCC correlation
    m_dot_in_eng    = m_dot_in     * 7936.64   # kg/s → lb/hr
    m_dot_rated_eng = m_dot_rated  * 7936.64   # kg/s → lb/hr
    PD_eng          = PD           * 39.3701   # m → in
    P_in_eng        = P_in         / 6894.76   # Pa → psia
    Design_EP_eng   = Design_EP    / 6894.76   # Pa → psia
    v_in_eng        = v_in         * 16.0185   # m^3/kg → ft^3/lb
    v_design_in_eng = v_design_in  * 16.0185   # m^3/kg → ft^3/lb
    Vol_dot_rated_eng = m_dot_rated_eng * v_design_in_eng

    if N_row == 1.0:
        # Find Throttle Flow Ratio (TFR)
        TFR = m_dot_in_eng / m_dot_rated_eng

        # Start from Base Efficiency (From Table 1)
        eta_base = 87.00

        # Efficiency Correction for Volume Flow - Poorer (From Table 1)
        delta_eta = 1005200.0 * N_parallel / Vol_dot_rated_eng / 100.0
        eta_base = eta_base - delta_eta * eta_base

        # Efficiency Correction for Governing Stage (From Fig 7)
        delta_eta = (-0.115 * PD_eng + 4.37) / 100.0  # Found in Fig 7. of SCC paper
        eta_base = eta_base + delta_eta * eta_base

        # Correction for Pressure Ratio (From Fig 6)
        x = Design_EP_eng / P_in_eng
        y = math.log(Vol_dot_rated_eng)
        delta_eta = (11.151 - 63.0 * x - 0.50091 * y + 2.83 * x * y) / 100.0
        eta_base = eta_base + delta_eta * eta_base

        # Correction for Governing Stage at Partial Load (From Fig 8)
        delta_eta = (-21.8085 + 21.8085 * TFR + 0.573908 * PD_eng - 0.573908 * TFR * PD_eng) / 100.0
        eta_base = eta_base + delta_eta * eta_base

        # Correction for Partial Load Operation (From Fig 9)
        x = TFR
        y = math.log(P_in_eng / Design_EP_eng)
        delta_eta = (
            -60.75 + 66.85 * x + 29.75 * x**2.0 - 35.85 * x**3.0
            + 17.50 * y - 20.02 * y * x - 0.525 * y * x**2.0 + 3.045 * y * x**3.0
        ) / 100.0
        eta_base = eta_base + delta_eta * eta_base

        # Correction for Mean of Loops (From Fig 12)
        x = TFR
        y = N_CV
        delta_eta = (-5.4 + 4.395 * x + 0.45 * N_CV - 0.36625 * x * y) / 100.0
        eta_base = eta_base + delta_eta * eta_base
        # end of 1 row governing stage efficiency calculations

    else:  # 2 row governing stage
        # Find Throttle Flow Ratio (TFR)
        m_dot_in_eng    = m_dot_in   * 7936.64  # converting from kg/s to lb/hr
        m_dot_rated_eng = m_dot_rated * 7936.64  # converting from kg/s to lb/hr
        TFR = m_dot_in_eng / m_dot_rated_eng

        # Start from Base Efficiency (From Table 1)
        eta_base = 84.00

        # Efficiency Correction for Volume Flow - Poorer (From Table 1)
        delta_eta = (1350000.0 * N_parallel / (m_dot_in_eng * v_in_eng)) / 100.0
        eta_base = eta_base - delta_eta * eta_base

        # Efficiency Correction for pressure ratio (From Fig 10)
        x = Design_EP_eng / P_in_eng
        y = math.log(Vol_dot_rated_eng)
        delta_eta = (25.665 - 145.0 * x - 1.33281 * y + 7.53 * x * y) / 100.0
        eta_base = eta_base + delta_eta * eta_base

        # Efficiency Correction for partial load (From Fig 11)
        x = 1.0 - TFR
        y = P_in_eng / Design_EP_eng
        delta_eta = (
            42.676909 * x - 89.391147 * x**2.0 + 9.0376638 * x**3.0
            - 26.221836 * x * y + 25.549385 * y * x**2.0 + 8.8283868 * y * x**3.0
            + 4.0479550 * y**2.0 * x - 1.4725197 * x**2.0 * y**2.0
            - 4.0183332 * y**2.0 * x**3.0 - 0.14502211 * y**3.0 * x
            - 0.18580363 * y**3.0 * x**2.0 + 0.42657518 * x**3.0 * y**3.0
        ) / 100.0
        eta_base = eta_base + delta_eta * eta_base

        # Efficiency Correction for mean of loops (Fig 12)
        x = TFR
        y = N_CV
        delta_eta = (-5.4 + 4.395 * x + 0.45 * N_CV - 0.36625 * x * y) / 100.0
        eta_base = eta_base + delta_eta * eta_base
        # end of 2 row governing stage efficiency calculations

    return eta_base / 100.0


# Alias for call sites that use the Fortran convention (case-insensitive name)
eta_SCC_HPT = eta_SCC_hpt


def eta_SCC_lpt(
    m_dot_in: float,
    m_dot_rated: float,
    N_parallel: float,
    P_in: float,
    T_in: float,
    v_in: float,
    v_design_in: float,
    Design_EP: float,
) -> float:
    """
    Spencer-Cotton-Cannon efficiency model for the Low Pressure Turbine.

    Inputs (SI units)
    -----------------
    m_dot_in     : float  Actual mass flow entering the turbine [kg/s]
    m_dot_rated  : float  Design mass flow for the turbine [kg/s]
    N_parallel   : float  Number of parallel sections in the turbine [-]
    P_in         : float  Inlet pressure (actual conditions) [Pa]
    T_in         : float  Inlet temperature (actual conditions) [K]
    v_in         : float  Inlet specific volume (actual conditions) [m^3/kg]
    v_design_in  : float  Inlet specific volume (design conditions) [m^3/kg]
    Design_EP    : float  Design exhaust pressure [Pa]

    Returns
    -------
    float
        Isentropic efficiency [-]
    """
    # AUTO UNITS CONVERSION IMPLEMENTED: SI → English units for SCC correlation
    m_dot_in_eng    = m_dot_in   * 7936.64   # Converting from kg/s to lb/hr
    m_dot_rated_eng = m_dot_rated * 7936.64  # Converting from kg/s to lb/hr
    P_in_eng        = P_in       / 6894.76   # Converting from Pa to psia
    v_in_eng        = v_in       * 16.0185   # Converting from m^3/kg to ft^3/lb
    v_design_in_eng = v_design_in * 16.0185  # Converting from m^3/kg to ft^3/lb
    T_in_eng        = (T_in - 273.15) * 9.0 / 5.0 + 32.0  # Converting from Kelvin to Fahrenheit
    Vol_dot_rated_eng = m_dot_rated_eng * v_design_in_eng  # Finding rated volumetric flow rate

    # Start with SCC base efficiency
    eta_base = 91.93  # reheat section, 3600 rpm without governing stage baseline efficiency

    # Efficiency correction for governing stage
    delta_eta = 1270000 * N_parallel / Vol_dot_rated_eng
    delta_eta = delta_eta / 100.0
    eta_base = eta_base - delta_eta * eta_base

    # Efficiency correction for initial pressure and temperature (FIG 14)
    # Converting from kJ/kg-K to J/kg-K
    s_in = fp.entropy("water", P=P_in / 1000.0, T=T_in) * 1000.0
    # Converting from kJ/kg to J/kg
    h_in = fp.enthalpy("water", P=P_in / 1000.0, T=T_in) * 1000.0
    s_in_eng = 0.0002388459 * s_in   # Convert from J/kg-K to BTU/lbm-R
    h_in_eng = h_in * 0.0004299226   # Convert from J/kg to BTU/lbm

    Fig14Values = np.array([
        [ 28.232252,    -92.390491,   -625.79590,   207.2301,    70.251642,  -22.516388  ],
        [ -0.047796308,   1.2844571,    0.38556961,  -0.039652999,-0.27180357,  0.064869467],
        [ -0.69791427e-3,-0.17037268e-2, 0.86563845e-3,-0.59510660e-3, 0.39705804e-3,-0.73533255e-4],
        [  0.12050837e-5,  0.26826382e-6,-0.67887771e-6, 0.52886157e-6,-0.24106229e-6,  0.37881801e-7],
        [ -0.50719109e-9,  0.26393497e-9, 0.38021911e-10,-0.10149993e-9, 0.47757232e-10,-0.70989561e-11],
    ])

    delta_eta = 0.0  # Reset DELTA_ETA to 0
    x = math.log10(P_in_eng)

    if s_in_eng > 2.0041:
        y = min(h_in_eng, 1154.0 + 80.0 * x + 88.0 * x**2.0)
    else:
        y = h_in_eng

    for j in range(5):
        for i in range(6):
            delta_eta = delta_eta + Fig14Values[j, i] * (x**i) * (y**j)

    delta_eta = delta_eta / 100.0  # Convert to a percentage
    eta_base = eta_base + delta_eta * eta_base

    return eta_base / 100.0  # Converting from percentage to decimal


# Alias for call sites that use the Fortran convention (case-insensitive name)
eta_SCC_LPT = eta_SCC_lpt


# ---------------------------------------------------------------------------
# Viscosity of steam (empirical exponential correlation)
# ---------------------------------------------------------------------------

def viscosity_steam(T: float) -> float:
    """
    Dynamic viscosity of steam using an empirical exponential correlation.

    Inputs
    ------
    T : float  Temperature [K]

    Returns
    -------
    float
        Dynamic viscosity [Pa·s]
    """
    A = 0.00000000001856
    B = 4209.0
    C = 0.04527
    D = -0.00003376
    return A * math.exp(B / T + C * T + D * T**2.0)


# ---------------------------------------------------------------------------
# Specific heat at constant volume for water/steam
# ---------------------------------------------------------------------------

def f_cv_water(P: float, T: float) -> float:
    """
    Specific heat at constant volume for water/steam.
    Uses a numerical enthalpy derivative with special logic to avoid
    crossing the saturation boundary (vapor dome).

    Inputs
    ------
    P : float  Pressure [Pa]
    T : float  Temperature [K]

    Returns
    -------
    float
        Specific heat at constant volume  [J/(kg·K)]

    Notes
    -----
    The Fortran uses FIT_TD (T, density) to differentiate u at constant density.
    Here we approximate using u = h - P/rho evaluated at (T ± 1K, P), which
    is exact for ideal gases and a close approximation for real steam above the
    saturation boundary.
    """
    P_kPa = P / 1000.0
    h      = fp.enthalpy("water", P=P_kPa, T=T) * 1000.0  # J/kg
    rho    = fp.density("water",  P=P_kPa, T=T)
    T_sat  = fp.temperature("water", P=P_kPa, Q=0.0)
    h_sat_f = fp.enthalpy("water", P=P_kPa, Q=0.0) * 1000.0  # J/kg

    # Helper: internal energy approximation at constant P via u = h - P/rho
    def u_at_T(T_eval: float) -> float:
        h_e   = fp.enthalpy("water", T=T_eval, P=P_kPa) * 1000.0  # J/kg
        rho_e = fp.density("water",  T=T_eval, P=P_kPa)
        return h_e - P / rho_e  # J/kg

    if abs(T - T_sat) > 1.0:  # No worries about running into vapor dome
        dT = 1.0
        return (u_at_T(T + dT) - u_at_T(T - dT)) / (2.0 * dT)
    else:
        h_sat_g = fp.enthalpy("water", P=P_kPa, Q=1.0) * 1000.0  # J/kg
        if h <= h_sat_f:   # Flow is subcooled but very close to vapor dome
            T_high = T_sat
            T_low  = T - 1.0
        else:              # Flow is superheated but very close to the vapor dome
            T_low  = T_sat
            T_high = T + 1.0
        return (u_at_T(T_high) - u_at_T(T_low)) / (T_high - T_low)


# ---------------------------------------------------------------------------
# Darcy–Weisbach friction factor (implicit Colebrook equation)
# ---------------------------------------------------------------------------

def FricFactor_IC(Rough: float, Reynold: float, guess: float) -> float:
    """
    Solve the implicit Colebrook equation for the Darcy friction factor using
    the Newton–Raphson / secant method.

    Inputs
    ------
    Rough   : float  Relative roughness (ε/D) [-]
    Reynold : float  Reynolds number [-]
    guess   : float  Initial guess for the friction factor [-]

    Returns
    -------
    float
        Darcy friction factor [-]
    """
    # Rough is relative roughness [--]
    if Reynold < 2750.0:
        return 64.0 / max(Reynold, 1.0)

    Acc = 0.00001
    X = 1.0 / math.sqrt(max(guess, 1.0e-10))
    TestOld = X + 2.0 * math.log10(Rough / 3.7 + 2.51 * X / Reynold)
    Xold = X
    X = X * 0.7
    NumTries = 0

    while True:
        NumTries += 1
        Test = X + 2.0 * math.log10(Rough / 3.7 + 2.51 * X / Reynold)
        if abs(Test - TestOld) <= Acc:
            return 1.0 / (X * X)
        if NumTries > 20:
            # Could not find friction factor solution — return best current estimate
            return 1.0 / (X * X)
        Slope = (Test - TestOld) / (X - Xold)
        Xold = X
        TestOld = Test
        X = max((Slope * X - Test) / Slope, 1.0e-5)


# ---------------------------------------------------------------------------
# Pipe convection coefficient (dynamic pipe, single-phase and two-phase)
# ---------------------------------------------------------------------------

def convection_dynamicpipe(
    P_steam: float,
    h_steam: float,
    D_pipe: float,
    m_dot: float,
    ff_guess: float,
    T_metal: float,
) -> float:
    """
    Convection heat transfer coefficient between steam and pipe wall.

    Uses the Gnielinski correlation for turbulent single-phase superheated steam,
    and the Shah (2016) / Chato (1962) correlations for two-phase condensing conditions.

    Inputs
    ------
    P_steam  : float  Steam pressure [Pa]
    h_steam  : float  Steam specific enthalpy [J/kg]
    D_pipe   : float  Pipe inner diameter [m]
    m_dot    : float  Mass flow rate [kg/s]
    ff_guess : float  Previous iteration friction factor guess [-]
    T_metal  : float  Pipe wall (metal) temperature [K]

    Returns
    -------
    float
        Convection coefficient between steam and pipe wall [W/(m^2·K)]
    """
    P_kPa = P_steam / 1000.0
    T_steam   = fp.temperature("water", P=P_kPa, h=h_steam / 1000.0)
    rho_steam = fp.density("water",     P=P_kPa, h=h_steam / 1000.0)
    x         = fp.quality("water",     P=P_kPa, h=h_steam / 1000.0)
    mu_steam  = fp.viscosity("water",   P=P_kPa, h=h_steam / 1000.0) / 1000000.0   # converting from microPa-s to Pa-s
    k_steam   = fp.conductivity("water",P=P_kPa, h=h_steam / 1000.0)

    T_sat     = fp.temperature("water", P=P_kPa, Q=1.0)
    h_sat_g   = fp.enthalpy("water",    P=P_kPa, Q=1.0) * 1000.0
    rho_sat_g = fp.density("water",     P=P_kPa, Q=1.0)
    mu_G      = fp.viscosity("water",   P=P_kPa, Q=1.0) / 1000000.0   # converting from microPa-s to Pa-s
    k_G       = fp.conductivity("water",P=P_kPa, Q=1.0)

    h_sat_f   = fp.enthalpy("water",    P=P_kPa, Q=0.0) * 1000.0
    rho_sat_f = fp.density("water",     P=P_kPa, Q=0.0)
    mu_L      = fp.viscosity("water",   P=P_kPa, Q=0.0) / 1000000.0   # converting from microPa-s to Pa-s
    k_L       = fp.conductivity("water",P=P_kPa, Q=0.0)

    Area   = 3.14 / 4.0 * D_pipe**2.0
    cp_L   = 4200.0
    P_crit = 22064000.0  # Pa

    if h_steam >= h_sat_g:  # steam in pipe is superheated
        vel = m_dot / (rho_steam * Area)
        # use single phase vapor coefficient (Gnielinski Correlation)
        Re = rho_steam * vel * D_pipe / mu_steam
        if Re > 2300.0:  # The flow is turbulent
            cp = f_cp_water(P=P_steam, T=T_steam)
            Pr = mu_steam * cp / k_steam
            ff = FricFactor_IC(0.0, Re, ff_guess)
            Nu = ((ff / 8.0) * (Re - 1000.0) * Pr) / (1.0 + 12.7 * (ff / 8.0)**0.5 * (Pr**(2.0 / 3.0) - 1.0))
            h_bar = Nu * k_steam / D_pipe
        else:
            Nu    = 3.66  # Fully-developed Nusselt Number for laminar flow with a uniform wall temperature
            h_bar = Nu * k_steam / D_pipe
    else:  # steam in pipe is saturated, need two phase coefficient
        vel = m_dot / (rho_sat_g * Area)
        Re  = rho_sat_g * vel * D_pipe / mu_G
        if Re >= 35000.0:
            Z      = (1.0 / x - 1.0)**0.8 * (P_steam / P_crit)**0.4
            G_tot  = m_dot / Area
            J_g    = x * G_tot / math.sqrt(9.81 * D_pipe * rho_sat_g * (rho_sat_f - rho_sat_g))
            J_g_I  = 0.98 * (Z + 0.263)**(-0.62)
            if J_g >= J_g_I:  # First Heat Transfer Region
                Re_L  = rho_sat_f * vel * D_pipe / mu_L
                Pr_L  = mu_L * cp_L / k_L
                h_L   = 0.23 * Re_L**0.8 * Pr_L**0.4 * k_L / D_pipe
                h_bar = h_L * (1.0 + 3.8 / (Z**0.95)) * (mu_L / (14.0 * mu_G))**(0.0058 + 0.557 * P_steam / P_crit)
            else:
                J_g_III = 0.95 * (1.254 + 2.27 * Z**1.249)**(-1.0)
                if J_g <= J_g_III:  # Third Heat Transfer Region
                    Re_L  = rho_sat_f * vel * D_pipe / mu_L
                    h_NU  = 1.32 * Re_L**(1.0 / 3.0) * ((rho_sat_f * (rho_sat_f - rho_sat_g) * 9.81 * k_L**3.0) / mu_L**2.0)**(1.0 / 3.0)
                    h_bar = h_NU
                else:
                    Re_L       = rho_sat_f * vel * D_pipe / mu_L
                    Pr_L       = mu_L * cp_L / k_L
                    mu_g_steam = viscosity_steam(T_steam)
                    h_L   = 0.23 * Re_L**0.8 * Pr_L**0.4 * k_L / D_pipe
                    h_I   = h_L * (1.0 + 3.8 / (Z**0.95)) * (mu_L / (14.0 * mu_G))**(0.0058 + 0.557 * P_steam / P_crit)
                    h_NU  = 1.32 * Re_L**(1.0 / 3.0) * ((rho_sat_f * (rho_sat_f - rho_sat_g) * 9.81 * k_L**3.0) / mu_L**2.0)**(1.0 / 3.0)
                    h_bar = h_I + h_NU
        else:
            # Reynolds number too low for Shah 2016 Correlation
            # Using Chate 1962 Correlation for film condensation in a horizontal pipe
            if abs(T_sat - T_metal) > 0.01:
                h_fg  = (h_sat_g - h_sat_f) + 3.0 / 8.0 * cp_L * (T_sat - T_metal)
                h_bar = 0.555 * ((9.81 * rho_sat_f * (rho_sat_f - rho_sat_g) * k_L**3.0 * h_fg) / (mu_L * abs(T_sat - T_metal) * D_pipe))**0.25
            else:
                h_bar = 0.0

    return h_bar


# ---------------------------------------------------------------------------
# Butterfly valve mass flow (ISA 75.01 liquid/gas/two-phase formulation)
# ---------------------------------------------------------------------------

def valve_massflow(CV: float, P_in: float, h_in: float, P_out: float) -> float:
    """
    Mass flow rate through a butterfly valve based on inlet/outlet pressures.

    Supports incompressible (subcooled liquid), compressible (superheated steam),
    and two-phase flow using the ISA 75.01 / Masoneilan formulation.

    Inputs
    ------
    CV    : float  Flow coefficient [gpm / sqrt(psi)]
    P_in  : float  Inlet pressure [Pa]
    h_in  : float  Inlet specific enthalpy [J/kg]
    P_out : float  Outlet pressure [Pa]

    Returns
    -------
    float
        Mass flow rate [kg/s]
    """
    P_kPa  = P_in / 1000.0
    T_in   = fp.temperature("water", P=P_kPa, h=h_in / 1000.0)
    rho_in = fp.density("water",     P=P_kPa, h=h_in / 1000.0)
    h_sat_f = fp.enthalpy("water", P=P_kPa, Q=0.0) * 1000.0
    T_sat   = fp.temperature("water", P=P_kPa, Q=0.0)
    h_sat_g = fp.enthalpy("water", P=P_kPa, Q=1.0) * 1000.0

    KV      = CV / 1.156  # converting english flow coefficient to metric flow coefficient
    DELTA_P = P_in - P_out

    if DELTA_P > 0.0:
        if h_in <= h_sat_f:  # flow is incompressible
            P_crit  = 22064000.0  # Critical Pressure of Water [Pa]
            P_ref   = 87726.1     # vapor pressure of steam at 100 degrees C [Pa]
            T_ref   = 373.0       # [K]
            # check if flow is choked or not choked
            lnP1P2 = 8.314 * (1.0 / T_ref - 1.0 / T_in)  # Clausius-Clapeyron Equation to find vapor pressure of water
            P_v    = P_ref * math.exp(lnP1P2)
            F_L    = 0.62
            F_F    = 0.96 - 0.28 * math.sqrt(P_v / P_crit)
            N1     = 0.1
            rho_ref = 1000.0  # reference density of water at room temp and pressure [kg/m^3]
            if DELTA_P < (F_L**2.0 * (P_in - F_F * P_v)):  # Flow is not choked
                Q_dot = KV * N1 / math.sqrt(rho_in / rho_ref / (DELTA_P / 1000.0))
                m_dot = Q_dot * rho_in / 3600.0
            else:  # Flow is choked
                Q_dot = KV * N1 * F_L / math.sqrt(rho_in / rho_ref / (P_in / 1000.0 - F_F * P_v / 1000.0))
                m_dot = Q_dot / rho_in / 3600.0

        elif h_in > h_sat_g:  # flow is compressible
            x_T = 0.35
            N6  = 3.16
            x   = DELTA_P / P_in
            if T_in > T_sat + 2.0:
                specheat_p = f_cp_water(P_in, T_in)
                specheat_v = f_cv_water(P_in, T_in)
            else:
                specheat_p = f_cp_water(P_in, T_sat + 2.0)
                specheat_v = f_cv_water(P_in, T_sat + 2.0)
            F_y = specheat_p / specheat_v / 1.4
            if x < F_y * x_T:  # Flow is not choked
                Y     = 1.0 - x / (3 * F_y * x_T)
                m_dot = KV * N6 * Y * math.sqrt(x * P_in / 1000.0 * rho_in) / 3600.0  # compressible flow equation, Pressure in kPa, mass flow in kg/hr
            else:  # Flow is choked
                Y     = 0.667
                m_dot = KV * N6 * math.sqrt(F_y * x_T * P_in / 1000.0 * rho_in) / 3600.0

        else:  # flow is two phase
            qual   = fp.quality("water", P=P_kPa, h=h_in / 1000.0)
            rho_in = fp.density("water", P=P_kPa, Q=1.0)  # use vapor density for compressible part
            # Solve for flow as if compressible
            x_T = 0.35
            N6  = 3.16
            x   = DELTA_P / P_in
            specheat_p = f_cp_water(P_in, T_in + 2.0)
            specheat_v = f_cv_water(P_in, T_in + 2.0)
            F_y = specheat_p / specheat_v / 1.4
            if x < F_y * x_T:  # Flow is not choked
                Y       = 1.0 - x / (3 * F_y * x_T)
                m_dot_c = KV * N6 * Y * math.sqrt(x * P_in / 1000.0 * rho_in) / 3600.0
            else:  # Flow is choked
                Y       = 0.667
                m_dot_c = KV * N6 * math.sqrt(F_y * x_T * P_in / 1000.0 * rho_in) / 3600.0
            # Solve for flow as if incompressible
            rho_in  = 1000.0
            P_crit  = 22064000.0  # Critical Pressure of Water [Pa]
            P_ref   = 87726.1     # vapor pressure of steam at 100 degrees C [Pa]
            T_ref   = 373.0
            # check if flow is choked or not choked
            lnP1P2  = 8.314 * (1.0 / T_ref - 1.0 / T_in)  # Clausius-Clapeyron Equation to find vapor pressure of water
            P_v     = P_ref * math.exp(lnP1P2)
            F_L     = 0.62
            F_F     = 0.96 - 0.28 * math.sqrt(P_v / P_crit)
            N1      = 0.1
            rho_ref = 1000.0  # reference density of water at room temp and pressure
            if DELTA_P < (F_L**2.0 * (P_in - F_F * P_v)):  # Flow is not choked
                Q_dot     = KV * N1 / math.sqrt(rho_in / rho_ref / (DELTA_P / 1000.0))
                m_dot_inc = Q_dot * rho_in / 3600.0
            else:  # Flow is choked
                Q_dot     = KV * N1 * F_L / math.sqrt(rho_in / rho_ref / (P_in / 1000.0 - F_F * P_v / 1000.0))
                m_dot_inc = Q_dot / rho_in / 3600.0
            # use quality to estimate flow through valve
            m_dot = m_dot_c * qual + m_dot_inc * (1.0 - qual)

    else:  # Pressure Gradient going backwards
        m_dot = 0.0

    return m_dot


# ---------------------------------------------------------------------------
# HTF property wrappers (delegate to Incompressible database)
# ---------------------------------------------------------------------------

def specheat(fnumd: str, T: float, P: float = 0.0) -> float:
    """
    Specific heat of an HTF identified by a fluid string.
    Delegates to the Incompressible property database.

    Inputs
    ------
    fnumd : str    Fluid identifier string (e.g., "Salt (60 NaNO3, 40 KNO3)")
    T     : float  Temperature [K]
    P     : float  Pressure [Pa] (unused for most HTF fluids)

    Returns
    -------
    float
        Specific heat [J/(kg·K)]

    Notes
    -----
    Incompressible.specheat returns kJ/(kg·K); multiplied by 1000 here for SI.
    """
    return _incompressible.specheat(fnumd, T, P) * 1000.0


def density(fnumd: str, T: float, P: float = 0.0) -> float:
    """
    Density of an HTF identified by a fluid string.
    Delegates to the Incompressible property database.

    Inputs
    ------
    fnumd : str    Fluid identifier string
    T     : float  Temperature [K]
    P     : float  Pressure [Pa] (required for Argon and Hydrogen)

    Returns
    -------
    float
        Density [kg/m^3]
    """
    return _incompressible.density(fnumd, T, P)
