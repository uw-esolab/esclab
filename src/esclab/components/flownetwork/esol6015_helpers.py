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

NOTE: specheat, Viscosity, and Density are NOT implemented here;
      use esclab.components.esol_properties.Incompressible instead.
"""

import numpy as np
from eeslib import fluid_properties as fp

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
