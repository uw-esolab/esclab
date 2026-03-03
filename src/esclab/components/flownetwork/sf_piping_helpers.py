"""
Solar field piping helper functions.

Converted from solar_field_modules.f90 (SF_piping_functions, Header_functions,
and Solar_Position modules).
All functions operate in SI units:
  - Lengths / diameters: [m]
  - Pressures: [Pa]
  - Temperatures: [K]
  - Mass flow rates: [kg/s]
  - Returned pressure drops: [Pa]
"""

import math
import numpy as np
from esclab.components.esol_properties import Incompressible as Inc


# ---------------------------------------------------------------------------
# Friction factor (from solar_field_modules.f90 :: FricFactor_piping)
# ---------------------------------------------------------------------------

def FricFactor_piping(Rough: float, Reynold: float) -> float:
    """Compute Darcy friction factor using Colebrook-White iteration.

    Uses an iterative method to solve the implicit friction factor function.
    For more on this method, refer to Fox, et al., 2006 Introduction to Fluid Mechanics.

    Parameters
    ----------
    Rough : float
        Relative roughness (ε/D) [-]
    Reynold : float
        Reynolds number [-]

    Returns
    -------
    float
        Darcy friction factor [-]
    """
    Acc = 0.01  # convergence tolerance (looser than FricFactor_IC)

    X = 33.33333  # 1. / 0.03
    TestOld = X + 2.0 * math.log10(Rough / 3.7 + 2.51 * X / max(Reynold, 1e-10))
    Xold = X
    X = 28.5714  # 1. / (0.03 + 0.005)
    NumTries = 0

    while True:
        NumTries += 1
        Test = X + 2.0 * math.log10(Rough / 3.7 + 2.51 * X / max(Reynold, 1e-10))
        if abs(Test - TestOld) <= Acc:
            return 1.0 / (X * X)

        if NumTries > 20:
            # Could not find friction factor solution; return last estimate
            return 1.0 / (X * X)

        Slope = (Test - TestOld) / (X - Xold)
        Xold = X
        TestOld = Test
        X = max((Slope * X - Test) / Slope, 1e-5)


# ---------------------------------------------------------------------------
# Pressure drop model (from solar_field_modules.f90 :: PressureDrop)
# ---------------------------------------------------------------------------

def PressureDrop(
    Fluid: float,
    m_dot: float,
    T: float,
    P: float,
    D: float,
    Rough: float,
    L_pipe: float,
    Nexp: float = 0.0,
    Ncon: float = 0.0,
    Nels: float = 0.0,
    Nelm: float = 0.0,
    Nell: float = 0.0,
    Ngav: float = 0.0,
    Nglv: float = 0.0,
    Nchv: float = 0.0,
    Nlw: float = 0.0,
    Nlcv: float = 0.0,
    Nbja: float = 0.0,
) -> float:
    """Compute total pressure drop across a length of pipe with fittings.

    Derived from the pressure drop calculations presented in:
    "Parabolic Trough Solar System Piping Model" (Kelly, Nexant Inc., NREL/SR-550-40165, 2006).

    This function should be called multiple times, once for each section under
    consideration (e.g., HCE, header sections, field piping).

    Parameters
    ----------
    Fluid : float
        Fluid identifier passed to SF_props (Incompressible)
    m_dot : float
        Mass flow rate [kg/s]
    T : float
        Fluid temperature [K]
    P : float
        Fluid pressure [Pa]
    D : float
        Pipe inner diameter [m]
    Rough : float
        Pipe absolute roughness [m]
    L_pipe : float
        Pipe length for pressure drop [m]
    Nexp : float
        Number of expansions [-]
    Ncon : float
        Number of contractions [-]
    Nels : float
        Number of standard elbows [-]
    Nelm : float
        Number of medium elbows [-]
    Nell : float
        Number of long elbows [-]
    Ngav : float
        Number of gate valves [-]
    Nglv : float
        Number of globe valves [-]
    Nchv : float
        Number of check valves [-]
    Nlw : float
        Number of loop weldolets [-]
    Nlcv : float
        Number of loop control valves [-]
    Nbja : float
        Number of ball joint assemblies [-]

    Returns
    -------
    float
        Total pressure drop [Pa]
    """
    pi = math.pi
    g = 9.80665

    # Calculate fluid properties and characteristics
    rho = Inc.density(Fluid, T=T, P=P)
    mu = Inc.viscosity(Fluid, T=T, P=P)
    nu = mu / rho
    v_dot = m_dot / rho                          # fluid volumetric flow rate [m^3/s]
    u_fluid = v_dot / (pi * (D / 2.0) ** 2)     # fluid mean velocity [m/s]

    # Dimensionless numbers
    Re = u_fluid * D / nu
    if Re < 2300.0:
        f = 64.0 / max(Re, 1.0)
    else:
        f = FricFactor_piping(Rough / D, Re)

    # Calculation of pressure loss from pipe length
    HL_pm = f * u_fluid * u_fluid / (2.0 * D * g)
    DP_pipe = HL_pm * rho * g * L_pipe

    # Calculation of pressure loss from fittings (K-factor method)
    DP_exp = 0.25 * rho * u_fluid * u_fluid * Nexp
    DP_con = 0.25 * rho * u_fluid * u_fluid * Ncon
    DP_els = 0.9 * D / f * HL_pm * rho * g * Nels
    DP_elm = 0.75 * D / f * HL_pm * rho * g * Nelm
    DP_ell = 0.6 * D / f * HL_pm * rho * g * Nell
    DP_gav = 0.19 * D / f * HL_pm * rho * g * Ngav
    DP_glv = 10.0 * D / f * HL_pm * rho * g * Nglv
    DP_chv = 2.5 * D / f * HL_pm * rho * g * Nchv
    DP_lw = 1.8 * D / f * HL_pm * rho * g * Nlw
    DP_lcv = 10.0 * D / f * HL_pm * rho * g * Nlcv
    DP_bja = 8.69 * D / f * HL_pm * rho * g * Nbja

    return sum([DP_pipe, DP_exp, DP_con, DP_els, DP_elm, DP_ell,
                DP_gav, DP_glv, DP_chv, DP_lw, DP_lcv, DP_bja])


# ---------------------------------------------------------------------------
# HTF pipe temperature derivative array (from solar_field_modules.f90 :: pipe_dTdt)
# ---------------------------------------------------------------------------

def pipe_dTdt(
    n_nodes: int,
    T: np.ndarray,
    Vol: float,
    mass_flow: float,
    mc_mult: float,
    Fluid_ID: float,
    heat_loss: float,
    L_cv: float,
) -> np.ndarray:
    """Compute nodal temperature time-derivatives for a 1-D pipe model.

    Parameters
    ----------
    n_nodes : int
        Number of nodes in the pipe
    T : np.ndarray, shape (n_nodes,)
        Array of nodal temperatures [K]
    Vol : float
        Volume of a typical control volume [m^3]
    mass_flow : float
        Mass flow rate through pipe [kg/s]
    mc_mult : float
        Thermal capacitance multiplier to account for piping and supports [-]
    Fluid_ID : float
        Fluid identifier for property lookups
    heat_loss : float
        Heat loss per unit length [W/m]
    L_cv : float
        Length of each control volume [m]

    Returns
    -------
    np.ndarray, shape (n_nodes,)
        Array of nodal temperature time-derivatives [K/s]
    """
    dt = np.zeros(n_nodes)
    dtdt_bar = np.zeros(n_nodes - 1)

    # Loop through control volumes to compute CV temperature rate
    for n in range(n_nodes - 1):
        # Heat loss
        q_out = heat_loss * L_cv
        # Compute CV average temperature and properties
        T_ave = (T[n] + T[n + 1]) / 2.0
        rho = Inc.density(Fluid_ID, T=T_ave, P=0.0)
        c = Inc.specheat(Fluid_ID, T=T_ave, P=0.0)  # [J/kg-K]
        # Compute CV temperature rate
        dtdt_bar[n] = 1.0 / Vol / rho / c / mc_mult * (mass_flow * c * (T[n] - T[n + 1]) - q_out)

    # Compute nodal temperature rates
    dt[0] = 0.0                     # First node is set by input to type
    dt[n_nodes - 1] = dtdt_bar[n_nodes - 2]  # Cannot average last node
    for n in range(1, n_nodes - 1):  # Remaining middle nodes
        dt[n] = (dtdt_bar[n - 1] + dtdt_bar[n]) / 2.0

    return dt


# ---------------------------------------------------------------------------
# Row shadow efficiency (from solar_field_modules.f90 :: Solar_Position)
# ---------------------------------------------------------------------------

def Row_shadow(phi: float, row_distance: float, w_ap: float) -> float:
    """Compute row-to-row shadow efficiency for a parabolic trough field.

    Parameters
    ----------
    phi : float
        Collector tracking / row shadow angle [rad]
    row_distance : float
        Row-to-row distance [m]
    w_ap : float
        Collector aperture width [m]

    Returns
    -------
    float
        Row shadow efficiency factor [-], clamped to [0, 1]
    """
    eta_row = abs(math.cos(phi)) * row_distance / w_ap

    if eta_row > 1.0:
        eta_row = 1.0
    elif eta_row < 0.0:
        eta_row = 0.0

    return eta_row


# ---------------------------------------------------------------------------
# Enthalpy of Dowtherm A (from solar_field_modules.f90 :: Header_functions)
# ---------------------------------------------------------------------------

def H_Dowtherm_A(T: float) -> float:
    """Compute specific enthalpy of Dowtherm A.

    Reference: Dowtherm A heat transfer fluid TDS, Applied Thermal Fluids.
    https://www.appliedthermalfluids.com/wp-content/uploads/2018/02/Dowtherm-A-heat-transfer-fluid-TDS.pdf

    Parameters
    ----------
    T : float
        Temperature [K]

    Returns
    -------
    float
        Specific enthalpy [J/kg]
    """
    Td = T - 273.15  # [C]
    return (-12.7078 + 1.481714 * Td + 0.0014292857 * Td ** 2) * 1000  # [J/kg]


# ---------------------------------------------------------------------------
# Header geometry functions (from solar_field_modules.f90 :: Header_functions)
# ---------------------------------------------------------------------------

def diams_inlet(n_cv: int, n_loop: int, L_row: float, L_exp: float, diam_file: str) -> np.ndarray:
    """Build array of inner diameters for each inlet-header control volume.

    Reads inlet and return header diameters from the geometry file.  The file
    is expected to have one row per loop-pair with two whitespace-delimited
    columns: inlet_diameter [m] and return_diameter [m].

    Parameters
    ----------
    n_cv : int
        Number of inlet-header control volumes [-]
    n_loop : int
        Total number of loops in the sector [-]
    L_row : float
        Row-to-row distance [m]
    L_exp : float
        Expansion loop length [m]
    diam_file : str
        Path to the header geometry file

    Returns
    -------
    np.ndarray, shape (n_cv,)
        Inner diameter of each inlet-header control volume [m]
    """
    Diam = np.zeros(n_loop // 2)

    # Get specified geometries
    with open(diam_file) as fh:
        for cc in range(n_loop // 2):
            vals = fh.readline().strip().split()
            Diam[cc] = float(vals[0])  # inlet diameter in first column

    Diam_cv = np.zeros(n_cv)
    cc = 1
    loop_count = 0
    # Loop through all control volumes
    for n in range(n_cv):
        D_curr = Diam[loop_count]
        # Check if control volume is part of an expansion loop
        if cc > 1:
            Diam_cv[n] = D_curr
            cc += 1
            if cc == 4:
                cc = 1
                loop_count += 1
        # Else control volume is the segment of header between sf loop inlets
        else:
            Diam_cv[n] = D_curr
            cc += 1
            loop_count += 1

    return Diam_cv


def vols_inlet(n_cv: int, n_loop: int, L_row: float, L_exp: float, diam_file: str) -> np.ndarray:
    """Build array of volumes for each inlet-header control volume.

    This function creates an array with the volumes of each control volume
    in an inlet header.

    Parameters
    ----------
    n_cv : int
        Number of inlet-header control volumes [-]
    n_loop : int
        Total number of loops in the sector [-]
    L_row : float
        Row-to-row distance [m]
    L_exp : float
        Expansion loop length [m]
    diam_file : str
        Path to the header geometry file

    Returns
    -------
    np.ndarray, shape (n_cv,)
        Volume of each inlet-header control volume [m^3]
    """
    Diam = np.zeros(n_loop // 2)

    # Get specified geometries
    with open(diam_file) as fh:
        for cc in range(n_loop // 2):
            vals = fh.readline().strip().split()
            Diam[cc] = float(vals[0])  # inlet diameter in first column

    pi = 3.141529
    vols = np.zeros(n_cv)
    cc = 1
    loop_count = 0
    # Loop through all control volumes
    for n in range(n_cv):
        D_curr = Diam[loop_count]
        # Check if control volume is part of an expansion loop
        if cc > 1:
            vols[n] = (D_curr / 2) ** 2 * pi * (L_row * 2 + L_exp * 2) / 2
            cc += 1
            if cc == 4:
                cc = 1
                loop_count += 1
        # Else control volume is the pipe between two sf loop returns
        else:
            vols[n] = (D_curr / 2) ** 2 * pi * L_row * 2
            cc += 1
            loop_count += 1

    return vols


def diams_return(n_cv: int, n_loop: int, L_row: float, L_exp: float, diam_file: str) -> np.ndarray:
    """Build array of inner diameters for each return-header control volume.

    Reads inlet and return header diameters from the geometry file.  The file
    is expected to have one row per loop-pair with two whitespace-delimited
    columns: inlet_diameter [m] and return_diameter [m].

    Parameters
    ----------
    n_cv : int
        Number of return-header control volumes [-]
    n_loop : int
        Total number of loops in the sector [-]
    L_row : float
        Row-to-row distance [m]
    L_exp : float
        Expansion loop length [m]
    diam_file : str
        Path to the header geometry file

    Returns
    -------
    np.ndarray, shape (n_cv,)
        Inner diameter of each return-header control volume [m]
    """
    Diam = np.zeros(n_loop // 2)

    # Get specified geometries
    with open(diam_file) as fh:
        for cc in range(n_loop // 2):
            vals = fh.readline().strip().split()
            Diam[cc] = float(vals[1])  # return diameter in second column

    Diam_cv = np.zeros(n_cv)
    cc = 1
    loop_count = 0
    # Loop through all control volumes
    for n in range(n_cv):
        D_curr = Diam[loop_count]
        # Check if control volume is part of an expansion loop
        if cc > 1:
            Diam_cv[n] = D_curr
            cc += 1
            if cc == 4:
                cc = 1
                loop_count += 1
        # Else control volume is the segment of header between sf loop inlets
        else:
            Diam_cv[n] = D_curr
            cc += 1
            loop_count += 1

    return Diam_cv


def vols_return(n_cv: int, n_loop: int, L_row: float, L_exp: float, diam_file: str) -> np.ndarray:
    """Build array of volumes for each return-header control volume.

    This function creates an array with the volumes of each control volume
    in a return header.

    Parameters
    ----------
    n_cv : int
        Number of return-header control volumes [-]
    n_loop : int
        Total number of loops in the sector [-]
    L_row : float
        Row-to-row distance [m]
    L_exp : float
        Expansion loop length [m]
    diam_file : str
        Path to the header geometry file

    Returns
    -------
    np.ndarray, shape (n_cv,)
        Volume of each return-header control volume [m^3]
    """
    Diam = np.zeros(n_loop // 2)

    # Get specified geometries
    with open(diam_file) as fh:
        for cc in range(n_loop // 2):
            vals = fh.readline().strip().split()
            Diam[cc] = float(vals[1])  # return diameter in second column

    pi = 3.141529
    vols = np.zeros(n_cv)
    cc = 1
    loop_count = 0
    for n in range(n_cv):
        D_curr = Diam[loop_count]

        if cc > 1:
            vols[n] = (D_curr / 2) ** 2 * pi * (L_row * 2 + L_exp * 2) / 2
            cc += 1
            if cc == 4:
                cc = 1
                loop_count += 1
        else:
            vols[n] = (D_curr / 2) ** 2 * pi * L_row * 2
            cc += 1
            loop_count += 1

    return vols


# ---------------------------------------------------------------------------
# Header temperature derivative functions (from solar_field_modules.f90 :: Header_functions)
# ---------------------------------------------------------------------------

def dT_dt_inlet(
    m_dots: np.ndarray,
    t: np.ndarray,
    Vols: np.ndarray,
    L_cv: np.ndarray,
    mc: float,
    t_bar: np.ndarray,
    n_nodes: int,
    fluid: float,
    heat_loss: float,
) -> np.ndarray:
    """Compute nodal temperature time-derivatives for the inlet header.

    This function computes the nodal temperature rate of change for the
    inlet header using enthalpy-based energy balance on each control volume.

    Parameters
    ----------
    m_dots : np.ndarray, shape (n_nodes-1,)
        Mass flow rates through each inlet-header control volume [kg/s]
    t : np.ndarray, shape (n_nodes,)
        Nodal temperatures [K]
    Vols : np.ndarray, shape (n_nodes-1,)
        Control volume volumes [m^3]
    L_cv : np.ndarray, shape (n_nodes-1,)
        Control volume lengths [m]
    mc : float
        Thermal capacitance multiplier (header wall + fluid) [-]
    t_bar : np.ndarray, shape (n_nodes-1,)
        Control volume average temperatures [K]
    n_nodes : int
        Number of header nodes [-]
    fluid : float
        Fluid identifier for property lookups
    heat_loss : float
        Heat loss coefficient per unit length [W/m]

    Returns
    -------
    np.ndarray, shape (n_nodes,)
        Nodal temperature time-derivatives [K/s]
    """
    dT_bar = np.zeros(n_nodes - 1)
    dT = np.zeros(n_nodes)

    # Compute control volume temperature rate of change
    for n in range(n_nodes - 1):
        c = Inc.specheat(fluid=fluid, T=t_bar[n], P=0.0)  # [J/kg-K]
        rho = Inc.density(fluid=fluid, T=t_bar[n], P=0.0)
        h1 = H_Dowtherm_A(t[n])
        h2 = H_Dowtherm_A(t[n + 1])
        dT_bar[n] = 1.0 / mc / (Vols[n] * rho * c) * (m_dots[n] * (h1 - h2) - heat_loss * L_cv[n])

    # Compute nodal temperature rate of change
    dT[0] = 0.0
    for n in range(1, n_nodes - 1):
        dT[n] = 0.5 * (dT_bar[n - 1] + dT_bar[n])
    dT[n_nodes - 1] = dT_bar[n_nodes - 2]

    return dT


def dT_dt_return(
    m_dot_l: np.ndarray,
    m_dot_r: np.ndarray,
    m_dot_out: np.ndarray,
    T_l: np.ndarray,
    T_r: np.ndarray,
    T: np.ndarray,
    T_bar: np.ndarray,
    Vols: np.ndarray,
    L_cv: np.ndarray,
    mc: float,
    n_nodes: int,
    n_loop: int,
    inds_header: np.ndarray,
    fluid: float,
    heat_loss: float,
) -> np.ndarray:
    """Compute nodal temperature time-derivatives for the return header.

    This function computes the nodal temperature rate of change for the
    return header, accounting for HTF contributions from solar field loops
    at each mixing junction node.

    Parameters
    ----------
    m_dot_l : np.ndarray, shape (n_loop/2,)
        Mass flow from left-side loops at each junction [kg/s]
    m_dot_r : np.ndarray, shape (n_loop/2,)
        Mass flow from right-side loops at each junction [kg/s]
    m_dot_out : np.ndarray, shape (n_nodes-1,)
        Mass flow rate leaving each return-header control volume [kg/s]
    T_l : np.ndarray, shape (n_loop/2,)
        Outlet temperature of left-side loops at each junction [K]
    T_r : np.ndarray, shape (n_loop/2,)
        Outlet temperature of right-side loops at each junction [K]
    T : np.ndarray, shape (n_nodes,)
        Nodal temperatures [K]
    T_bar : np.ndarray, shape (n_nodes-1,)
        Control volume average temperatures [K]
    Vols : np.ndarray, shape (n_nodes-1,)
        Control volume volumes [m^3]
    L_cv : np.ndarray, shape (n_nodes-1,)
        Control volume lengths [m]
    mc : float
        Thermal capacitance multiplier (header wall + fluid) [-]
    n_nodes : int
        Number of header nodes [-]
    n_loop : int
        Total number of loops in sector [-]
    inds_header : np.ndarray, shape (n_loop/2,)
        0-based CV indices where SF loop outlets join the return header
    fluid : float
        Fluid identifier for property lookups
    heat_loss : float
        Heat loss coefficient per unit length [W/m]

    Returns
    -------
    np.ndarray, shape (n_nodes,)
        Nodal temperature time-derivatives [K/s]
    """
    dT_bar = np.zeros(n_nodes - 1)
    dT = np.zeros(n_nodes)

    # Compute control volume temperature rate of changes
    # jj starts at index 1 (0-based), matching Fortran jj=2 (1-based) for the second junction onward
    jj = 1
    for n in range(n_nodes - 1):
        rho = Inc.density(fluid=fluid, T=T_bar[n], P=0.0)
        c = Inc.specheat(fluid=fluid, T=T_bar[n], P=0.0)  # [J/kg-K]
        h1 = H_Dowtherm_A(T[n])
        h2 = H_Dowtherm_A(T[n + 1])
        if n == 0:
            dT_bar[n] = 1.0 / mc / (Vols[n] * rho * c) * (m_dot_out[n] * (h1 - h2))
        else:
            if n == inds_header[jj]:
                h_r = H_Dowtherm_A(T_r[jj])
                h_l = H_Dowtherm_A(T_l[jj])
                dT_bar[n] = 1.0 / mc / (Vols[n] * rho * c) * (
                    m_dot_l[jj] * h_l + m_dot_r[jj] * h_r
                    + m_dot_out[n - 1] * h1 - m_dot_out[n] * h2
                    - heat_loss * L_cv[n]
                )
                jj += 1
            else:
                dT_bar[n] = 1.0 / Vols[n] / rho / c * (
                    m_dot_out[n] * (h1 - h2) - heat_loss * L_cv[n]
                )

    # Compute nodal temperature rate of changes
    dT[0] = 0.0
    for n in range(1, n_nodes - 1):
        dT[n] = 0.5 * (dT_bar[n - 1] + dT_bar[n])
    dT[n_nodes - 1] = dT_bar[n_nodes - 2]

    return dT
