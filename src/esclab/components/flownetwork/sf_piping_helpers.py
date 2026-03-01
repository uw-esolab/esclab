"""
Solar field piping helper functions.

Converted from solar_field_modules.f90 (SF_piping_functions and related modules).
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
        c = Inc.specheat(Fluid_ID, T=T_ave, P=0.0) * 1000.0  # [J/kg-K], spec_SF returns kJ/kg-K
        # Compute CV temperature rate
        dtdt_bar[n] = 1.0 / Vol / rho / c / mc_mult * (mass_flow * c * (T[n] - T[n + 1]) - q_out)

    # Compute nodal temperature rates
    dt[0] = 0.0                     # First node is set by input to type
    dt[n_nodes - 1] = dtdt_bar[n_nodes - 2]  # Cannot average last node
    for n in range(1, n_nodes - 1):  # Remaining middle nodes
        dt[n] = (dtdt_bar[n - 1] + dtdt_bar[n]) / 2.0

    return dt
