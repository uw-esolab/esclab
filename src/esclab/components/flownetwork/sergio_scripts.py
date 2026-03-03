"""Python conversion of the TRNSYS SergioScripts Fortran module.

Contains helper subroutines for shell-and-tube heat exchanger calculations
including correlation functions, 0D/1D sub-solvers, and numerical utilities.

Original Fortran module: SergioScripts
Author (original):    Sergio Alcalde Morales
Converted by:         GitHub Copilot, March 02, 2026
"""

import numpy as np

from esclab.components.esol_properties import Incompressible as Inc


# ---------------------------------------------------------------------------
# CORRELATIONS
# ---------------------------------------------------------------------------

def gnielinski(T_fluid, m_tube, n_tube, D_tube, fluid, pressure):
    """Gnielinski correlation for internal pipe flow heat-transfer coefficient.

    Parameters
    ----------
    T_fluid     : float   - bulk fluid temperature [K]
    m_tube      : float   - total mass-flow rate through all tubes [kg/s]
    n_tube      : float   - number of tubes in the bundle [-]
    D_tube      : float   - tube inner diameter [m]
    fluid       : str     - fluid identifier passed to Incompressible properties
    pressure    : float   - fluid pressure [Pa]

    Returns
    -------
    h_tube  : float   - convection coefficient [W/m²·K]
    vel_tube: float   - mean fluid velocity [m/s]
    """
    rho_tube = Inc.density(fluid=fluid, T=T_fluid, P=pressure)
    visc_tube = Inc.viscosity(fluid=fluid, T=T_fluid, P=pressure)
    cp_tube = Inc.specheat(fluid=fluid, T=T_fluid, P=pressure)
    # TODO-NEEDS REVIEW: (MJW) the value k=.1 is too low if the fluid is salt. If therminol (or oil), should probably be 0.07
    k_film = 0.1  # Constant in the whole code
    Pr_tube = visc_tube * cp_tube / k_film

    m_dot_tube = m_tube / n_tube
    vel_tube = m_dot_tube / (rho_tube * (np.pi / 4.0) * D_tube ** 2)

    Re_tube = (rho_tube * vel_tube * D_tube) / visc_tube
    if Re_tube > 4000.0:  # Gnielinski Method Valid
        # Steel properties:
        f = 0.035
        # Nusselt number
        Nu_tube = ((f / 8.0) * (Re_tube - 1000.0) * Pr_tube) / (
            1.0 + (12.7 * np.sqrt(f / 8.0) * ((Pr_tube ** (2.0 / 3.0)) - 1.0))
        )
    else:  # Gnielinski Method Not Valid, estimate Nusselt Number using Nellis Heat Transfer Book pp. 662, uniform wall flux
        Nu_tube = 4.36

    h_tube = Nu_tube * k_film / D_tube
    return h_tube, vel_tube


def zukauskas(fluid_shell, config, D_shell, S_T, S_L, n_tube, D_tube,
              T_shell, m_shell, L_Shell, n_baffles, pressure):
    """Zukauskas correlation for external cross-flow over a tube bundle.

    Parameters
    ----------
    fluid_shell : str    - fluid identifier for shell-side fluid
    config      : bool   - tube layout flag (True=inline, False=staggered)
    D_shell     : float  - shell inner diameter [m]
    S_T         : float  - transverse pitch [m]
    S_L         : float  - longitudinal pitch [m]
    n_tube      : float  - number of tubes in the bundle [-]
    D_tube      : float  - tube outer diameter [m]
    T_shell     : float  - shell-side fluid temperature [K]
    m_shell     : float  - shell-side mass-flow rate [kg/s]
    L_Shell     : float  - shell length [m]
    n_baffles   : int    - number of baffles [-]
    pressure    : float  - shell-side pressure [Pa]

    Returns
    -------
    h_shell  : float  - shell-side convection coefficient [W/m²·K]
    vel_shell: float  - representative shell-side velocity [m/s]
    """
    rho_shell = Inc.density(fluid=fluid_shell, T=T_shell, P=pressure)
    visc_shell = Inc.viscosity(fluid=fluid_shell, T=T_shell, P=pressure)
    cp_shell = Inc.specheat(fluid=fluid_shell, T=T_shell, P=pressure)
    k_film = 0.5  # W/m-K
    Pr_shell = visc_shell * cp_shell / k_film
    Pr_film = Pr_shell

    gap_baffle = L_Shell / (n_baffles + 1)
    vel_shell = 3.0 * m_shell / rho_shell / (np.pi * D_shell * gap_baffle / 2)
    Re_shell = (rho_shell * vel_shell * D_tube) / visc_shell

    # Zukauskas coefficients as a function of the alignment of the tubes and Reynolds number
    if config:
        if Re_shell < 100:
            C = 0.8
            m = 0.4
        elif Re_shell < 1000:  # Approximate as a single cylinder
            C = 0.51
            m = 0.50
        elif Re_shell < 2e5:
            C = 0.27
            m = 0.63
        else:
            C = 0.021
            m = 0.84
    else:  # Staggered configuration
        if Re_shell < 100:
            C = 0.9
            m = 0.4
        elif Re_shell < 1000:  # Approximate as a single cylinder
            C = 0.51
            m = 0.5
        elif Re_shell < 2e5:
            if (S_T / S_L) > 2:
                C = 0.4
                m = 0.60
            else:
                C = 0.35 * (S_T / S_L) ** 0.2
                m = 0.60
        else:
            C = 0.022
            m = 0.84

    # Zukauskas coefficients as a function of the number of tubes
    if Re_shell > 1000 and n_tube < 16:
        if config:
            C2 = 0.6233 + 0.1007 * n_tube - 0.0095 * n_tube ** 2 + (2.8849e-4) * n_tube ** 3
        else:
            C2 = 0.5385 + 0.1292 * n_tube - 0.0123 * n_tube ** 2 + (3.7849e-4) * n_tube ** 3
    else:
        C2 = 1.0

    Nu_shell = C * Re_shell ** m * Pr_shell ** 0.36 * (Pr_shell / Pr_film) ** 0.25
    Nu_shell = Nu_shell * C2
    h_shell = Nu_shell * k_film / D_tube  # Lc = D_tubes
    if h_shell == 0.0:
        h_shell = 0.01  # done to avoid division by zero

    return h_shell, vel_shell


def taborek(fluid_shell, D_shell, D_tube, T_shell, m_shell, Lp, pressure):
    """Taborek correlation for counterflow in a tube bundle.  Lp = S_L.

    Parameters
    ----------
    fluid_shell : str    - fluid identifier for shell-side fluid
    D_shell     : float  - shell inner diameter [m]
    D_tube      : float  - tube outer diameter [m]
    T_shell     : float  - shell-side fluid temperature [K]
    m_shell     : float  - shell-side mass-flow rate [kg/s]
    Lp          : float  - longitudinal pitch S_L [m]
    pressure    : float  - shell-side pressure [Pa]

    Returns
    -------
    h_U : float  - U-section convection coefficient [W/m²·K]
    """
    rho_U = Inc.density(fluid=fluid_shell, T=T_shell, P=pressure)
    visc_U = Inc.viscosity(fluid=fluid_shell, T=T_shell, P=pressure)
    cp_U = Inc.specheat(fluid=fluid_shell, T=T_shell, P=pressure)
    k_film = 0.5
    Pr_shell = visc_U * cp_U / k_film

    vel_U = m_shell / (rho_U * (np.pi / 4) * D_shell ** 2)
    Re_U = (rho_U * vel_U * D_tube) / visc_U

    # TODO-NEEDS CONVERSION REVIEW: Pr_U used in Nu_U_t but never assigned; Pr_shell computed but not used in Nu_U_t
    Pr_U = Pr_shell  # assumed Pr_U == Pr_shell (Pr_U was unassigned in original Fortran)
    Nu_U_t = 0.023 * Re_U ** 0.8 * Pr_U ** 0.4
    PR = Lp / D_tube
    C_i = -0.1315 + 1.5658 * PR - 0.415 * PR
    Nu_U = Nu_U_t * C_i

    h_U = Nu_U * k_film / D_tube
    return h_U


# ---------------------------------------------------------------------------
# ZERO-D HEAT TRANSFER
# ---------------------------------------------------------------------------

def ZeroD_Eq(n_nodes, h_tube, h_shell, T_shell, T_tube, k_tube, Delta_x,
             d_in, d_out, q_HTF_tube, q_HTF_shell):
    """Calculate the heat transfer from the shell to the tubes by using resistances.

    Modifies q_HTF_tube and q_HTF_shell in-place.

    Parameters
    ----------
    n_nodes     : int             - number of nodes
    h_tube      : ndarray[n_nodes] – tube-side convection coefficients [W/m²·K]
    h_shell     : ndarray[n_nodes] – shell-side convection coefficients [W/m²·K]
    T_shell     : ndarray[n_nodes] – shell-side temperatures [K]
    T_tube      : ndarray[n_nodes] – tube-side temperatures [K]
    k_tube      : float            - tube-wall conductivity [W/m·K]
    Delta_x     : float            - axial node length [m]
    d_in        : float            - tube inner diameter [m]
    d_out       : float            - tube outer diameter [m]
    q_HTF_tube  : ndarray[n_nodes] – (out) heat flux at tube side [W]
    q_HTF_shell : ndarray[n_nodes] – (out) heat flux at shell side [W]
    """
    A_out = d_out * np.pi * Delta_x
    A_in = d_in * np.pi * Delta_x
    # Convection resistance in the shell part (h_shell obtained in shell subroutine):
    R_conv_shell = 1.0 / A_out / h_shell
    # Convection resistance in the tube part (h_tube obtained in tube code)
    R_conv_tube = 1.0 / A_in / h_tube
    # Conduction resistance
    r_in = d_in / 2.0
    r_out = d_out / 2.0
    R_cond = np.log(r_out / r_in) / (2.0 * np.pi * k_tube * Delta_x)

    # Fouling resistance: Taken from the datasheet
    R_fouling_out = (0.0005 * 0.3048 ** 2 * 3.41 * 0.5556) / A_out  # [K/W]
    R_fouling_in = (0.0005 * 0.3048 ** 2 * 3.41 * 0.5556) / A_in    # [K/W]

    # Heat flux from the shell to the Heat Transfer Fluid (q_HTF [W])
    q_HTF_tube[:] = h_tube
    q_HTF_shell[:] = h_tube

    for i in range(1, n_nodes + 1):
        q_HTF_shell[i - 1] = (T_shell[i - 1] - T_tube[n_nodes - i]) / (
            R_conv_shell[i - 1] + R_conv_tube[n_nodes - i] + R_cond + R_fouling_in + R_fouling_out
        )
        q_HTF_tube[n_nodes - i] = (T_shell[i - 1] - T_tube[n_nodes - i]) / (
            R_conv_shell[i - 1] + R_conv_tube[n_nodes - i] + R_cond + R_fouling_in + R_fouling_out
        )


# ---------------------------------------------------------------------------
# E-NTU INITIALISATION
# ---------------------------------------------------------------------------

def eNTU(p_in_shell, p_in_tube, n_nodes, shell_p, n_tube, r_out, r_in,
         L_tube, D_shell, T_tube_in, T_shell_in, m_tube, m_shell,
         fluid_shell, config, S_T, S_L, fluid_tube, n_HEX, n_baffles,
         T_shell, T_tube, h_shell, h_tube):
    """Subroutine eNTU is useful to initialise problems of Shell-tube HEX.

    Fills T_shell, T_tube, h_shell, h_tube (length n_nodes) in-place
    with a steady-state e-NTU solution.

    Parameters
    ----------
    p_in_shell  : float  - shell-side inlet pressure [Pa]
    p_in_tube   : float  - tube-side inlet pressure [Pa]
    n_nodes     : int    - number of axial nodes
    shell_p     : int    - number of shell passes per HEX
    n_tube      : float  - number of tubes
    r_out       : float  - tube outer radius [m]
    r_in        : float  - tube inner radius [m]
    L_tube      : float  - tube length [m]
    D_shell     : float  - shell inner diameter [m]
    T_tube_in   : float  - tube-side inlet temperature [K]
    T_shell_in  : float  - shell-side inlet temperature [K]
    m_tube      : float  - tube-side mass-flow rate [kg/s]
    m_shell     : float  - shell-side mass-flow rate [kg/s]
    fluid_shell : str    - shell-side fluid identifier
    config      : bool   - tube layout flag (True=inline, False=staggered)
    S_T         : float  - transverse pitch [m]
    S_L         : float  - longitudinal pitch [m]
    fluid_tube  : str    - tube-side fluid identifier
    n_HEX       : int    - number of HEX in series
    n_baffles   : int    - number of baffles
    T_shell     : ndarray[n_nodes] – (out) shell-side temperature profile [K]
    T_tube      : ndarray[n_nodes] – (out) tube-side temperature profile [K]
    h_shell     : ndarray[n_nodes] – (out) shell-side convection profile [W/m²·K]
    h_tube      : ndarray[n_nodes] – (out) tube-side convection profile [W/m²·K]
    """
    # Specific heat at a representative temperature
    cp_tube_st = Inc.specheat(fluid=fluid_tube, T=T_tube_in, P=p_in_tube)
    cp_shell_st = Inc.specheat(fluid=fluid_shell, T=T_shell_in, P=p_in_shell)

    if cp_tube_st < cp_shell_st:
        c_r = cp_tube_st / cp_shell_st
    else:
        c_r = cp_shell_st / cp_tube_st  # Specific heat of the limiting fluid

    D_out = 2.0 * r_out
    D_in = 2.0 * r_in
    D_tube = D_out
    th_tube = r_out - r_in
    A_out = np.pi * D_out * L_tube * n_tube
    A_in = np.pi * D_in * L_tube * n_tube

    shell_passes = shell_p * n_HEX

    # Outlet temperatures are supposed
    T_tube_out = 400.0
    T_shell_out = 500.0

    error = 1.0

    while error > 0.0000001:
        # Caloric temperatures
        Tcal_tube = (T_tube_out + T_tube_in) / 2
        Tcal_shell = (T_shell_out + T_shell_in) / 2

        # For this first calculation, properties are obtained at caloric temperatures
        k_film = 0.1  # Constant in all calculation

        h_shell_i, vel_shell = zukauskas(
            fluid_shell, config, D_shell, S_T, S_L, n_tube, D_tube,
            Tcal_shell, m_shell, L_tube, n_baffles, p_in_shell
        )

        h_tube_i, vel_tube = gnielinski(
            Tcal_tube, m_tube, n_tube, D_tube, fluid_tube, p_in_tube
        )

        # Fouling resistances (From datasheet in IS units)
        FR_tube = 0.0005 * 0.092903 * 3.41 / 1.8
        FR_shell = 0.0005 * 0.092903 * 3.41 / 1.8

        # Global coefficient U
        U_global = (1 / A_out) * (
            1 / (1 / h_tube_i / A_in + 1 / h_shell_i / A_out + th_tube / k_film / A_out
                 + FR_shell / A_out + FR_tube / A_in)
        )
        # Number of transfer units
        NTU = U_global * A_out / c_r

        # Effectiveness for 1 shell pass
        epsi_1 = 2 * (
            1 + c_r + ((1 + c_r ** 2) ** 0.5)
            * (1 + np.exp(-NTU * (1 + c_r ** 2) ** 0.5))
            / (1 - np.exp(-NTU * (1 + c_r ** 2) ** 0.5))
        ) ** (-1)

        # Effectiveness for all shell passes
        efectivity = (
            ((1 - epsi_1 * c_r) / (1 - epsi_1)) ** shell_passes - 1
        ) * (
            ((1 - epsi_1 * c_r) / (1 - epsi_1)) ** shell_passes - c_r
        ) ** (-1)

        # Outlet temperature
        DeltaT_max = T_shell_in - T_tube_in
        T_s_out = T_shell_out
        T_t_out = T_tube_out

        if cp_shell_st < cp_tube_st:
            T_shell_out = T_shell_in - efectivity * DeltaT_max
            T_tube_out = T_tube_in + c_r * (T_shell_in - T_shell_out)
        else:
            T_tube_out = T_tube_in + efectivity * DeltaT_max
            T_shell_out = T_shell_in - c_r * (T_tube_out - T_tube_in)

        error_shell = np.sqrt((T_s_out - T_shell_out) ** 2)
        error_tube = np.sqrt((T_tube_out - T_t_out) ** 2)

        error = error_tube + error_shell

    # Here initial profile of temperatures is defined (Linear profile)
    T_shell[0] = T_shell_in
    T_tube[0] = T_tube_in

    for i in range(2, n_nodes - 1 + 1):
        T_shell[i - 1] = T_shell[i - 2] + (T_shell_out - T_shell_in) * (1.0 / (n_nodes - 1))
        T_tube[i - 1] = T_tube[i - 2] + (T_tube_out - T_tube_in) * (1.0 / (n_nodes - 1))

    T_shell[n_nodes - 1] = T_shell_out
    T_tube[n_nodes - 1] = T_tube_out

    for i in range(1, n_nodes + 1):
        h_shell[i - 1] = h_shell_i
        h_tube[i - 1] = h_tube_i


# ---------------------------------------------------------------------------
# SHELL SUBROUTINE
# ---------------------------------------------------------------------------

def shell(config, shell_p, Tenv, P_in_shell, gap_baffle, n_nodes,
          k_steel, k_ins, S_T, S_L, m_shell, n_baffles, n_HEX,
          th_baffle, th_ins, T_shell, D_shell, L_baffle, A_baffle,
          shell_passes, L_Shell, n_tubes, D_tube, q_HTF_shell,
          v_out, Delta_t, fluid_shell,
          T_shell_new, h_shell, vel_shell):
    """Subroutine shell: Transient energy balance in a HEX in the shell part.

    Modifies T_shell_new, h_shell, vel_shell in-place.

    Parameters
    ----------
    config       : bool              - tube layout flag
    shell_p      : int               - number of shell passes
    Tenv         : float             - ambient temperature [K]
    P_in_shell   : float             - shell-side inlet pressure [Pa]
    gap_baffle   : float             - baffle spacing [m]
    n_nodes      : int               - number of axial nodes
    k_steel      : float             - steel thermal conductivity [W/m·K]
    k_ins        : float             - insulation thermal conductivity [W/m·K]
    S_T          : float             - transverse pitch [m]
    S_L          : float             - longitudinal pitch [m]
    m_shell      : float             - shell-side mass-flow rate [kg/s]
    n_baffles    : int               - number of baffles
    n_HEX        : int               - number of HEX in series
    th_baffle    : float             - baffle thickness fraction [-]
    th_ins       : float             - insulation thickness [m]
    T_shell      : ndarray[n_nodes]  - shell-side temperatures from last step [K]
    D_shell      : float             - shell inner diameter [m]
    L_baffle     : float             - baffle length [m]
    A_baffle     : float             - baffle cross-sectional area [m²]
    shell_passes : int               - total shell passes
    L_Shell      : float             - shell length [m]
    n_tubes      : float             - number of tubes [-]
    D_tube       : float             - tube outer diameter [m]
    q_HTF_shell  : ndarray[n_nodes]  - heat flux term from ZeroD_Eq [W]
    v_out        : float             - external wind velocity [m/s]
    Delta_t      : float             - timestep [s]
    fluid_shell  : str               - shell-side fluid identifier
    T_shell_new  : ndarray[n_nodes]  - (out) updated shell temperatures [K]
    h_shell      : ndarray[n_nodes]  - (out) updated h-values [W/m²·K]
    vel_shell    : ndarray[n_nodes]  - (out) shell-side velocities [m/s]

    Returns
    -------
    p_shell_out : float – shell-side outlet pressure [Pa]
    """
    # Volume of the tubes at each control volume [m3]
    V_tubes = (np.pi * D_tube ** 2 / 4.0) * gap_baffle * n_tubes
    # Control volume [m3]
    Vc = (np.pi * D_shell ** 2 / 4.0) * gap_baffle / shell_p - V_tubes
    g = 9.8  # Gravity
    # Total number of regions, distinguished
    n_regions = 2 * shell_passes - 1
    u = 1  # Region identification: Fluid between baffles: u=1; Fluid in the U of each HEX: u = -1

    # Properties of the fluid evaluated at the temperature of each node in the previous time step
    cp_shell = np.zeros(n_nodes)
    rho_shell = np.zeros(n_nodes)
    mu_shell = np.zeros(n_nodes)
    k_shell = np.zeros(n_nodes)
    Pr_shell = np.zeros(n_nodes)
    Tw = np.ones(n_nodes)
    Tww = np.ones(n_nodes)

    for i in range(1, n_nodes + 1):
        cp_shell[i - 1] = Inc.specheat(fluid=fluid_shell, T=T_shell[i - 1], P=P_in_shell)
        rho_shell[i - 1] = Inc.density(fluid=fluid_shell, T=T_shell[i - 1], P=P_in_shell)
        mu_shell[i - 1] = Inc.viscosity(fluid=fluid_shell, T=T_shell[i - 1], P=P_in_shell)
        k_shell[i - 1] = 0.5  # Here I impose constant conductivity since we don't have correlations for this fluid
        Pr_shell[i - 1] = cp_shell[i - 1] * mu_shell[i - 1] / k_shell[i - 1]
        Tw[i - 1] = Tenv + 5.0  # Here we suppose the outside wall temperature of the HEX

    # Minimum velocity of the fluid through the shell
    vel_shell[:] = m_shell / rho_shell / (D_shell * gap_baffle)
    # Representative velocity through the shell part.
    vel_shell[:] = 2.0 * vel_shell

    # Parameter and resistances for the first node (k=1) are calculated here:

    # HL_1: Heat transfer from node 1 to node 2 through the baffle
    # Convective Resistance at node 1
    Re_k = rho_shell[0] * L_baffle * vel_shell[0] / mu_shell[0]  # Reynolds at node 1 evaluated at the temperature of node 1 in the previous step time
    Nu_plate_k = 0.0296 * Re_k ** 0.8 * Pr_shell[0] ** 0.333
    h_k = Nu_plate_k * k_shell[0] / L_baffle
    R1conv_i = 1 / A_baffle / h_k
    # Convective Resistance at node 2
    Re_k1 = rho_shell[1] * L_baffle * vel_shell[1] / mu_shell[1]  # Reynolds at node 2 evaluated at the temperature of node 2 in the previous step time
    Nu_plate_k1 = 0.0296 * Re_k1 ** 0.8 * Pr_shell[1] ** 0.333
    h_k1 = Nu_plate_k1 * k_shell[1] / L_baffle
    R1conv_i1 = 1 / A_baffle / h_k1
    # Conductivity resistance through the baffle:
    R1cond_i = th_baffle / A_baffle / k_steel

    # Total resistance for HL_1
    RT_HL1 = R1cond_i + R1conv_i1 + R1conv_i

    # Initialization of some variables of the code:
    Tww[:] = Tw
    error = 1.0

    while error > 0.00001:
        k = 2  # First node calculated in this while loop (1-based index)

        for i in range(1, n_regions + 1):

            if u > 0:  # Baffle region
                vol = Vc
                A_gap = np.pi * D_shell * gap_baffle / 2  # Transfer area from shell to environment
                nodes_region = n_baffles - 1  # Number of nodes in the baffle region
            else:  # U region
                vol = 2 * Vc  # We suppose that the control volume
                A_gap = np.pi * D_shell * gap_baffle  # Transfer area from the shell to the environment for the U-part: I assume it is double the area of the baffles region.
                nodes_region = 1  # Number of nodes in the U region

            for z in range(1, nodes_region + 1):
                if k > n_nodes:
                    break

                RT_HL4 = RT_HL1  # RT_HL4 is the resistance from the node to the previous node through the baffle

                # Heat losses at each node: Resistance calculation

                # HL_1: Heat flux from the node k to node k+1 through the baffle
                Re_k = rho_shell[k - 1] * L_baffle * vel_shell[k - 1] / mu_shell[k - 1]  # Reynolds at node k evaluated with properties in the previous time step
                Nu_plate_k = 0.0296 * Re_k ** 0.8 * Pr_shell[k - 1] ** 0.333
                h_k = Nu_plate_k * k_shell[k - 1] / L_baffle
                R1conv_i = 1 / A_baffle / h_k  # Convective resistance in the baffle side at node k
                if k > n_nodes - 1:
                    R1conv_i = 1000000000.0  # This resistance will be infty in the last node
                else:  # Variables for convective resistance in the side part of the baffle in node k+1
                    Re_k1 = rho_shell[k] * L_baffle * vel_shell[k] / mu_shell[k]
                    Nu_plate_k1 = 0.0296 * Re_k1 ** 0.8 * Pr_shell[k] ** 0.333
                    h_k1 = Nu_plate_k1 * k_shell[k] / L_baffle
                    R1conv_i1 = 1 / A_baffle / h_k1
                R1cond_i = th_baffle / A_baffle / k_steel
                RT_HL1 = R1cond_i + R1conv_i1 + R1conv_i

                # HL_2: Heat losses from shell to environment
                # Air properties are constant:
                mu_air = 184.6e-7
                k_air = 26.3e-3
                rho_air = 1.1614
                visc_air = 15.89e-6
                alpha_air = 22.5e-6
                beta_air = 1 / Tenv
                # Thermal resistances
                R_I = th_ins / A_gap / k_ins

                # Radiative resistance
                emis = 0.9
                sigma = 5.67e-8
                h_rad = emis * sigma * (Tenv ** 2 + Tw[k - 1] ** 2) * (Tenv + Tw[k - 1])
                R_II = (1 / A_gap) * (1 / h_rad)
                # External convective resistance
                Pr_ext = 0.7
                if v_out < 0.1:  # Natural convection outside
                    Ra_ext = beta_air * (Tw[k - 1] - Tenv) * g * D_shell ** 3 / visc_air / alpha_air
                    Nu_ext = (0.60 + 0.387 * Ra_ext ** 0.16 / (1 + (0.559 / Pr_ext) ** (9.0 / 16.0)) ** 0.2963) ** 2
                else:  # Cylinder in cross flow
                    Re_ext = rho_air * v_out * D_shell / mu_air
                    # Equation (7.54) Incropera
                    Nu_ext = (
                        0.3 + (0.62 * Re_ext ** 0.5 * Pr_ext ** 0.333)
                        / (1 + (0.4 / Pr_ext) ** 0.666) ** 0.25
                    ) * (1 + (Re_ext / 282000.0) ** 0.625) ** 0.8
                hext = k_air * Nu_ext / D_shell
                R_III = 1 / (A_gap * hext)
                RT_paralel = (R_II * R_III) / (R_III + R_II)
                RT_HL2 = RT_paralel + R_I  # Thermal resistance between shell and the environment. It depends on Tw, so an iterative process is needed.

                # HL_3: Heat flux from shell to the tubes. This is an input obtained from the profile of temperatures at shell and tube part in the previous time step
                HL_3 = q_HTF_shell[k - 1] * n_tubes

                # Resolution of the energy balance:
                if k > (n_nodes - 1):
                    k = k - 1
                    num = (rho_shell[k - 1] * vol * cp_shell[k - 1] * T_shell[k - 1] / Delta_t
                           + m_shell * cp_shell[k - 2] * T_shell_new[k - 2]
                           - HL_3 + T_shell[k] / RT_HL1 + Tenv / RT_HL2
                           + T_shell_new[k - 2] / RT_HL4)
                    den = (rho_shell[k - 1] * vol * cp_shell[k - 1] / Delta_t
                           + m_shell * cp_shell[k - 1]
                           + 1 / RT_HL1 + 1 / RT_HL2 + 1 / RT_HL4)
                    k = k + 1
                else:
                    num = (rho_shell[k - 1] * vol * cp_shell[k - 1] * T_shell[k - 1] / Delta_t
                           + m_shell * cp_shell[k - 2] * T_shell_new[k - 2]
                           - HL_3 + T_shell[k] / RT_HL1 + Tenv / RT_HL2
                           + T_shell_new[k - 2] / RT_HL4)
                    den = (rho_shell[k - 1] * vol * cp_shell[k - 1] / Delta_t
                           + m_shell * cp_shell[k - 1]
                           + 1 / RT_HL1 + 1 / RT_HL2 + 1 / RT_HL4)

                T_shell_new[k - 1] = num / den

                q_loss = (T_shell_new[k - 1] - Tenv) / RT_HL2
                Tww[k - 1] = q_loss * RT_paralel + Tenv  # Here, wall temperature is recalculated and later, compared with previous wall temperature

                k = k + 1

            u = u * (-1)  # At the end of each region, u value changes and the next region is solved

        # Here wall temperature is compared with the one obtained in the previous iteration
        # error = (Tw-Tww)**2
        # error = max(error) *****VER PORQUÉ NO COINCIDEN TAMAÑOS
        error = 0.000000001
        Tw[:] = Tww
        u = 1

    # T_shell in the last node
    T_shell_new[k - 1] = T_shell_new[k - 2] - (T_shell_new[k - 3] - T_shell_new[k - 2])

    for i in range(1, n_nodes + 1):
        h_shell_i, vel_shell[i - 1] = zukauskas(
            fluid_shell, config, D_shell, S_T, S_L, n_tubes, D_tube,
            T_shell_new[i - 1], m_shell, L_Shell, n_baffles, P_in_shell
        )
        h_shell[i - 1] = h_shell_i

    p_shell_out = P_in_shell  # Constant pressure is assumed in the shell
    # do i = n_baffles+1, n_nodes, n_baffles
    #    call taborek(fluid_shell,D_shell,D_tube,T_shell_new(i),m_shell,S_T,h_shell_i, p_in_shell)
    #    h_shell(i) = h_shell_i
    # end do

    return p_shell_out


# ---------------------------------------------------------------------------
# TUBE 1D INCOMPRESSIBLE SUBROUTINE
# ---------------------------------------------------------------------------

def tube_oneD_inc(q_HTF_tube, T_tube, n_nodes, gap_baffle, d_in,
                  p_tube_in, m_tube, n_tube, fluid_tube, Delta_t,
                  T_tube_new, h_tube, vel_tube):
    """Subroutine tube_oneD_inc: Dynamic solution of the 1D Navier-Stokes equation
    for an incompressible fluid through a tube.

    Modifies T_tube_new, h_tube, vel_tube in-place.

    Parameters
    ----------
    q_HTF_tube  : ndarray[n_nodes] – heat flux from ZeroD_Eq [W]
    T_tube      : ndarray[n_nodes] – tube temperatures from last timestep [K]
    n_nodes     : int              - number of nodes
    gap_baffle  : float            - node/baffle spacing [m]
    d_in        : float            - tube inner diameter [m]
    p_tube_in   : float            - tube-side inlet pressure [Pa]
    m_tube      : float            - tube-side mass-flow rate [kg/s]
    n_tube      : float            - number of tubes [-]
    fluid_tube  : str              - tube-side fluid identifier
    Delta_t     : float            - timestep [s]
    T_tube_new  : ndarray[n_nodes] – (out) updated tube temperatures [K]
    h_tube      : ndarray[n_nodes] – (out) updated convection coefficients [W/m²·K]
    vel_tube    : ndarray[n_nodes] – (out) updated tube velocities [m/s]

    Returns
    -------
    p_tube_out : float – tube-side outlet pressure [Pa]
    """
    J = n_nodes
    L = n_nodes * gap_baffle
    D = d_in
    rh = D / 4.0
    lambda_ = 0.003  # Friction coefficient assumed as constant *****

    # Preparing spatial derivative
    z = np.array([gap_baffle / 2.0 + i * gap_baffle for i in range(J)])
    Dz = z[1] - z[0]

    # Forming derivation matrix D1_D1z (with Finite Differences ( https://web.media.mit.edu/~crtaylor/calculator.html )
    D1_D1z = np.zeros((n_nodes, n_nodes))
    for i in range(3, J - 1):  # i = 3..J-2 (1-based) → range(3, J-1) gives 3,4,...,J-2 then [i-1] for 0-based
        D1_D1z[i - 1, i - 3] = 1.0
        D1_D1z[i - 1, i - 2] = -8.0
        D1_D1z[i - 1, i] = 8.0
        D1_D1z[i - 1, i + 1] = -1.0
    D1_D1z = D1_D1z / 12.0 / Dz

    D1_D1z[0, 0:5] = np.array([-25.0, 48.0, -36.0, 16.0, -3.0]) / 12.0 / Dz   # For the first node
    D1_D1z[1, 0:5] = np.array([-3.0, -10.0, 18.0, -6.0, 1.0]) / 12.0 / Dz     # For the second node
    D1_D1z[J - 1, J - 5:J] = np.array([3.0, -16.0, 36.0, -48.0, 25.0]) / 12.0 / Dz   # For the last node
    D1_D1z[J - 2, J - 5:J] = np.array([-1.0, 6.0, -18.0, 10.0, 3.0]) / 12.0 / Dz    # For the second to last node

    # Steady-State initial:

    # Inlet conditions:
    T_in = T_tube[0] - 273.15  # Own correlations in Celsius
    h_in = (1.181138093142442e-11 * T_in ** 5 - 1.340092005730253e-08 * T_in ** 4
            + 5.872766536462754e-06 * T_in ** 3 + 1.537140422596891e-04 * T_in ** 2
            + 1.624749057644918 * T_in - 19.52) * 1e3  # [J/kg]
    rho_in = 1 / ((5.161496408907718e-19 * (h_in * 1e-3) ** 5
                   - 4.960573299627946e-16 * (h_in * 1e-3) ** 4
                   + 3.134480901443810e-13 * (h_in * 1e-3) ** 3
                   + 7.807498873110475e-11 * (h_in * 1e-3) ** 2
                   + 4.315796357756615e-07 * (h_in * 1e-3)
                   + 9.381445088449658e-04))
    p_in = p_tube_in
    v_in = m_tube / rho_in / (np.pi * D ** 2 / 4.0)

    # Q_HTF BC:
    Q = q_HTF_tube / (Dz * np.pi * D ** 2 / 4.0)

    # Initial condition of the variables
    T_arr = T_tube - 273.15  # Own correlations in Celsius
    h_old = (1.181138093142442e-11 * T_arr ** 5 - 1.340092005730253e-08 * T_arr ** 4
             + 5.872766536462754e-06 * T_arr ** 3 + 1.537140422596891e-04 * T_arr ** 2
             + 1.624749057644918 * T_arr - 19.52) * 1e3  # [kJ/kg]
    rho_old = 1.0 / ((5.161496408907718e-19 * (h_old * 1e-3) ** 5
                      - 4.960573299627946e-16 * (h_old * 1e-3) ** 4
                      + 3.134480901443810e-13 * (h_old * 1e-3) ** 3
                      + 7.807498873110475e-11 * (h_old * 1e-3) ** 2
                      + 4.315796357756615e-07 * (h_old * 1e-3)
                      + 9.381445088449658e-04))

    p_old = np.full(n_nodes, p_tube_in)

    v_old = m_tube / rho_old / (np.pi * D ** 2 / 4.0)

    # Iterative calculation
    lim = 1.0

    h = h_old.copy()
    v = v_old.copy()
    p = p_old.copy()

    while lim > 1e-5:

        h_old = (1.181138093142442e-11 * (T_tube - 273.15) ** 5
                 - 1.340092005730253e-08 * (T_tube - 273.15) ** 4
                 + 5.872766536462754e-06 * (T_tube - 273.15) ** 3
                 + 1.537140422596891e-04 * (T_tube - 273.15) ** 2
                 + 1.624749057644918 * (T_tube - 273.15) - 19.52) * 1e3  # [kJ/kg]
        rho_old = 1.0 / ((5.161496408907718e-19 * (h_old * 1e-3) ** 5
                          - 4.960573299627946e-16 * (h_old * 1e-3) ** 4
                          + 3.134480901443810e-13 * (h_old * 1e-3) ** 3
                          + 7.807498873110475e-11 * (h_old * 1e-3) ** 2
                          + 4.315796357756615e-07 * (h_old * 1e-3)
                          + 9.381445088449658e-04))
        p_old[:] = p_tube_in
        v_old = m_tube / rho_old / (np.pi * D ** 2 / 4.0)

        Dv_Dt = (v - v_old) / Delta_t

        # Newton Raphson for Mass + Momentum + Energy
        dV_dh = ((2.580748204453859e-18 * (h * 1e-3) ** 4
                  - 1.984229319851178e-15 * (h * 1e-3) ** 3
                  + 9.403442704331430e-13 * (h * 1e-3) ** 2
                  + 1.561499774622095e-10 * (h * 1e-3)
                  + 4.315796357756615e-07) * 1e-3)
        V_spec = ((5.161496408907718e-19 * (h * 1e-3) ** 5
                   - 4.960573299627946e-16 * (h * 1e-3) ** 4
                   + 3.134480901443810e-13 * (h * 1e-3) ** 3
                   + 7.807498873110475e-11 * (h * 1e-3) ** 2
                   + 4.315796357756615e-07 * (h * 1e-3)
                   + 9.381445088449658e-04))
        rho = 1 / V_spec
        d_rho_dh = (-rho ** 2) * dV_dh

        # Mass equation
        Dh_Dt = (h - h_old) / Delta_t
        Dh_Dz = D1_D1z @ h
        Dv_Dz = D1_D1z @ v

        F1 = d_rho_dh * Dh_Dt + rho * (D1_D1z @ v) + v * d_rho_dh * Dh_Dz

        J11 = np.diag(rho) @ D1_D1z + diag(d_rho_dh * Dh_Dz)
        J13 = (diag(d_rho_dh) / Delta_t
               + diag(d_rho_dh * Dv_Dz)
               + diag(v * d_rho_dh) @ D1_D1z)

        # Momentum equation
        Dp_Dz = D1_D1z @ p

        F2 = V_spec * Dp_Dz + Dv_Dt + v * Dv_Dz + lambda_ * (v ** 2) / (8.0 * rh)

        J23 = diag(dV_dh * Dp_Dz)

        # J21 = eye(J)/Delta_t + (diag(Dv_Dz) + matmul(diag(v),D1_D1z)) + diag(2*lambda*v/8.d0/rh)
        # Built step by step to save memory
        J21 = np.zeros((J, J))

        # eye(J)/Delta_t:
        for i in range(J):
            J21[i, i] = 1.0 / Delta_t

        # diag(Dv_Dz) + matmul(diag(v), D1_D1z)
        for i in range(J):
            J21[i, i] += Dv_Dz[i] + 2.0 * lambda_ * v[i] / (8.0 * rh)

        # matmul(diag(v), D1_D1z)
        for i in range(J):
            for jj in range(J):
                J21[i, jj] += v[i] * D1_D1z[i, jj]

        J22 = diag(V_spec) @ D1_D1z

        # Energy equation
        Dh_Dz = D1_D1z @ h
        Dp_Dt = (p - p_old) / Delta_t

        F3 = Dh_Dt + v * Dh_Dz + V_spec * (-Q - Dp_Dt)

        J31 = diag(Dh_Dz)
        J32 = -diag(V_spec / Delta_t)

        # J33 = eye(J)/Delta_t + matmul(diag(v),D1_D1z) + diag(dV_dh*(-Q - Dp_Dt))

        # Initialization
        J33 = np.zeros((J, J))

        # eye(J)/Delta_t
        for i in range(J):
            J33[i, i] = 1.0 / Delta_t

        # Operación matmul(diag(v), D1_D1z)
        for i in range(J):
            for jj in range(J):
                J33[i, jj] += v[i] * D1_D1z[i, jj]

        # diag(dV_dh*(-Q - Dp_Dt))
        for i in range(J):
            J33[i, i] += dV_dh[i] * (-Q[i] - Dp_Dt[i])

        # BC
        # For setting velocity BC:
        J11[0, 0] = 1.0
        J11[0, 1:J] = 0.0
        J13[0, 0:J] = 0.0
        F1[0] = v[0] - v_in

        # For setting pressure BC:
        J21[0, 0:J] = 0.0
        J22[0, 0] = 1.0
        J22[0, 1:J] = 0.0
        J23[0, 0:J] = 0.0
        F2[0] = p[0] - p_in

        # For setting enthalpy BC:
        J31[0, 0:J] = 0.0
        J32[0, 0:J] = 0.0
        J33[0, 0] = 1.0
        J33[0, 1:J] = 0.0
        F3[0] = h[0] - h_in

        F = np.empty(3 * n_nodes)
        F[0:J] = F1
        F[J:2 * J] = F2
        F[2 * J:3 * J] = F3

        mat1 = zeros(J, J)

        # Solving linear system (Newton Raphson): **** CHECK JACOB MATRIX!!! *******
        Jacob = np.zeros((3 * n_nodes, 3 * n_nodes))

        Jacob[0:J, 0:J] = J11
        Jacob[0:J, J:2 * J] = mat1
        Jacob[0:J, 2 * J:3 * J] = J13

        Jacob[J:2 * J, 0:J] = J21
        Jacob[J:2 * J, J:2 * J] = J22
        Jacob[J:2 * J, 2 * J:3 * J] = J23

        Jacob[2 * J:3 * J, 0:J] = J31
        Jacob[2 * J:3 * J, J:2 * J] = J32
        Jacob[2 * J:3 * J, 2 * J:3 * J] = J33

        lth = 3 * J

        DX = gauss_elimination(Jacob, F, lth)

        # DX = matmul(Jinv,F)

        v = v - DX[0:J]
        p = p - DX[J:2 * J]
        h = h - DX[2 * J:3 * J]
        DX = DX ** 2
        lim = np.max(DX)

        lim = np.sqrt(lim)

    T_tube_new[:] = ((-1.629510314864758e-13 * (h * 1e-3) ** 5
                      + 3.147528019331875e-10 * (h * 1e-3) ** 4
                      - 1.359541381535188e-07 * (h * 1e-3) ** 3
                      - 1.822021668868121e-04 * (h * 1e-3) ** 2
                      + 0.619357049700913 * (h * 1e-3)
                      + 11.960141056529313) + 273.15)
    # coeff = T_tube(1)/T_tube_new(1)
    # T_tube_new = coeff*T_tube_new
    h_tube[:] = T_tube_new

    for i in range(1, J + 1):
        h_tube_i, vel_tube[i - 1] = gnielinski(
            T_tube_new[i - 1], m_tube, n_tube, d_in, fluid_tube, p[i - 1]
        )
        h_tube[i - 1] = h_tube_i

    p_tube_out = p[J - 1]

    return p_tube_out


# ---------------------------------------------------------------------------
# UTILITY / MATHEMATICAL FUNCTIONS
# ---------------------------------------------------------------------------

def diag(vector):
    """Returns a square matrix with a given vector on its diagonal."""
    n = len(vector)
    matrix = np.zeros((n, n))
    for i in range(n):
        matrix[i, i] = vector[i]
    return matrix


def eye_mat(n):
    """Returns a square matrix with ones in the diagonal and zeros in the rest."""
    matrix = np.zeros((n, n))
    for i in range(n):
        matrix[i, i] = 1.0
    return matrix


def frobenius(a, n):
    """Returns the Frobenius norm of a square matrix."""
    f = 0.0
    for i in range(1, n + 1):
        for j in range(1, n + 1):
            f = f + a[i - 1, j - 1] ** 2
    f = np.sqrt(f)
    return f


def norm2(v):
    """Returns the 2-norm (Euclidean norm) of a vector."""
    norm = 0.0
    for i in range(len(v)):
        norm = norm + v[i] ** 2
    norm = np.sqrt(norm)
    return norm


def norm2_diff(v1, v2, n):
    """Returns the 2-norm of the difference between two vectors."""
    norm = 0.0
    for i in range(1, n + 1):
        norm = norm + (v1[i - 1] - v2[i - 1]) ** 2
    norm = np.sqrt(norm)
    return norm


def op_norm(a, n):
    """Returns the operator norm (spectral norm) of a square matrix using power iteration."""
    max_iter = 50
    # init b vector
    b = np.full(n, np.sqrt(float(n)))  # convert to float

    for i in range(1, max_iter + 1):
        temp_v = a @ b  # multiply a and b
        temp_norm = norm2(temp_v)  # calc norm of product
        for j in range(1, n + 1):
            temp_v[j - 1] = temp_v[j - 1] / temp_norm  # renormalize b

        # difference between b and temp_v is small therefore found eigenvector
        if norm2_diff(b, temp_v, n) <= 0.001:
            break
        else:  # reset b vector
            b = temp_v

    # a @ b = lambda * b therefore
    op = norm2(a @ b)
    return op


def kappa(a, a_inv, n):
    """Returns the condition number of matrix a given its inverse a_inv."""
    k = op_norm(a, n) * op_norm(a_inv, n)
    return k


def update_flow(prev_flow, solved_flow, min_flow, lr):
    """Applies a learning-rate blending update for flow convergence.

    Parameters
    ----------
    prev_flow   : ndarray – previous iteration flow values
    solved_flow : ndarray – newly solved flow values
    min_flow    : float   - minimum flow threshold
    lr          : float   - learning rate

    Returns
    -------
    flow_out : ndarray – updated flow values
    """
    n = len(prev_flow)
    flow_out = np.empty(n)
    for i in range(1, n + 1):
        # if the flow is very low apply learning rate to exponent of flow
        if solved_flow[i - 1] < min_flow:
            temp = (np.log10(abs(prev_flow[i - 1])) * (1.0 - lr)
                    + np.log10(abs(solved_flow[i - 1])) * lr)
            flow_out[i - 1] = np.sign(solved_flow[i - 1]) * 10 ** temp
        else:  # apply learning rate to flow
            flow_out[i - 1] = prev_flow[i - 1] + (solved_flow[i - 1] - prev_flow[i - 1]) * lr
        # if the flow is less than the min_flow, make it the min_flow
        if flow_out[i - 1] < min_flow:
            flow_out[i - 1] = min_flow
    return flow_out


def zeros(m, n):
    """Returns an m×n matrix of zeros (integer type, matching Fortran allocatable integer array)."""
    return np.zeros((m, n), dtype=float)


# ---------------------------------------------------------------------------
# DEBUG OUTPUT UTILITIES
# ---------------------------------------------------------------------------

def find_66_iter(time, unit, iter):
    """This function writes iteration info to a txt file.

    Parameters
    ----------
    time : float – simulation time
    unit : int   - TRNSYS unit number
    iter : int   - iteration number
    """
    with open('find_66.txt', 'a') as f:
        f.write(f"end of timestep {time:16.9E} {unit:5d} {iter:5d}\n")


def find_66_matrix(time, unit, iter, m):
    """This function writes the contents of a double precision matrix to a txt file.

    Parameters
    ----------
    time : float     - simulation time
    unit : int       - TRNSYS unit number
    iter : int       - iteration number
    m    : ndarray   - matrix to write
    """
    dim = m.shape
    with open('find_66.txt', 'a') as f:
        # write info about matrix
        f.write(f"matrix {time:16.9E} {unit:5d} {iter:5d} {dim[0]:5d} {dim[1]:5d}\n")
        # for every row
        for i in range(dim[0]):
            row_str = '| '
            for j in range(dim[1]):
                row_str += f"{m[i, j]:18.9E} "
            row_str += '|'
            f.write(row_str + '\n')


def find_66_outputs(time, unit, iter, diff, arr):
    """This function writes the contents of a double precision array to a txt file.

    Parameters
    ----------
    time : float     - simulation time
    unit : int       - TRNSYS unit number
    iter : int       - iteration number
    diff : float     - convergence metric
    arr  : ndarray   - output array to write
    """
    final_str = f"{arr[0]:16.9E}"
    for i in range(1, len(arr)):
        final_str += f", {arr[i]:16.9E}"

    with open('find_66.txt', 'a') as f:
        f.write(f"set outputs {time:16.9E} {unit:5d} {iter:5d} {diff:16.9E}\n")
        f.write(final_str + '\n')


# ---------------------------------------------------------------------------
# MATRIX INVERSE (Gauss-Jordan elimination)  - Taken from Type 4049
# ---------------------------------------------------------------------------

def matrixinv(a, n):
    """Subroutine to calculate the inverse of a matrix using Gauss-Jordan elimination.
    # Thank you friend: https://www.webpages.uidaho.edu/~gabrielp/ME549-CE546/matrix-inverse.pdf

    The inverse of matrix a (n×n) is returned as b (n×n).
    Input matrix a is modified in place (matching original Fortran intent(inout) behaviour).

    Parameters
    ----------
    a : ndarray (n×n) – input matrix (modified in-place as in Fortran)
    n : int           - matrix dimension

    Returns
    -------
    b : ndarray (n×n) – inverse of a
    """
    # build the identity matrix
    b = np.zeros((n, n))
    for i in range(1, n + 1):
        for j in range(1, n + 1):
            b[i - 1, j - 1] = 0.0
        b[i - 1, i - 1] = 1.0

    for i in range(1, n + 1):  # this is the big loop over all the columns of a(n,n)
        # in case the entry a(i,i) is zero, we need to find a good pivot; this pivot
        # is chosen as the largest value on the column i from a(j,i) with j = 1,n
        big = abs(a[i - 1, i - 1])
        irow = i
        for j in range(i, n + 1):
            if abs(a[j - 1, i - 1]) > big:
                big = abs(a[j - 1, i - 1])
                irow = j
        # interchange lines i with irow for both a() and b() matrices
        if big > abs(a[i - 1, i - 1]):
            for k in range(1, n + 1):
                dum = a[i - 1, k - 1]  # matrix a()
                a[i - 1, k - 1] = a[irow - 1, k - 1]
                a[irow - 1, k - 1] = dum
                dum = b[i - 1, k - 1]  # matrix b()
                b[i - 1, k - 1] = b[irow - 1, k - 1]
                b[irow - 1, k - 1] = dum
        # divide all entries in line i from a(i,j) by the value a(i,i);
        # same operation for the identity matrix
        dum = a[i - 1, i - 1]
        for j in range(1, n + 1):
            a[i - 1, j - 1] = a[i - 1, j - 1] / dum
            b[i - 1, j - 1] = b[i - 1, j - 1] / dum
        # make zero all entries in the column a(j,i); same operation for indent()
        for j in range(i + 1, n + 1):
            dum = a[j - 1, i - 1]
            for k in range(1, n + 1):
                a[j - 1, k - 1] = a[j - 1, k - 1] - dum * a[i - 1, k - 1]
                b[j - 1, k - 1] = b[j - 1, k - 1] - dum * b[i - 1, k - 1]

    # subtract appropriate multiple of row j from row j-1
    for i in range(1, n):
        for j in range(i + 1, n + 1):
            dum = a[i - 1, j - 1]
            for l in range(1, n + 1):
                a[i - 1, l - 1] = a[i - 1, l - 1] - dum * a[j - 1, l - 1]
                b[i - 1, l - 1] = b[i - 1, l - 1] - dum * b[j - 1, l - 1]

    return b


# ---------------------------------------------------------------------------
# GAUSS ELIMINATION  (with partial pivoting)
# ---------------------------------------------------------------------------

def gauss_elimination(A, b, n):
    """Wrapper: calls gauss_elimination_pivot to solve the linear system A·x = b.

    Parameters
    ----------
    A : ndarray (n×n) – coefficient matrix (modified in-place)
    b : ndarray (n,)  - right-hand-side vector (modified in-place)
    n : int           - system dimension

    Returns
    -------
    x : ndarray (n,) – solution vector
    """
    # Llamada a la subrutina para resolver el sistema
    return gauss_elimination_pivot(n, A, b)


def gauss_elimination_pivot(n, A, b):
    """Solver using Gaussian elimination with partial pivoting.
    # Subrutina para resolver el sistema usando eliminación de Gauss con pivoteo parcial

    Parameters
    ----------
    n : int           - system dimension
    A : ndarray (n×n) – coefficient matrix (modified in-place)
    b : ndarray (n,)  - right-hand-side vector (modified in-place)

    Returns
    -------
    x : ndarray (n,) – solution vector
    """
    x = np.zeros(n)

    # Eliminación de Gauss con pivoteo parcial
    for k in range(1, n):
        # Encontrar el pivote máximo en la columna k
        max_row = k
        for i in range(k + 1, n + 1):
            if abs(A[i - 1, k - 1]) > abs(A[max_row - 1, k - 1]):
                max_row = i

        # Intercambiar filas en A y b si es necesario
        if max_row != k:
            A[k - 1, :], A[max_row - 1, :] = A[max_row - 1, :].copy(), A[k - 1, :].copy()
            b[k - 1], b[max_row - 1] = b[max_row - 1], b[k - 1]

        # Proceso de eliminación
        for i in range(k + 1, n + 1):
            factor = A[i - 1, k - 1] / A[k - 1, k - 1]
            A[i - 1, k - 1:n] = A[i - 1, k - 1:n] - factor * A[k - 1, k - 1:n]
            b[i - 1] = b[i - 1] - factor * b[k - 1]

    # Sustitución hacia atrás
    x[n - 1] = b[n - 1] / A[n - 1, n - 1]
    for i in range(n - 1, 0, -1):
        x[i - 1] = (b[i - 1] - np.sum(A[i - 1, i:n] * x[i:n])) / A[i - 1, i - 1]

    return x
