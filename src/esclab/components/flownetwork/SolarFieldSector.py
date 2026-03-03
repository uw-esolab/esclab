"""Solar Field Sector component model (Type 4034)."""

import math
import subprocess
import time

import numpy as np

from esclab.simulate import Component
from esclab.components.esol_properties import Incompressible
from esclab.components.flownetwork.sf_piping_helpers import (
    PressureDrop,
    Row_shadow,
    diams_inlet,
    diams_return,
    dT_dt_inlet,
    dT_dt_return,
    vols_inlet,
    vols_return,
)
from esclab.components.flownetwork.nn_functions import (
    load_NN,
    dt_dtime_NN,
)

Inc = Incompressible()


class SolarFieldSector(Component):
    """
    Object: Solar Field Sector
    Simulation Studio Model: ESOL4034-SolarFieldSector

    Author: Matt Tuman
    Date:    November 27, 2023
    last modified: November 27, 2023
    Ported by: GitHub Copilot, March 02, 2026

    A solar field sector component that models the thermal-hydraulic behavior of a parabolic
    trough solar field sector, including receiver tube state modeling via neural networks,
    RK-4 time integration for loop and header temperatures, and a detailed defocusing scheme.

    Parameters
    ----------
    eta_defocus_1 : float
        Defocus efficiency setpoint 1 [-].
    eta_defocus_2 : float
        Defocus efficiency setpoint 2 [-].
    eta_defocus_3 : float
        Defocus efficiency setpoint 3 [-].
    eta_defocus_4 : float
        Defocus efficiency setpoint 4 [-].
    IAM_a0 : float
        Incidence angle modifier coefficient a0 [-].
    IAM_a1 : float
        Incidence angle modifier coefficient a1 [-].
    IAM_a2 : float
        Incidence angle modifier coefficient a2 [-].
    eta_tracking : float
        Tracking efficiency [-].
    eta_soil : float
        Soiling efficiency [-].
    eta_reflect : float
        Reflectivity efficiency [-].
    SF_avail : float
        Solar field availability [-].
    Distance_SCA : float
        Distance between SCAs [m].
    Ave_focal_length : float
        Average focal length [m].
    W_ap : float
        Aperture width [m].
    row_distance : float
        Row-to-row distance [m].
    L_exp_loop : float
        Expansion loop length [m].
    n_SCA : int
        Number of SCAs per loop [-].
    n_loop : int
        Number of loops per sector [-].
    fluid_ID : float
        HTF fluid identifier [-].
    mc_receiver_mult : float
        Receiver thermal mass multiplier [-].
    mc_header_mult : float
        Header thermal mass multiplier [-].
    T_nominal : float
        Nominal operating temperature [K].
    Time_lim_df : float
        Time limit for defocus mode 4 before full defocus [s].
    sf_label : int
        Sector index (1-based) used to key shared class-level storage [-].
    max_loop : int
        Maximum number of loops across all sectors [-].
    n_sectors : int
        Total number of sectors [-].
    L_tot : float
        Total loop length [m].
    n_nodes_per_loop : int
        Number of thermal nodes per loop (including inlet and outlet) [-].
    T_init_SF : float
        Initial solar field temperature [K].
    T_init_in_header : float
        Initial inlet header temperature [K].
    T_init_return_header : float
        Initial return header temperature [K].
    D_receiver : float
        Receiver inner diameter [m].
    Roughness_pipe : float
        Pipe roughness [m].
    Prob_pristine : float
        Probability of pristine receiver state [-].
    Prob_brokenGlass : float
        Probability of broken-glass state [-].
    Prob_lostVacuum : float
        Probability of lost-vacuum state [-].
    m_dot_std : float
        Standard deviation of loop mass-flow variation [-].
    n_H2 : int
        Number of hydrogen state parameters [-].
    New_Loops : int
        Flag to regenerate loop definition files (1=yes) [-].
    inlet_header_heat_loss : float
        Inlet header heat loss coefficient [-].
    return_header_heat_loss : float
        Return header heat loss coefficient [-].
    GeomFile : str
        Label string for header geometry file.
    NNBase : str
        Label string for neural-network base path.

    Inputs
    ------
    MassFlow : float
        Inlet mass flow rate [kg/s].
    Pressure : float
        Inlet pressure [Pa].
    Temperature : float
        Inlet fluid temperature [K].
    ControlSignal : int
        Control signal [-].
    MassCounter : float
        Incoming mass counter [m^3].
    ANI : float
        Direct normal irradiance [W/m^2].
    T_amb : float
        Ambient temperature [K].
    T_sky : float
        Sky temperature [K].
    Wind : float
        Wind speed [m/s].
    Theta : float
        Incidence angle [rad].
    Phi : float
        Row shadow angle [rad].
    T_tracking : float
        Tracking temperature setpoint [K].

    Outputs
    -------
    MassFlow_out : float
        Outlet mass flow rate [kg/s].
    Pressure_out : float
        Outlet pressure [Pa].
    Temperature_out : float
        Outlet fluid temperature [K].
    T4Ave : float
        Average temperature at sensor 4 location [K].
    Mass_Counter_out : float
        Total HTF mass in sector [kg].
    defocusing_out : float
        Encoded defocusing information [-].
    defocus_groups_1 through defocus_groups_8 : float
        Per-group defocus status encoding [-].
    temperature_groups_1 through temperature_groups_8 : float
        Per-group temperature status encoding [-].
    """

    # Shared sector data: class-level dict keyed by integer sf_label (1-based).
    # Each entry holds the shared numpy arrays for that sector.
    _sectors = {}

    # *** Model Parameters ***
    eta_defocus_1 = Component.Parameter()          # [-]
    eta_defocus_2 = Component.Parameter()          # [-]
    eta_defocus_3 = Component.Parameter()          # [-]
    eta_defocus_4 = Component.Parameter()          # [-]
    IAM_a0 = Component.Parameter()                 # [-]
    IAM_a1 = Component.Parameter()                 # [-]
    IAM_a2 = Component.Parameter()                 # [-]
    eta_tracking = Component.Parameter()           # [-]
    eta_soil = Component.Parameter()               # [-]
    eta_reflect = Component.Parameter()            # [-]
    SF_avail = Component.Parameter()               # [-]
    Distance_SCA = Component.Parameter()           # [m]
    Ave_focal_length = Component.Parameter()       # [m]
    W_ap = Component.Parameter()                   # [m]
    row_distance = Component.Parameter()           # [m]
    L_exp_loop = Component.Parameter()             # [m]
    n_SCA = Component.Parameter()                  # [-]
    n_loop = Component.Parameter()                 # [-]
    fluid_ID = Component.Parameter()               # [-]
    mc_receiver_mult = Component.Parameter()       # [-]
    mc_header_mult = Component.Parameter()         # [-]
    T_nominal = Component.Parameter()              # [K]
    Time_lim_df = Component.Parameter()            # [s]
    sf_label = Component.Parameter()               # [-]
    max_loop = Component.Parameter()               # [-]
    n_sectors = Component.Parameter()              # [-]
    L_tot = Component.Parameter()                  # [m]
    n_nodes_per_loop = Component.Parameter()       # [-]
    T_init_SF = Component.Parameter()              # [K]
    T_init_in_header = Component.Parameter()       # [K]
    T_init_return_header = Component.Parameter()   # [K]
    D_receiver = Component.Parameter()             # [m]
    Roughness_pipe = Component.Parameter()         # [m]
    Prob_pristine = Component.Parameter()          # [-]
    Prob_brokenGlass = Component.Parameter()       # [-]
    Prob_lostVacuum = Component.Parameter()        # [-]
    m_dot_std = Component.Parameter()              # [-]
    n_H2 = Component.Parameter()                   # [-]
    New_Loops = Component.Parameter()              # [-]
    inlet_header_heat_loss = Component.Parameter() # [-]
    return_header_heat_loss = Component.Parameter()# [-]
    # Label parameters (equivalent to getLabel(CurrentUnit, 1) and getLabel(CurrentUnit, 2))
    GeomFile = Component.Parameter()               # header geometry file label
    NNBase = Component.Parameter()                 # neural network base path label
    # TODO-NEEDS CONVERSION REVIEW: H2_std and H2_mean are variable-length arrays (n_H2 elements each)
    # read from parameter positions 42..42+n_H2-1 (H2_std) and 42+n_H2..42+2*n_H2-1 (H2_mean).
    # In Python these are stored as instance lists in initialize() rather than as Component.Parameter() members.

    # *** Model Inputs ***
    MassFlow = Component.Input()                   # [kg/s]
    Pressure = Component.Input()                   # [Pa]
    Temperature = Component.Input()                # [K]
    ControlSignal = Component.Input()              # [-]
    MassCounter = Component.Input()                # [m^3]
    ANI = Component.Input()                        # [W/m^2]
    T_amb = Component.Input()                      # [K]
    T_sky = Component.Input()                      # [K]
    Wind = Component.Input()                       # [m/s]
    Theta = Component.Input()                      # [rad]
    Phi = Component.Input()                        # [rad]
    T_tracking = Component.Input()                 # [K]

    # *** Model Outputs ***
    MassFlow_out = Component.Output()              # [kg/s]
    Pressure_out = Component.Output()              # [Pa]
    Temperature_out = Component.Output()           # [K]
    T4Ave = Component.Output()                     # [K]
    Mass_Counter_out = Component.Output()          # [kg]
    defocusing_out = Component.Output()            # [-]
    defocus_groups_1 = Component.Output()          # [-]
    defocus_groups_2 = Component.Output()          # [-]
    defocus_groups_3 = Component.Output()          # [-]
    defocus_groups_4 = Component.Output()          # [-]
    defocus_groups_5 = Component.Output()          # [-]
    defocus_groups_6 = Component.Output()          # [-]
    defocus_groups_7 = Component.Output()          # [-]
    defocus_groups_8 = Component.Output()          # [-]
    temperature_groups_1 = Component.Output()      # [-]
    temperature_groups_2 = Component.Output()      # [-]
    temperature_groups_3 = Component.Output()      # [-]
    temperature_groups_4 = Component.Output()      # [-]
    temperature_groups_5 = Component.Output()      # [-]
    temperature_groups_6 = Component.Output()      # [-]
    temperature_groups_7 = Component.Output()      # [-]
    temperature_groups_8 = Component.Output()      # [-]

    def initialize(self):
        """One-time setup equivalent to getIsFirstCallofSimulation() + getIsStartTime()."""
        # Clear any previously allocated shared sector data (equivalent to deallocate_memory_sf if already allocated)
        if SolarFieldSector._sectors:
            SolarFieldSector._sectors.clear()
            # deallocate_memory_sf from SF_data

        # getIsStartTime() block: read parameters and perform first-timestep initialization
        sf_label = int(self.sf_label.v)
        n_loop = int(self.n_loop.v)
        n_SCA = int(self.n_SCA.v)
        n_nodes_per_loop = int(self.n_nodes_per_loop.v)
        n_h2 = int(self.n_H2.v)
        max_loop = int(self.max_loop.v)
        n_sectors = int(self.n_sectors.v)

        # Read H2_std and H2_mean: variable-length parameters beyond index 41
        # TODO-NEEDS CONVERSION REVIEW: H2_std and H2_mean are provided as variable-length parameter
        # lists. In Python they are collected here as instance attributes.
        H2_std = list(getattr(self, f'_h2_std_raw', []))
        H2_mean = list(getattr(self, f'_h2_mean_raw', []))

        # Generate Loop Data Using Python Scripts (if desired)
        # TODO-NEEDS CONVERSION REVIEW: 0-based indexing — sf_label is 1-based in Fortran; used here
        # to build sector directory index string.
        num_keys = ['0', '1', '2', '3', '4', '5', '6', '7', '8', '9']
        sector_ind = num_keys[sf_label]  # TODO-NEEDS CONVERSION REVIEW: 0-based indexing — Fortran uses sf_label+1 on 1-based num_keys
        loop_base = 'LoopDef_Sector' + sector_ind

        if int(self.New_Loops.v) == 1:
            # Check if loop file directory exists and clear it
            subprocess.run(
                f'If exist {loop_base.strip()} del /f /s /q {loop_base.strip()} 1>nul',
                shell=True,
            )
            subprocess.run(
                f'If exist {loop_base.strip()} rmdir /s /q {loop_base.strip()}',
                shell=True,
            )
            subprocess.run(f'mkdir {loop_base.strip()}', shell=True)

            # Build python command string for loop generation
            python_command = f'python gen_loopFile.py {int(self.L_tot.v):5d} {n_nodes_per_loop:3d}'
            python_command += f' {self.Prob_pristine.v:4.2f}'
            python_command += f' {self.Prob_brokenGlass.v:4.2f}'
            python_command += f' {self.Prob_lostVacuum.v:4.2f}'
            python_command += f' {self.m_dot_std.v:8.4f}'
            python_command += f' {n_loop:5d}'
            python_command += f' {sf_label:5d}'
            python_command += f' {n_h2:5d}'
            for n in range(n_h2):
                python_command += f' {H2_std[n]:8.4f}'
            for n in range(n_h2):
                python_command += f' {H2_mean[n]:8.4f}'
            python_command += f' >{loop_base}/output.log'

            # Call python script to generate loop
            subprocess.run(python_command, shell=True)

        # Initialize neural network data and solar field global variables
        # if this is the first SF sector
        if sf_label == 1:
            # Directly allocate zeros here, replace allocate_memory_sf from SF_data
            # Initializes shared arrays for all sectors. In Python we use the class-level _sectors dict.
            for s in range(1, n_sectors + 1):
                SolarFieldSector._sectors[s] = {
                    # Per-sector shared arrays (Fortran: indexed by sf_label)
                    't_sf': np.zeros((n_nodes_per_loop, max_loop)),           # t_sf(n_nodes_per_loop, max_loop, s)
                    't_header_inlet': np.zeros(max_loop),                     # t_header_inlet(:, s)
                    't_header_return': np.zeros(max_loop + 2),                # t_header_return(:, s)
                    'defocus_mode': np.ones(max_loop),                        # defocus_mode(:, s)
                    'time_df': np.zeros((max_loop,)),                         # time_df(loop, s)
                    'mass_HTF_hold': 0.0,                                     # mass_HTF_hold(s)
                    't4Ave_hold': 0.0,                                        # t4Ave_hold(s)
                    't_bar_sf': np.zeros(n_nodes_per_loop),                   # t_bar_sf(:, s)
                    'num_cv_header': 0,                                       # num_cv_header(s)
                    'Vol_inlet': np.zeros(max_loop),                          # Vol_inlet(:, s)
                    'Vol_return': np.zeros(max_loop + 1),                     # Vol_return(:, s)
                    'D_inlet': np.zeros(max_loop),                            # D_inlet(:, s)
                    'D_return': np.zeros(max_loop + 1),                       # D_return(:, s)
                    'L_cv_inlet': np.zeros(max_loop),                         # L_cv_inlet(:, s)
                    'L_cv_return': np.zeros(max_loop + 1),                    # L_cv_return(:, s)
                    'inds_header_in': np.zeros(max_loop, dtype=int),          # inds_header_in(:, s)
                    'inds_left': np.zeros(max_loop // 2 + 1, dtype=int),      # inds_left(:, s)
                    'inds_right': np.zeros(max_loop // 2 + 1, dtype=int),     # inds_right(:, s)
                    'm_dot_var': np.ones(max_loop),                           # m_dot_var(loop, s)
                    'nCV_state': np.zeros((4, max_loop), dtype=int),          # nCV_state(:, loop, s)
                    'inds_pristine': np.zeros((n_nodes_per_loop - 1, max_loop), dtype=int),  # inds_pristine(:, loop, s)
                    'inds_lVacuum': np.zeros((n_nodes_per_loop - 1, max_loop), dtype=int),   # inds_lVacuum(:, loop, s)
                    'inds_bGlass': np.zeros((n_nodes_per_loop - 1, max_loop), dtype=int),    # inds_bGlass(:, loop, s)
                    'inds_H2': np.zeros((n_nodes_per_loop - 1, max_loop), dtype=int),        # inds_H2(:, loop, s)
                    'H2_pressure': np.zeros((n_nodes_per_loop - 1, max_loop)),               # H2_pressure(:, loop, s)
                    'r_number': np.zeros((n_SCA, max_loop)),                                 # r_number(:, loop, s)
                    'features': np.zeros((n_nodes_per_loop - 1, 7)),                         # features(:, :) per-sector workspace
                    'minMax2': np.zeros((16, 4)),                                            # minMax2(16,4) — shared NN scaling
                    'temperature_groups': np.zeros(8),                                       # temperature_groups(:, s)
                    # RK4 arrays for return header
                    'k1_rh': np.zeros(max_loop + 2),
                    'k2_rh': np.zeros(max_loop + 2),
                    'k3_rh': np.zeros(max_loop + 2),
                    'k4_rh': np.zeros(max_loop + 2),
                    # RK4 arrays for inlet header
                    'k1_ih': np.zeros(max_loop),
                    'k2_ih': np.zeros(max_loop),
                    'k3_ih': np.zeros(max_loop),
                    'k4_ih': np.zeros(max_loop),
                    # RK4 arrays for loops
                    'k1': np.zeros(n_nodes_per_loop),
                    'k2': np.zeros(n_nodes_per_loop),
                    'k3': np.zeros(n_nodes_per_loop),
                    'k4': np.zeros(n_nodes_per_loop),
                    # Intermediate arrays
                    't_hat': np.zeros(n_nodes_per_loop),
                    't_bar_hat': np.zeros(n_nodes_per_loop - 1),
                    't_hat_return': np.zeros(max_loop + 2),
                    't_bar_hat_return': np.zeros(max_loop + 1),
                    't_hat_inlet': np.zeros(max_loop),
                    't_bar_hat_inlet': np.zeros(max_loop - 1),
                    't_bar_inlet': np.zeros(max_loop - 1),
                    't_bar_return': np.zeros(max_loop + 1),
                    't_bar': np.zeros(n_nodes_per_loop - 1),
                    # Scalar shared state
                    'L_segment': 0.0,
                    'Vol': 0.0,
                    'defocus_groups': np.zeros(8),   # defocus_groups(8) — persisted from end-of-timestep block
                    'dni_array': np.zeros(n_nodes_per_loop - 1),
                    'dni_labels': np.zeros(n_nodes_per_loop - 1),
                    'm_dots_return': np.zeros(max_loop + 1),
                    'm_dots_in': np.zeros(max_loop),
                    'm_left': np.zeros(max_loop // 2 + 1),
                    'm_right': np.zeros(max_loop // 2 + 1),
                    't_hold_l': np.zeros(max_loop // 2 + 1),
                    't_hold_r': np.zeros(max_loop // 2 + 1),
                }

            # Compute L_segment and Vol (per Fortran: set once when sf_label==1)
            SolarFieldSector._sectors[1]['L_segment'] = self.L_tot.v / (n_nodes_per_loop - 1)
            L_segment = SolarFieldSector._sectors[1]['L_segment']
            Vol = L_segment * 3.1415 * (self.D_receiver.v / 2.0) ** 2
            SolarFieldSector._sectors[1]['Vol'] = Vol
            # Propagate L_segment and Vol to all sectors (shared global values in Fortran)
            for s in range(1, n_sectors + 1):
                SolarFieldSector._sectors[s]['L_segment'] = L_segment
                SolarFieldSector._sectors[s]['Vol'] = Vol

            # Random_Number from SF_data — populate r_number arrays
            for s in range(1, n_sectors + 1):
                SolarFieldSector._sectors[s]['r_number'] = np.random.random(
                    (n_SCA, max_loop)
                )

            # Load neural network weights and populate minMax2 scaling arrays.
            # Fortran: minMax2 = load_NN_toMod(NNbase) inside allocate_memory_sf.
            # In Python, load_NN returns the (16,4) minMax matrix; store it in
            # every sector's 'minMax2' entry (it is the same for all sectors).
            minMax2 = load_NN(self.NNBase.v)
            for s in range(1, n_sectors + 1):
                SolarFieldSector._sectors[s]['minMax2'] = minMax2

        # Worker arrays stored per instance (no sf_label dimension in Fortran — used only within a timestep)
        self._curr_inds = np.zeros(n_nodes_per_loop - 1, dtype=int)
        self._t_sf_hold = np.zeros(n_nodes_per_loop)
        self._inds_pristine_hold = np.zeros(n_nodes_per_loop - 1, dtype=int)
        self._inds_lVacuum_hold = np.zeros(n_nodes_per_loop - 1, dtype=int)
        self._inds_bGlass_hold = np.zeros(n_nodes_per_loop - 1, dtype=int)
        self._inds_H2_hold = np.zeros(n_nodes_per_loop - 1, dtype=int)
        self._nCV_state_hold = np.zeros(4, dtype=int)

        sec = SolarFieldSector._sectors[sf_label]

        # Determine control volumes of inlet and return headers
        cc = 0
        # TODO-NEEDS CONVERSION REVIEW: 0-based indexing — Fortran loop k from 1 to n_loop/2
        for k in range(1, n_loop // 2 + 1):
            if k == 1:
                cc += 1
            elif k % 2 == 0:
                cc += 1
            else:
                cc += 2
        n_node_header = cc
        n_cv_header = cc - 1
        sec['num_cv_header'] = n_cv_header

        sec['Vol_inlet'][0:n_cv_header] = vols_inlet(n_cv_header, n_loop, self.row_distance.v, self.L_exp_loop.v, self.GeomFile.v)
        sec['Vol_return'][0:n_cv_header+1] = vols_return(n_cv_header+1, n_loop, self.row_distance.v, self.L_exp_loop.v, self.GeomFile.v)
        sec['D_inlet'][0:n_cv_header] = diams_inlet(n_cv_header, n_loop, self.row_distance.v, self.L_exp_loop.v, self.GeomFile.v)
        sec['D_return'][0:n_cv_header+1] = diams_return(n_cv_header+1, n_loop, self.row_distance.v, self.L_exp_loop.v, self.GeomFile.v)

        # Compute control volume lengths of headers
        # TODO-NEEDS CONVERSION REVIEW: 0-based indexing — Fortran n from 1 to n_cv_header (1-based); Python uses n from 0 to n_cv_header-1
        cc = 1
        for n in range(n_cv_header):
            if cc > 1:
                sec['L_cv_return'][n] = (self.row_distance.v * 2.0 + self.L_exp_loop.v * 2.0) / 2.0
                sec['L_cv_inlet'][n] = (self.row_distance.v * 2.0 + self.L_exp_loop.v * 2.0) / 2.0
                cc += 1
                if cc == 4:
                    cc = 1
            else:
                sec['L_cv_return'][n] = self.row_distance.v * 2.0
                sec['L_cv_inlet'][n] = self.row_distance.v * 2.0
                cc += 1
        # Last return header CV length (Fortran: L_cv_return(n, sf_label) where n = n_cv_header+1 after loop exit is n_cv_header+1, but loop runs n from 1 to n_cv_header, so n after loop = n_cv_header)
        # TODO-NEEDS CONVERSION REVIEW: 0-based indexing — Fortran writes L_cv_return(n, sf_label) after loop where n==n_cv_header (1-based), so Python index is n_cv_header (0-based = n_cv_header)
        sec['L_cv_return'][n_cv_header] = self.row_distance.v * 2.0

        # Compute indices where mass flow is leaving inlet header
        # TODO-NEEDS CONVERSION REVIEW: 0-based indexing — Fortran arrays are 1-based; inds_header_in values stored here
        # are node indices (1-based in Fortran). Converting to 0-based here.
        cc = 1
        k = 0
        for n in range(1, n_node_header + 1):
            if cc == 1:
                sec['inds_header_in'][k] = n - 1  # TODO-NEEDS CONVERSION REVIEW: 0-based indexing
                cc += 1
                k += 1
            elif cc == 2:
                sec['inds_header_in'][k] = n - 1  # TODO-NEEDS CONVERSION REVIEW: 0-based indexing
                cc += 1
                k += 1
            else:
                cc = 1

        # Determine index mapping (solar field loop to return header)
        # TODO-NEEDS CONVERSION REVIEW: 0-based indexing — loop n from n_loop down to 1 (1-based);
        # stored values are loop indices converted to 0-based
        ll = 0
        rr = 0
        for n in range(n_loop, 0, -1):
            if n % 2 == 0:
                sec['inds_right'][rr] = n - 1  # TODO-NEEDS CONVERSION REVIEW: 0-based indexing
                rr += 1
            else:
                sec['inds_left'][ll] = n - 1   # TODO-NEEDS CONVERSION REVIEW: 0-based indexing
                ll += 1

        # Load in receiver states for each loop from loop definition file
        loop_file = 'LoopDef_Sector' + sector_ind + '/loops.txt'
        # TODO-NEEDS CONVERSION REVIEW: 0-based indexing — Fortran loops nn from 1 to n_loop; Python uses 0 to n_loop-1
        with open(loop_file, 'r') as fh:
            for nn in range(n_loop):
                fh.readline()  # skip Loop Number
                fh.readline()  # skip Label
                # total nodes and length of loop
                line = fh.readline().strip().split()
                n_node = int(line[0])
                L_tot_read = float(line[1])
                n_cv = n_node - 1
                fh.readline()  # skip blank/separator
                # Pristine Mapping
                sec['nCV_state'][0, nn] = int(fh.readline().strip())
                if sec['nCV_state'][0, nn] > 0:
                    vals = list(map(int, fh.readline().strip().split()))
                    # TODO-NEEDS CONVERSION REVIEW: 0-based indexing — Fortran indices are 1-based receiver CV indices
                    sec['inds_pristine'][:sec['nCV_state'][0, nn], nn] = [v - 1 for v in vals]
                fh.readline()  # skip separator
                # Lost Vacuum Mapping
                sec['nCV_state'][1, nn] = int(fh.readline().strip())
                if sec['nCV_state'][1, nn] > 0:
                    vals = list(map(int, fh.readline().strip().split()))
                    # TODO-NEEDS CONVERSION REVIEW: 0-based indexing
                    sec['inds_lVacuum'][:sec['nCV_state'][1, nn], nn] = [v - 1 for v in vals]
                fh.readline()  # skip separator
                # Broken Glass Mapping
                sec['nCV_state'][2, nn] = int(fh.readline().strip())
                if sec['nCV_state'][2, nn] > 0:
                    vals = list(map(int, fh.readline().strip().split()))
                    # TODO-NEEDS CONVERSION REVIEW: 0-based indexing
                    sec['inds_bGlass'][:sec['nCV_state'][2, nn], nn] = [v - 1 for v in vals]
                fh.readline()  # skip separator
                # Hydrogen Mapping
                sec['nCV_state'][3, nn] = int(fh.readline().strip())
                if sec['nCV_state'][3, nn] > 0:
                    vals = list(map(int, fh.readline().strip().split()))
                    # TODO-NEEDS CONVERSION REVIEW: 0-based indexing
                    sec['inds_H2'][:sec['nCV_state'][3, nn], nn] = [v - 1 for v in vals]
                fh.readline()  # skip separator
                # Read H2 pressures in collectors
                h2_count = sec['nCV_state'][3, nn]
                if h2_count > 0:
                    vals = list(map(float, fh.readline().strip().split()))
                    sec['H2_pressure'][:h2_count, nn] = vals
                else:
                    fh.readline()
                # Read Mass flow variation
                fh.readline()  # skip label
                sec['m_dot_var'][nn] = float(fh.readline().strip())

        # Define SF_Avail mask for each sector
        # TODO-NEEDS CONVERSION REVIEW: 0-based indexing — Fortran loop, n from 1 to n_SCA and 1 to n_loop
        for loop_idx in range(n_loop):
            for n in range(n_SCA):
                if sec['r_number'][n, loop_idx] <= self.SF_avail.v:
                    sec['r_number'][n, loop_idx] = 1.0
                else:
                    sec['r_number'][n, loop_idx] = 0.0

        # Initialize temperatures in solar field
        diff_T = self.T_init_return_header.v - self.T_init_in_header.v
        dT = diff_T / n_nodes_per_loop
        T_sf_lin = self.T_init_in_header.v
        # TODO-NEEDS CONVERSION REVIEW: 0-based indexing — Fortran n from 1 to n_nodes_per_loop
        for n in range(n_nodes_per_loop):
            sec['t_sf'][n, :] = T_sf_lin
            T_sf_lin += dT

        sec['t_sf'][:, :] = self.T_init_SF.v
        sec['t_hat'][:] = self.T_init_SF.v
        sec['t_bar'][:] = self.T_init_SF.v
        sec['t_bar_hat'][:] = sec['t_bar'][:]
        sec['t_bar_inlet'][:] = self.T_init_in_header.v
        sec['t_header_inlet'][:] = self.T_init_in_header.v
        sec['t_bar_return'][:] = self.T_init_return_header.v
        sec['t_header_return'][:] = self.T_init_return_header.v

        # Initialize Tracking Modes
        sec['defocus_mode'][:] = 1.0

        # HTF Mass in SF Computation
        mass_HTF_hold = 0.0

        # Inlet Header
        # TODO-NEEDS CONVERSION REVIEW: 0-based indexing — Fortran n from 1 to n_cv_header
        for n in range(n_cv_header):
            mass_HTF_hold += sec['Vol_inlet'][n] * Inc.density(
                fluid=self.fluid_ID.v, T=self.T_init_in_header.v, P=0.0
            )

        # Solar Field
        Vol = sec['Vol']
        # TODO-NEEDS CONVERSION REVIEW: 0-based indexing — Fortran loop from 1 to n_loop, n from 1 to n_nodes_per_loop-1
        for loop_idx in range(n_loop):
            for n in range(n_nodes_per_loop - 1):
                mass_HTF_hold += Vol * Inc.density(
                    fluid=self.fluid_ID.v, T=sec['t_sf'][n, loop_idx], P=0.0
                )

        # Return Header
        # TODO-NEEDS CONVERSION REVIEW: 0-based indexing — Fortran n from 1 to n_cv_header+1
        for n in range(n_cv_header + 1):
            mass_HTF_hold += sec['Vol_return'][n] * Inc.density(
                fluid=self.fluid_ID.v, T=self.T_init_return_header.v, P=0.0
            )

        sec['mass_HTF_hold'] = mass_HTF_hold

        # Set the Initial Values of the Outputs
        self.MassFlow_out.v = self.MassFlow.v                       # MassFlow
        self.Pressure_out.v = self.Pressure.v                       # Pressure
        self.Temperature_out.v = self.T_init_return_header.v        # Temperature
        self.T4Ave.v = self.T_init_SF.v                             # T4Ave
        self.Mass_Counter_out.v = sec['mass_HTF_hold']              # Mass Counter

        # Define Relevant Tracking Variables
        self._eta_defocus = np.array([
            self.eta_defocus_1.v, self.eta_defocus_2.v,
            self.eta_defocus_3.v, self.eta_defocus_4.v, 0.0,
        ])
        self._T_41 = self.T_tracking.v + 2.0   # [K]
        self._T_42 = self.T_tracking.v + 4.0   # [K]
        self._T_4d = self.T_tracking.v + 6.0   # [K]
        self._eta_tot = self.eta_tracking.v * self.eta_soil.v * self.eta_reflect.v

    def calculate(self):
        super().calculate()

        # getIsReReadParameters() equivalent: recompute derived scalar(s) each call
        self._eta_tot = self.eta_tracking.v * self.eta_soil.v * self.eta_reflect.v

        # Read inputs
        # (values already updated on self.XXX.v via connection mechanism before calculate() is called)

        sf_label = int(self.sf_label.v)
        sec = SolarFieldSector._sectors[sf_label]

        DTheta = self.model.settings.timestep * 3600.0  # [s]; settings.timestep is in hours

        n_loop = int(self.n_loop.v)
        n_SCA = int(self.n_SCA.v)
        n_nodes_per_loop = int(self.n_nodes_per_loop.v)
        n_cv_header = sec['num_cv_header']
        n_node_header = n_cv_header + 1

        # Perform Thermal Computations at the End of Each Timestep
        if self.model.is_converged:
            start = time.perf_counter()

            # Define Relevant Tracking Variables
            eta_defocus = np.array([
                self.eta_defocus_1.v, self.eta_defocus_2.v,
                self.eta_defocus_3.v, self.eta_defocus_4.v, 0.0,
            ])
            T_41 = self.T_tracking.v + 2.0   # [K]
            T_42 = self.T_tracking.v + 4.0   # [K]
            T_4d = self.T_tracking.v + 6.0   # [K]
            eta_tot = self.eta_tracking.v * self.eta_soil.v * self.eta_reflect.v

            # Specify individual mass flow
            m_dot = self.MassFlow.v / n_loop

            # Update inlet header temperatures and specify inlet loop temperatures
            # TODO-NEEDS CONVERSION REVIEW: 0-based indexing — t_header_inlet and t_sf are 0-based
            sec['t_header_inlet'][0] = self.Temperature.v
            cc = 0
            # TODO-NEEDS CONVERSION REVIEW: 0-based indexing — Fortran n from 1 to n_loop/2, cc starting at 1 (0-based: cc=0 here)
            for n in range(n_loop // 2):
                sec['t_sf'][0, cc] = sec['t_header_inlet'][sec['inds_header_in'][n]]
                sec['t_sf'][0, cc + 1] = sec['t_header_inlet'][sec['inds_header_in'][n]]
                cc += 2

            # Return Header Calculations
            # Update Return Header Node 1 Temperature
            # TODO-NEEDS CONVERSION REVIEW: 0-based indexing — Fortran t_sf(n_node, n_loop, sf), t_sf(n_node, n_loop-1, sf)
            # where n_node is 1-based last node. In Python: t_sf[n_nodes_per_loop-1, n_loop-1] etc.
            sec['t_header_return'][0] = (
                sec['t_sf'][n_nodes_per_loop - 1, n_loop - 1] * sec['m_dot_var'][n_loop - 1]
                + sec['t_sf'][n_nodes_per_loop - 1, n_loop - 2] * sec['m_dot_var'][n_loop - 2]
            ) / (sec['m_dot_var'][n_loop - 1] + sec['m_dot_var'][n_loop - 2])

            # Specify mass flow rates out of each return header control volume
            # TODO-NEEDS CONVERSION REVIEW: 0-based indexing — Fortran n from 1 to n_cv_header+1;
            # inds_header_in values are 0-based (converted in initialize())
            cc = 0
            jj = 1  # TODO-NEEDS CONVERSION REVIEW: 0-based indexing — Fortran jj starts at 2 (1-based); here jj=1 indexes the 2nd entry (0-based)
            for n in range(n_cv_header + 1):
                if n == 0:
                    sec['m_dots_return'][n] = m_dot * (
                        sec['m_dot_var'][sec['inds_right'][cc]]
                        + sec['m_dot_var'][sec['inds_left'][cc]]
                    )
                    cc += 1
                else:
                    if n == sec['inds_header_in'][jj]:  # if CV n has SF loops attached
                        sec['m_dots_return'][n] = (
                            sec['m_dots_return'][n - 1]
                            + m_dot * (
                                sec['m_dot_var'][sec['inds_right'][cc]]
                                + sec['m_dot_var'][sec['inds_left'][cc]]
                            )
                        )
                        jj += 1
                        cc += 1
                    else:
                        sec['m_dots_return'][n] = sec['m_dots_return'][n - 1]

            # Step Return Header Temperatures Through Time (RK4)
            # Update temporary arrays to avoid warning in function call
            # TODO-NEEDS CONVERSION REVIEW: 0-based indexing — Fortran 1:n_loop/2 slice → Python 0:n_loop//2
            m_left = m_dot * sec['m_dot_var'][sec['inds_left'][:n_loop // 2]]
            m_right = m_dot * sec['m_dot_var'][sec['inds_right'][:n_loop // 2]]
            t_hold_l = sec['t_sf'][n_nodes_per_loop - 1, sec['inds_left'][:n_loop // 2]]
            t_hold_r = sec['t_sf'][n_nodes_per_loop - 1, sec['inds_right'][:n_loop // 2]]
            t_hold = sec['t_header_return'][:n_node_header + 1].copy()
            vol_hold = sec['Vol_return'][:n_cv_header + 1].copy()
            L_cv_hold = sec['L_cv_return'][:n_cv_header + 1].copy()
            inds_hold = sec['inds_header_in'][:n_loop // 2].copy()

            # Compute control volume average temperatures
            # TODO-NEEDS CONVERSION REVIEW: 0-based indexing — Fortran t_hold(2:n_node_header+1) → Python t_hold[1:n_node_header+1]
            sec['t_bar_return'][:n_node_header] = (
                t_hold[1:n_node_header + 1] + t_hold[:n_node_header]
            ) / 2.0

            # K1
            sec['k1_rh'][:n_node_header + 1] = dT_dt_return(m_left[:n_loop // 2], m_right[:n_loop // 2], sec['m_dots_return'][:n_node_header], t_hold_l[:n_loop // 2], t_hold_r[:n_loop // 2], t_hold[:n_node_header + 1], sec['t_bar_return'][:n_node_header], vol_hold[:n_cv_header + 1], L_cv_hold[:n_cv_header + 1], self.mc_header_mult.v, n_node_header + 1, n_loop, inds_hold[:n_loop // 2], self.fluid_ID.v, self.return_header_heat_loss.v)
            sec['t_hat_return'][:n_node_header + 1] = (
                t_hold[:n_node_header + 1] + sec['k1_rh'][:n_node_header + 1] * DTheta / 2.0
            )
            sec['t_bar_hat_return'][:n_node_header] = (
                sec['t_hat_return'][1:n_node_header + 1]
                + sec['t_hat_return'][:n_node_header]
            ) / 2.0

            # K2
            sec['k2_rh'][:n_node_header + 1] = dT_dt_return(m_left[:n_loop // 2], m_right[:n_loop // 2], sec['m_dots_return'][:n_node_header], t_hold_l[:n_loop // 2], t_hold_r[:n_loop // 2], sec['t_hat_return'][:n_node_header + 1], sec['t_bar_hat_return'][:n_node_header], vol_hold[:n_cv_header + 1], L_cv_hold[:n_cv_header + 1], self.mc_header_mult.v, n_node_header + 1, n_loop, inds_hold[:n_loop // 2], self.fluid_ID.v, self.return_header_heat_loss.v)
            sec['t_hat_return'][:n_node_header + 1] = (
                t_hold[:n_node_header + 1] + sec['k2_rh'][:n_node_header + 1] * DTheta / 2.0
            )
            sec['t_bar_hat_return'][:n_node_header] = (
                sec['t_hat_return'][1:n_node_header + 1]
                + sec['t_hat_return'][:n_node_header]
            ) / 2.0

            # K3
            sec['k3_rh'][:n_node_header + 1] = dT_dt_return(m_left[:n_loop // 2], m_right[:n_loop // 2], sec['m_dots_return'][:n_node_header], t_hold_l[:n_loop // 2], t_hold_r[:n_loop // 2], sec['t_hat_return'][:n_node_header + 1], sec['t_bar_hat_return'][:n_node_header], vol_hold[:n_cv_header + 1], L_cv_hold[:n_cv_header + 1], self.mc_header_mult.v, n_node_header + 1, n_loop, inds_hold[:n_loop // 2], self.fluid_ID.v, self.return_header_heat_loss.v)
            sec['t_hat_return'][:n_node_header + 1] = (
                t_hold[:n_node_header + 1] + sec['k3_rh'][:n_node_header + 1] * DTheta
            )
            sec['t_bar_hat_return'][:n_node_header] = (
                sec['t_hat_return'][1:n_node_header + 1]
                + sec['t_hat_return'][:n_node_header]
            ) / 2.0

            # K4
            sec['k4_rh'][:n_node_header + 1] = dT_dt_return(m_left[:n_loop // 2], m_right[:n_loop // 2], sec['m_dots_return'][:n_node_header], t_hold_l[:n_loop // 2], t_hold_r[:n_loop // 2], sec['t_hat_return'][:n_node_header + 1], sec['t_bar_hat_return'][:n_node_header], vol_hold[:n_cv_header + 1], L_cv_hold[:n_cv_header + 1], self.mc_header_mult.v, n_node_header + 1, n_loop, inds_hold[:n_loop // 2], self.fluid_ID.v, self.return_header_heat_loss.v)

            sec['t_header_return'][:n_node_header + 1] = t_hold[:n_node_header + 1] + (
                sec['k1_rh'][:n_node_header + 1] / 6.0
                + sec['k2_rh'][:n_node_header + 1] / 3.0
                + sec['k3_rh'][:n_node_header + 1] / 3.0
                + sec['k4_rh'][:n_node_header + 1] / 6.0
            ) * DTheta

            # Inlet Header Calculations
            # Specify mass flow rates through each inlet header control volume
            # TODO-NEEDS CONVERSION REVIEW: 0-based indexing — Fortran n from 1 to n_cv_header, jj starts at 1
            jj = 0
            for n in range(n_cv_header):
                if n == sec['inds_header_in'][jj]:
                    sec['m_dots_in'][n] = m_dot * n_loop * (n_loop // 2 - jj) / (n_loop // 2)
                    jj += 1
                else:
                    sec['m_dots_in'][n] = sec['m_dots_in'][n - 1]

            # Step Inlet Header Temperatures Through Time (RK4)
            t_hold = sec['t_header_inlet'][:n_node_header].copy()
            L_cv_hold = sec['L_cv_inlet'][:n_cv_header].copy()
            vol_hold = sec['Vol_inlet'][:n_cv_header].copy()

            # Compute average temperature of control volume
            # TODO-NEEDS CONVERSION REVIEW: 0-based indexing — Fortran t_hold(2:n_node_header) → Python t_hold[1:n_node_header]
            sec['t_bar_inlet'][:n_cv_header] = (
                t_hold[1:n_node_header] + t_hold[:n_node_header - 1]
            ) / 2.0

            # K1
            sec['k1_ih'][:n_node_header] = dT_dt_inlet(sec['m_dots_in'][:n_cv_header], t_hold[:n_node_header], vol_hold[:n_cv_header], L_cv_hold[:n_cv_header], self.mc_header_mult.v, sec['t_bar_inlet'][:n_cv_header], n_node_header, self.fluid_ID.v, self.inlet_header_heat_loss.v)
            sec['t_hat_inlet'][:n_node_header] = t_hold[:n_node_header] + sec['k1_ih'][:n_node_header] * DTheta / 2.0
            sec['t_bar_hat_inlet'][:n_cv_header] = (
                sec['t_hat_inlet'][1:n_node_header] + sec['t_hat_inlet'][:n_node_header - 1]
            ) / 2.0

            # K2
            sec['k2_ih'][:n_node_header] = dT_dt_inlet(sec['m_dots_in'][:n_cv_header], sec['t_hat_inlet'][:n_node_header], vol_hold[:n_cv_header], L_cv_hold[:n_cv_header], self.mc_header_mult.v, sec['t_bar_hat_inlet'][:n_cv_header], n_node_header, self.fluid_ID.v, self.inlet_header_heat_loss.v)
            sec['t_hat_inlet'][:n_node_header] = t_hold[:n_node_header] + sec['k2_ih'][:n_node_header] * DTheta / 2.0
            sec['t_bar_hat_inlet'][:n_cv_header] = (
                sec['t_hat_inlet'][1:n_node_header] + sec['t_hat_inlet'][:n_node_header - 1]
            ) / 2.0

            # K3
            sec['k3_ih'][:n_node_header] = dT_dt_inlet(sec['m_dots_in'][:n_cv_header], sec['t_hat_inlet'][:n_node_header], vol_hold[:n_cv_header], L_cv_hold[:n_cv_header], self.mc_header_mult.v, sec['t_bar_hat_inlet'][:n_cv_header], n_node_header, self.fluid_ID.v, self.inlet_header_heat_loss.v)
            sec['t_hat_inlet'][:n_node_header] = t_hold[:n_node_header] + sec['k3_ih'][:n_node_header] * DTheta
            sec['t_bar_hat_inlet'][:n_cv_header] = (
                sec['t_hat_inlet'][1:n_node_header] + sec['t_hat_inlet'][:n_node_header - 1]
            ) / 2.0

            # K4
            sec['k4_ih'][:n_node_header] = dT_dt_inlet(sec['m_dots_in'][:n_cv_header], sec['t_hat_inlet'][:n_node_header], vol_hold[:n_cv_header], L_cv_hold[:n_cv_header], self.mc_header_mult.v, sec['t_bar_hat_inlet'][:n_cv_header], n_node_header, self.fluid_ID.v, self.inlet_header_heat_loss.v)

            sec['t_header_inlet'][:n_node_header] = t_hold[:n_node_header] + (
                sec['k1_ih'][:n_node_header] / 6.0
                + sec['k2_ih'][:n_node_header] / 3.0
                + sec['k3_ih'][:n_node_header] / 3.0
                + sec['k4_ih'][:n_node_header] / 6.0
            ) * DTheta

            # Solar Field Loop Calculations
            defocusing = 0.0
            defocus_groups = np.zeros(8)
            theta = self.Theta.v
            phi = self.Phi.v

            # TODO-NEEDS CONVERSION REVIEW: 0-based indexing — Fortran loop from 1 to n_loop
            for loop_idx in range(n_loop):
                loop = loop_idx + 1  # 1-based loop number for group/parity logic

                # Store loop temperatures in temporary array
                t_sf_hold = sec['t_sf'][:, loop_idx].copy()
                # TODO-NEEDS CONVERSION REVIEW: 0-based indexing — Fortran t_sf_hold(2:) → Python t_sf_hold[1:]
                sec['t_bar'][:] = (t_sf_hold[1:] + t_sf_hold[:n_nodes_per_loop - 1]) / 2.0

                if loop <= 13:
                    group_index = 0      # TODO-NEEDS CONVERSION REVIEW: 0-based indexing — Fortran group_index=1 → Python 0
                elif loop <= 26:
                    group_index = 1
                elif loop <= 39:
                    group_index = 2
                elif loop <= 52:
                    group_index = 3
                elif loop <= 65:
                    group_index = 4
                elif loop <= 78:
                    group_index = 5
                elif loop <= 91:
                    group_index = 6
                elif loop <= 104:
                    group_index = 7

                # Solar irradiation normal to collector
                sec['dni_array'][:] = self.ANI.v
                eta_row = Row_shadow(phi, self.row_distance.v, self.W_ap.v)

                # Defocusing Scheme
                defocus_groups[group_index] = 10.0 * defocus_groups[group_index] + 1.0

                # TODO-NEEDS CONVERSION REVIEW: 0-based indexing — sec['defocus_mode'][loop_idx] replaces defocus_mode(loop, sf_label)
                if sec['defocus_mode'][loop_idx] == 1.0:  # Collector 4 is not defocused
                    if t_sf_hold[n_nodes_per_loop - 1] > T_41:
                        sec['defocus_mode'][loop_idx] = 2.0

                elif sec['defocus_mode'][loop_idx] == 2.0:  # Collector 4 is defocused to setpoint 1
                    if t_sf_hold[n_nodes_per_loop - 1] > T_42:
                        sec['defocus_mode'][loop_idx] = 3.0
                    elif t_sf_hold[n_nodes_per_loop - 1] < T_41:
                        sec['defocus_mode'][loop_idx] = 1.0
                        defocus_groups[group_index] = (
                            defocus_groups[group_index]
                            - defocus_groups[group_index] % 10.0
                            + 1.0
                        )
                    else:
                        defocus_groups[group_index] = (
                            defocus_groups[group_index]
                            - defocus_groups[group_index] % 10.0
                            + 2.0
                        )

                elif sec['defocus_mode'][loop_idx] == 3.0:  # Collector 4 is defocused to setpoint 2
                    if t_sf_hold[n_nodes_per_loop - 1] > T_4d:
                        sec['defocus_mode'][loop_idx] = 4.0
                    elif t_sf_hold[n_nodes_per_loop - 1] < T_42:
                        sec['defocus_mode'][loop_idx] = 2.0
                        defocus_groups[group_index] = (
                            defocus_groups[group_index]
                            - defocus_groups[group_index] % 10.0
                            + 2.0
                        )
                    else:
                        defocus_groups[group_index] = (
                            defocus_groups[group_index]
                            - defocus_groups[group_index] % 10.0
                            + 3.0
                        )

                elif sec['defocus_mode'][loop_idx] == 4.0:  # Collector 4 is fully defocused
                    if t_sf_hold[n_nodes_per_loop - 1] < T_4d:
                        sec['defocus_mode'][loop_idx] = 3.0
                        defocus_groups[group_index] = (
                            defocus_groups[group_index]
                            - defocus_groups[group_index] % 10.0
                            + 3.0
                        )
                        sec['time_df'][loop_idx] = 0.0
                    else:
                        sec['time_df'][loop_idx] += DTheta
                        defocus_groups[group_index] = (
                            defocus_groups[group_index]
                            - defocus_groups[group_index] % 10.0
                            + 4.0
                        )

                    if sec['time_df'][loop_idx] > self.Time_lim_df.v:
                        sec['defocus_mode'][loop_idx] = 5.0
                        defocus_groups[group_index] = (
                            defocus_groups[group_index]
                            - defocus_groups[group_index] % 10.0
                            + 5.0
                        )

                else:  # All collectors are defocused
                    if t_sf_hold[n_nodes_per_loop - 1] < T_4d:
                        sec['defocus_mode'][loop_idx] = 3.0
                        defocus_groups[group_index] = (
                            defocus_groups[group_index]
                            - defocus_groups[group_index] % 10.0
                            + 3.0
                        )
                        sec['time_df'][loop_idx] = 0.0

                # End Loss Efficiency
                L_SCA = self.L_tot.v / n_SCA
                EndGain_pred = max(self.Ave_focal_length.v * math.tan(theta) - self.Distance_SCA.v, 0.0) / L_SCA
                EndLoss_pred = 1.0 - self.Ave_focal_length.v * math.tan(theta) / L_SCA
                if loop % 2 == 0:  # Loops going north
                    eta_endLoss = np.array([
                        EndLoss_pred,
                        EndLoss_pred + EndGain_pred,
                        EndLoss_pred + EndGain_pred,
                        EndLoss_pred,
                    ])
                else:  # Loops going south
                    eta_endLoss = np.array([
                        EndLoss_pred + EndGain_pred,
                        EndLoss_pred,
                        EndLoss_pred,
                        EndLoss_pred + EndGain_pred,
                    ])

                n_cv = n_nodes_per_loop - 1

                # If all collectors are defocused
                if sec['defocus_mode'][loop_idx] == 5.0:
                    sec['dni_array'][:] = 0.0
                else:
                    # TODO-NEEDS CONVERSION REVIEW: 0-based indexing — Fortran cc from 1 to n_SCA
                    for cc in range(1, n_SCA + 1):
                        eta_IAM = (
                            self.IAM_a0.v
                            + self.IAM_a1.v * theta / math.cos(theta)
                            + self.IAM_a2.v * theta ** 2 / math.cos(theta)
                        )
                        # TODO-NEEDS CONVERSION REVIEW: 0-based indexing — Fortran slice 1+(cc-1)*n_cv//n_sca:cc*n_cv//n_sca
                        # → Python (cc-1)*n_cv//n_SCA : cc*n_cv//n_SCA
                        sl = slice((cc - 1) * n_cv // n_SCA, cc * n_cv // n_SCA)
                        if cc == 4 or cc == 3:
                            sec['dni_array'][sl] = (
                                self.ANI.v * eta_IAM * eta_tot * eta_row
                                * sec['r_number'][cc - 1, loop_idx]  # TODO-NEEDS CONVERSION REVIEW: 0-based indexing
                                * eta_defocus[int(sec['defocus_mode'][loop_idx]) - 1]  # TODO-NEEDS CONVERSION REVIEW: 0-based indexing into eta_defocus
                                * eta_endLoss[cc - 1]  # TODO-NEEDS CONVERSION REVIEW: 0-based indexing
                            )
                        else:
                            sec['dni_array'][sl] = (
                                self.ANI.v * eta_IAM * eta_tot * eta_row
                                * sec['r_number'][cc - 1, loop_idx]  # TODO-NEEDS CONVERSION REVIEW: 0-based indexing
                                * eta_endLoss[cc - 1]  # TODO-NEEDS CONVERSION REVIEW: 0-based indexing
                            )

                # Update features matrix for NN
                # TODO-NEEDS CONVERSION REVIEW: 0-based indexing — Fortran state from 1 to 4 → Python 0 to 3
                for state in range(4):
                    n = sec['nCV_state'][state, loop_idx]
                    if state == 0:
                        curr_inds = sec['inds_pristine'][:n, loop_idx].copy()
                    elif state == 1:
                        curr_inds = sec['inds_lVacuum'][:n, loop_idx].copy()
                    elif state == 2:
                        curr_inds = sec['inds_bGlass'][:n, loop_idx].copy()
                    else:
                        curr_inds = sec['inds_H2'][:n, loop_idx].copy()

                    if n > 0:
                        minMax2 = sec['minMax2']
                        # HTF Temperature
                        # TODO-NEEDS CONVERSION REVIEW: 0-based indexing — Fortran minMax2(row,col) is 1-based; Python [row-1, col-1]
                        sec['features'][curr_inds, 0] = (
                            sec['t_bar'][curr_inds] - minMax2[2, state]
                        ) / (minMax2[3, state] - minMax2[2, state]) + 0.1
                        # Mass flow rate
                        sec['features'][curr_inds, 1] = (
                            m_dot * sec['m_dot_var'][loop_idx] - minMax2[0, state]
                        ) / (minMax2[1, state] - minMax2[0, state]) + 0.1
                        # Ambient Temperature
                        sec['features'][curr_inds, 2] = (
                            self.T_amb.v - minMax2[6, state]
                        ) / (minMax2[7, state] - minMax2[6, state]) + 0.1
                        # Wind Velocity
                        sec['features'][curr_inds, 3] = (
                            self.Wind.v - minMax2[10, state]
                        ) / (minMax2[11, state] - minMax2[10, state]) + 0.1
                        # Incident Solar Energy
                        sec['features'][curr_inds, 4] = (
                            sec['dni_array'][curr_inds] * self.W_ap.v - minMax2[12, state]
                        ) / (minMax2[13, state] - minMax2[12, state]) + 0.1
                        # Annulus Pressure
                        if state == 3:  # If hydrogen in annulus
                            for mm in range(n):
                                sec['features'][curr_inds[mm], 5] = (
                                    sec['H2_pressure'][mm, loop_idx] - minMax2[4, state]
                                ) / (minMax2[5, state] - minMax2[4, state]) + 0.1
                        else:
                            # If not hydrogen state, pressure in annulus is assumed to be unchanging
                            sec['features'][curr_inds, 5] = 0.1
                        # Sky Temperature
                        sec['features'][curr_inds, 6] = (
                            self.T_sky.v - minMax2[8, state]
                        ) / (minMax2[9, state] - minMax2[8, state]) + 0.1

                # Step through time with RK-4
                # Assign temporary arrays for function call
                inds_pristine_hold = sec['inds_pristine'][:, loop_idx].copy()
                inds_lVacuum_hold = sec['inds_lVacuum'][:, loop_idx].copy()
                inds_bGlass_hold = sec['inds_bGlass'][:, loop_idx].copy()
                inds_H2_hold = sec['inds_H2'][:, loop_idx].copy()
                nCV_state_hold = sec['nCV_state'][:, loop_idx].copy()

                # K1
                sec['k1'][:] = dt_dtime_NN(
                    t_sf_hold, sec['t_bar'], sec['features'],
                    m_dot * sec['m_dot_var'][loop_idx],
                    self.mc_receiver_mult.v, sec['L_segment'], sec['Vol'],
                    n_nodes_per_loop, nCV_state_hold,
                    inds_pristine_hold, inds_lVacuum_hold,
                    inds_bGlass_hold, inds_H2_hold,
                    self.fluid_ID.v,
                )
                t_hat = t_sf_hold + sec['k1'] * DTheta / 2.0
                t_bar_hat = (t_hat[1:] + t_hat[:n_nodes_per_loop - 1]) / 2.0
                for state in range(4):
                    n = nCV_state_hold[state]
                    if state == 0:
                        curr_inds = inds_pristine_hold[:n].copy()
                    elif state == 1:
                        curr_inds = inds_lVacuum_hold[:n].copy()
                    elif state == 2:
                        curr_inds = inds_bGlass_hold[:n].copy()
                    else:
                        curr_inds = inds_H2_hold[:n].copy()
                    if n > 0:
                        minMax2 = sec['minMax2']
                        # HTF Temperature
                        # TODO-NEEDS CONVERSION REVIEW: 0-based indexing
                        sec['features'][curr_inds, 0] = (
                            t_bar_hat[curr_inds] - minMax2[2, state]
                        ) / (minMax2[3, state] - minMax2[2, state]) + 0.1

                # K2
                sec['k2'][:] = dt_dtime_NN(
                    t_hat, t_bar_hat, sec['features'],
                    m_dot * sec['m_dot_var'][loop_idx],
                    self.mc_receiver_mult.v, sec['L_segment'], sec['Vol'],
                    n_nodes_per_loop, nCV_state_hold,
                    inds_pristine_hold, inds_lVacuum_hold,
                    inds_bGlass_hold, inds_H2_hold,
                    self.fluid_ID.v,
                )
                t_hat = t_sf_hold + sec['k2'] * DTheta / 2.0
                t_bar_hat = (t_hat[1:] + t_hat[:n_nodes_per_loop - 1]) / 2.0
                for state in range(4):
                    n = nCV_state_hold[state]
                    if state == 0:
                        curr_inds = inds_pristine_hold[:n].copy()
                    elif state == 1:
                        curr_inds = inds_lVacuum_hold[:n].copy()
                    elif state == 2:
                        curr_inds = inds_bGlass_hold[:n].copy()
                    else:
                        curr_inds = inds_H2_hold[:n].copy()
                    if n > 0:
                        minMax2 = sec['minMax2']
                        # HTF Temperature
                        # TODO-NEEDS CONVERSION REVIEW: 0-based indexing
                        sec['features'][curr_inds, 0] = (
                            t_bar_hat[curr_inds] - minMax2[2, state]
                        ) / (minMax2[3, state] - minMax2[2, state]) + 0.1

                # K3
                sec['k3'][:] = dt_dtime_NN(
                    t_hat, t_bar_hat, sec['features'],
                    m_dot * sec['m_dot_var'][loop_idx],
                    self.mc_receiver_mult.v, sec['L_segment'], sec['Vol'],
                    n_nodes_per_loop, nCV_state_hold,
                    inds_pristine_hold, inds_lVacuum_hold,
                    inds_bGlass_hold, inds_H2_hold,
                    self.fluid_ID.v,
                )
                t_hat = t_sf_hold + sec['k3'] * DTheta
                t_bar_hat = (t_hat[1:] + t_hat[:n_nodes_per_loop - 1]) / 2.0
                for state in range(4):
                    n = nCV_state_hold[state]
                    if state == 0:
                        curr_inds = inds_pristine_hold[:n].copy()
                    elif state == 1:
                        curr_inds = inds_lVacuum_hold[:n].copy()
                    elif state == 2:
                        curr_inds = inds_bGlass_hold[:n].copy()
                    else:
                        curr_inds = inds_H2_hold[:n].copy()
                    if n > 0:
                        minMax2 = sec['minMax2']
                        # HTF Temperature
                        # TODO-NEEDS CONVERSION REVIEW: 0-based indexing
                        sec['features'][curr_inds, 0] = (
                            t_bar_hat[curr_inds] - minMax2[2, state]
                        ) / (minMax2[3, state] - minMax2[2, state]) + 0.1

                # K4
                sec['k4'][:] = dt_dtime_NN(
                    t_hat, t_bar_hat, sec['features'],
                    m_dot * sec['m_dot_var'][loop_idx],
                    self.mc_receiver_mult.v, sec['L_segment'], sec['Vol'],
                    n_nodes_per_loop, nCV_state_hold,
                    inds_pristine_hold, inds_lVacuum_hold,
                    inds_bGlass_hold, inds_H2_hold,
                    self.fluid_ID.v,
                )

                # Step Through Time
                # TODO-NEEDS CONVERSION REVIEW: 0-based indexing — Fortran t_sf(2:, loop, sf) → Python t_sf[1:, loop_idx]
                sec['t_sf'][1:, loop_idx] = t_sf_hold[1:] + (
                    sec['k1'][1:] / 6.0
                    + sec['k2'][1:] / 3.0
                    + sec['k3'][1:] / 3.0
                    + sec['k4'][1:] / 6.0
                ) * DTheta

            # Persist defocus_groups so outputs can be set in subsequent iterations
            sec['defocus_groups'] = defocus_groups.copy()

            finish = time.perf_counter()
            return  # End of is_converged block

        # HTF Mass in SF Computation (Only do once because temperatures aren't changing within timestep)
        if self.model.iteration == 0:
            sec['mass_HTF_hold'] = 0.0

            # Inlet Header
            # TODO-NEEDS CONVERSION REVIEW: 0-based indexing — Fortran n from 1 to n_cv_header
            for n in range(n_cv_header):
                sec['mass_HTF_hold'] += sec['Vol_inlet'][n] * Inc.density(
                    fluid=self.fluid_ID.v,
                    T=(sec['t_header_inlet'][n] + sec['t_header_inlet'][n + 1]) / 2.0,
                    P=0.0,
                )

            # Solar Field
            sec['t_bar_sf'][:] = 0.0
            sec['t4Ave_hold'] = 0.0
            # TODO-NEEDS CONVERSION REVIEW: 0-based indexing — Fortran CEILING((real(n_nodes_per_loop)-1.0)/real(n_SCA)/2.0)
            # → Python math.ceil; ind_4 converted from 1-based to 0-based
            ind_4 = (n_nodes_per_loop - 1) - math.ceil(
                (float(n_nodes_per_loop) - 1.0) / float(n_SCA) / 2.0
            )  # TODO-NEEDS CONVERSION REVIEW: 0-based indexing

            sec['temperature_groups'][:] = 0.0

            # TODO-NEEDS CONVERSION REVIEW: 0-based indexing — Fortran loop from 1 to n_loop
            for loop_idx in range(n_loop):
                loop = loop_idx + 1  # 1-based for group logic

                t_sf_hold = sec['t_sf'][:, loop_idx].copy()
                sec['t_bar'][:] = (t_sf_hold[1:] + t_sf_hold[:n_nodes_per_loop - 1]) / 2.0

                if loop <= 13:
                    group_index = 0  # TODO-NEEDS CONVERSION REVIEW: 0-based indexing
                elif loop <= 26:
                    group_index = 1
                elif loop <= 39:
                    group_index = 2
                elif loop <= 52:
                    group_index = 3
                elif loop <= 65:
                    group_index = 4
                elif loop <= 78:
                    group_index = 5
                elif loop <= 91:
                    group_index = 6
                elif loop <= 104:
                    group_index = 7

                # Temperature levels
                t_ind4 = sec['t_sf'][ind_4, loop_idx]
                if t_ind4 <= 370.0:                        # ~200 degrees F
                    sec['temperature_groups'][group_index] = 10.0 * sec['temperature_groups'][group_index] + 1.0
                elif t_ind4 <= 427.0:                      # ~308 degrees F
                    sec['temperature_groups'][group_index] = 10.0 * sec['temperature_groups'][group_index] + 2.0
                elif t_ind4 <= 450.0:                      # ~350 degrees F
                    sec['temperature_groups'][group_index] = 10.0 * sec['temperature_groups'][group_index] + 3.0
                elif t_ind4 <= 483.0:                      # ~400 degrees F
                    sec['temperature_groups'][group_index] = 10.0 * sec['temperature_groups'][group_index] + 4.0
                elif t_ind4 <= 538.0:                      # ~508 degrees F
                    sec['temperature_groups'][group_index] = 10.0 * sec['temperature_groups'][group_index] + 5.0
                elif t_ind4 <= 566.0:                      # ~560 degrees F
                    sec['temperature_groups'][group_index] = 10.0 * sec['temperature_groups'][group_index] + 6.0
                elif t_ind4 <= 594.0:                      # ~609 degrees F
                    sec['temperature_groups'][group_index] = 10.0 * sec['temperature_groups'][group_index] + 7.0
                elif t_ind4 <= 616.0:                      # ~650 degrees F
                    sec['temperature_groups'][group_index] = 10.0 * sec['temperature_groups'][group_index] + 8.0
                elif t_ind4 > 616.0:                       # > 650 degrees F
                    sec['temperature_groups'][group_index] = 10.0 * sec['temperature_groups'][group_index] + 9.0

                # TODO-NEEDS CONVERSION REVIEW: 0-based indexing — Fortran n from 1 to n_nodes_per_loop-1
                for n in range(n_nodes_per_loop - 1):
                    # HTF mass in control volume
                    sec['mass_HTF_hold'] += sec['Vol'] * Inc.density(
                        fluid=self.fluid_ID.v, T=sec['t_bar'][n], P=0.0
                    )
                    # Average loop temperature for pressure drop calculations
                    sec['t_bar_sf'][n] += sec['t_bar'][n]
                    # Average temperature at sensor 4
                    if n == ind_4:
                        sec['t4Ave_hold'] += t_sf_hold[n]

            # Average across all loops
            sec['t4Ave_hold'] /= n_loop
            sec['t_bar_sf'][:] /= n_loop

            # Return Header
            # TODO-NEEDS CONVERSION REVIEW: 0-based indexing — Fortran n from 1 to n_cv_header+1
            for n in range(n_cv_header + 1):
                sec['mass_HTF_hold'] += sec['Vol_return'][n] * Inc.density(
                    fluid=self.fluid_ID.v,
                    T=(sec['t_header_return'][n] + sec['t_header_return'][n + 1]) / 2.0,
                    P=0.0,
                )

        # Perform all hydraulic calculations here (allow new solution for every iteration)
        m_dot_htf = self.MassFlow.v / n_loop

        # Pressure Drop accounting for Inlet, Outlet, and Cross-Over-Piping of solar field loop
        DP_IOCOP = PressureDrop(self.fluid_ID.v, m_dot_htf, (sec['t_bar_sf'][0] + sec['t_bar_sf'][n_nodes_per_loop - 2]) / 2.0, 1.0, self.D_receiver.v, self.Roughness_pipe.v, (40.0 + self.row_distance.v), 0.0, 0.0, 2.0, 0.0, 0.0, 2.0, 0.0, 0.0, 2.0, 1.0, 0.0)

        # Pressure Drop Across Solar Field Loop
        DP_loop = 0.0
        SCA_ind_orig = (n_nodes_per_loop - 1) / n_SCA
        SCA_ind = SCA_ind_orig
        # TODO-NEEDS CONVERSION REVIEW: 0-based indexing — Fortran n from 1 to n_nodes_per_loop-1
        for n in range(n_nodes_per_loop - 1):
            # Account for extra fittings on the first HCE
            if n == 0:
                x1 = 10.0
                x2 = 3.0
            # Account for the two ball joints that connects SCAs
            elif n >= SCA_ind:
                x1 = 0.0
                x2 = 2.0
                SCA_ind += SCA_ind_orig
            # Remaining loop control volumes do not have ball joint connections
            else:
                x1 = 0.0
                x2 = 0.0
            DP_loop += PressureDrop(self.fluid_ID.v, m_dot_htf, sec['t_bar_sf'][n], 1.0, self.D_receiver.v, self.Roughness_pipe.v, sec['L_segment'], 0.0, 0.0, x1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, x2)

        # Pressure Drop Across Inlet Header
        DP_inlet_header = 0.0
        cc = 1
        loop_count = 1
        m_dot_header = self.MassFlow.v
        # TODO-NEEDS CONVERSION REVIEW: 0-based indexing — Fortran n from 1 to n_cv_header
        x1 = 0.0
        x2 = 0.0
        for n in range(n_cv_header):
            # Check if there is a pipe contraction
            if n > 0:
                if sec['D_inlet'][n] < sec['D_inlet'][n - 1]:
                    x1 = 1.0
                else:
                    x1 = 0.0

            # Check if control volume is part of an expansion loop with two long elbows
            if cc > 1:
                # If it is the first of two CV's modeling the expansion loop, change the mass flow
                if cc == 2:
                    m_dot_header -= 2.0 * m_dot_htf
                cc += 1
                x2 = 2.0
                if cc == 4:
                    cc = 1
                    loop_count += 1
            # Else control volume is the segment of header between sf loop inlets
            else:
                x2 = 0.0
                cc += 1
                loop_count += 1
                m_dot_header -= 2.0 * m_dot_htf

            DP_inlet_header += PressureDrop(self.fluid_ID.v, m_dot_header, (sec['t_header_inlet'][n] + sec['t_header_inlet'][n + 1]) / 2.0, 1.0, sec['D_inlet'][n], self.Roughness_pipe.v, sec['L_cv_inlet'][n], 0.0, x1, 0.0, 0.0, x2, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)

        # Pressure Drop Across Return Header
        DP_return_header = 0.0
        cc = 1
        loop_count = 1
        m_dot_header = 0.0
        x1 = 0.0
        x2 = 0.0
        # TODO-NEEDS CONVERSION REVIEW: 0-based indexing — Fortran n from 1 to n_cv_header
        for n in range(n_cv_header):
            # Check if there is a pipe contraction
            if n > 0:
                if sec['D_inlet'][n] > sec['D_inlet'][n - 1]:
                    x1 = 1.0
                else:
                    x1 = 0.0

            # Check if control volume is part of an expansion loop with two long elbows
            if cc > 1:
                # If it is the first of two CV's modeling the expansion loop, change the mass flow
                if cc == 2:
                    m_dot_header += 2.0 * m_dot_htf
                cc += 1
                x2 = 2.0
                if cc == 4:
                    cc = 1
                    loop_count += 1
            # Else control volume is the segment of header between sf loop inlets
            else:
                x2 = 0.0
                cc += 1
                loop_count += 1
                m_dot_header += 2.0 * m_dot_htf

            DP_return_header += PressureDrop(self.fluid_ID.v, m_dot_header, (sec['t_header_inlet'][n] + sec['t_header_inlet'][n + 1]) / 2.0, 1.0, sec['D_return'][n], self.Roughness_pipe.v, sec['L_cv_return'][n], x1, 0.0, 0.0, 0.0, x2, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)

        # Pressure Drop Across Entire Solar Field
        DP_tot = DP_return_header + DP_inlet_header + DP_loop + DP_IOCOP

        defocusing = 0.0
        # TODO-NEEDS CONVERSION REVIEW: 0-based indexing — Fortran loop from 1 to n_loop
        for loop_idx in range(n_loop):
            defocusing += sec['defocus_mode'][loop_idx]

        # Set the Outputs from this Model
        self.MassFlow_out.v = self.MassFlow.v                                     # MassFlow
        self.Pressure_out.v = self.Pressure.v - DP_tot                            # Pressure
        # TODO-NEEDS CONVERSION REVIEW: 0-based indexing — Fortran t_header_return(n_cv_header+1, sf) → Python [n_cv_header]
        self.Temperature_out.v = sec['t_header_return'][n_cv_header]              # Temperature
        self.T4Ave.v = sec['t4Ave_hold']                                          # T4Ave
        self.Mass_Counter_out.v = sec['mass_HTF_hold']                            # HTF Mass Counter
        self.defocusing_out.v = defocusing - n_loop                               # defocusing information

        self.defocus_groups_1.v = sec['defocus_groups'][0]                         # defocusing status group 1
        self.defocus_groups_2.v = sec['defocus_groups'][1]                         # defocusing status group 2
        self.defocus_groups_3.v = sec['defocus_groups'][2]                         # defocusing status group 3
        self.defocus_groups_4.v = sec['defocus_groups'][3]                         # defocusing status group 4
        self.defocus_groups_5.v = sec['defocus_groups'][4]                         # defocusing status group 5
        self.defocus_groups_6.v = sec['defocus_groups'][5]                         # defocusing status group 6
        self.defocus_groups_7.v = sec['defocus_groups'][6]                         # defocusing status group 7
        self.defocus_groups_8.v = sec['defocus_groups'][7]                         # defocusing status group 8
        self.temperature_groups_1.v = sec['temperature_groups'][0]
        self.temperature_groups_2.v = sec['temperature_groups'][1]
        self.temperature_groups_3.v = sec['temperature_groups'][2]
        self.temperature_groups_4.v = sec['temperature_groups'][3]
        self.temperature_groups_5.v = sec['temperature_groups'][4]
        self.temperature_groups_6.v = sec['temperature_groups'][5]
        self.temperature_groups_7.v = sec['temperature_groups'][6]
        self.temperature_groups_8.v = sec['temperature_groups'][7]
