"""Cooling Tower component model (Type 162)."""

# This component models the performance of a cooling tower as outlined "Effectiveness Models for Cooling Tower and
#  Cooling Coils", submitted to the ASME Journal or Heat and Mass Transfer, 1987.
#
# Written By: Solar Energy Laboratory - University of Wisconsin, Madison and Thermal Energy System Specialists, LLC
#
# Revision History:
#  05/1993 - JWT: removal of English units!
#  03/2004 - TPM: for TRNSYS 16
#  01/2008 - DEB: BUGFIX: OUT(19) was getting reset to TI every iteration if there was more than one instance
#                  of Type51 in a simulation.
#  04/2009 - MJD: set initial output array to appropriate input(s) or zero(s)
#  07/2009 - TPM: conversion to version 17 coding standards
#  03/2010 - DEB: BUGFIX: added a check to make sure that the water outlet temperature is less than the water
#                  inlet temperature in performance map.
#  07/2012 - MJD: corrected the final call to MoistAirProperties to be mode 6 instead of mode 1.
#  11/2012 - DEB: BUGFIXES: 1. handle the situation in which there is water flow but all the fan speeds are set to
#                  zero. 2. remove the line that reset the sump temperature to its initial value every time a new
#                  instance of Type51 is encountered. 3. REVISIT: undid MJD fix from 07/2012 as it generates a warning
#                  at every timestep.
#                 MODIFICATIONS: 1. added a control mode in which the tower computes the flow required to meet a user
#                  specified water outlet temperature. 2. added outputs for the fan relative speed. 3. added inputs
#                  for sump UA and sump auxiliary energy. Modified the sump temperature diffeq to account for the new
#                  energy terms.
#  08/2013 - DEB: ignore sump losses and sump energy gain input if the sump volume is set to -1.
#  02/2014 - DEB: changed the air mixing call to be mode 7 instead of mode 6 (reset w to wsat at h instead of resetting
#                  h to hsat at w). Also removed the psych warning from this case because the user doesn't really need
#                  to worry if condensation occurs. Removed the same error check from a number of other psych calls.
#  09/2014 - JWT: Added a check on the calculation of the effectiveness
#  11/2014 - DEB: modification of the 11/2012 control mode so that the cell fan control signal inputs (which were
#                  ignored in this mode) can now be used to tell the model whether each cell is available at the current
#                  time. The tower now determines the flow required to meet the water outlet temperature and then splits
#                  the flow equally between the available cells instead of between all cells.
#  02/2015 - DEB: 1. implemented the VARIANT keyword. VARIANT 0 corresponds to the PIO list and algorithms before the
#                  11/2012 modifications and VARIANT 1 corresponds to the PIO list and algorithms post 11/2012 mods.
#                  2. added SSR code.
#  03/2016 - DEB: bugfix in the fan power calculation. Fan power was being accumulated as it looped through cells but
#                  wasn't being reset to zero until the fan turned off.
#  05/2016 - DEB: 1. renumbered Type from 51 to 162 because of the changes in the parameter and input list. Had to
#                  remove the VARIANT coding as it is not implemented in the Studio. 2. Changed the units of parameter 6
#                  and output 3 from [kW] to [kJ/h] for consistency with other HVAC Types. 3. Fixed some of the format
#                  statements.
#  02/2017 - DEB: another bugfix in the fan power calculation. Fan power was being accumulated during iterations to
#                  find the necessary speed to match a desired temperature but wasn't being reset before recomputing at
#                  a different speed iteration.
#  04/2017 - DEB:  unit codes were incorrect for inputs 6 and 7.
#  09/2019 - JWT: Fixed a case with high air flow rates and low water flow rates that was causing a floating point issue.
#                  It turns out that the equation simplifies for values of mstar > 5 with nearly identical solutions.
#
# REVISIT: in testing this type, check and make sure that once it finds the total tower flow rate that it wants and
# tries to split that flow between available cells, that it does something logical if there aren't enough cells available.
# ----------------------------------------------------------------------------------------------------------------------
# Copyright 2019 Solar Energy Laboratory, University of Wisconsin-Madison and Thermal Energy System Specialists, LLC.
#    All rights reserved.

import numpy as np

from eeslib.fluid_properties import humid_air  
from esclab.simulate import Component

ndmax = 50
nMaxCells = 10

# Data constants
cpw = 4.186
rhow = 1000.
trefw = 0.
patm = 1.0
tol = 1.e-03
imax = 10
iu = 1
# emode determines how errors are handled from calls to the psych subroutine. If emode is: 0 - no errors
#  messages will be printed, 1 - error messages will be printed only once per simulation, 2 - error messages
#  will be printed every timestep that they occur.
emode = 2


def ts(hs):
    """
    Correlation for saturation temperature in terms of saturation enthalpy.
    The reference state for dry air enthalpies is: h=0.0 at 0 deg C
    The reference state for liquid water enthalpies is: h=0.0 at 0 deg C
    Correlation is in SI units and is good from 9.473 kJ/kg to 355.137 kJ/kg (0 to 55 C).
    (Correlation is for a total system pressure of 1 atmosphere.)
    """
    return (-5.79013 + 6.64030e-01 * hs - 5.07802e-03 * hs**2 + 2.80381e-05 * hs**3
            - 9.47051e-08 * hs**4 + 1.72758e-10 * hs**5 - 1.29547e-13 * hs**6)


class CoolingTower(Component):
    """
    Cooling Tower model using the NTU-effectiveness method (TRNSYS Type 162).

    Models the performance of a cooling tower as outlined in "Effectiveness Models for Cooling Tower and
    Cooling Coils", ASME Journal of Heat and Mass Transfer, 1987.
    """
    trnsys_type = "162"

    # Parameters
    # 1 | mode          | 0 = performance map, 1 = user-supplied coefficients, 2 = curve-fit from data
    mode = Component.Parameter()
    # 2 | controlMode   | 0 = relative fan speed provided as inputs, 1 = fan speed internally computed to temperature seek
    controlMode = Component.Parameter()
    # 3 | igeom         | 1 = counterflow, 2 = crossflow
    igeom = Component.Parameter()
    # 4 | maxcel        | number of possible cells to operate
    maxcel = Component.Parameter()
    # 5 | fanmax        | maximum tower air flow rate (m**3 dry air/hr)
    fanmax = Component.Parameter()
    # 6 | pwrmax        | maximum individual tower fan power (kJ/h)
    pwrmax = Component.Parameter()
    # 7 | fanoff        | natural convective air flow rate (m**3 dry air/hr)
    fanoff = Component.Parameter()
    # 8 | vsump         | volume of sump (m**3)
    vsump = Component.Parameter()
    # 9 | uaSump        | overall thermal loss coefficient from the sump [kJ/h.K]
    uaSump = Component.Parameter()
    # 10 | tstart       | initial sump temperature (C)
    tstart = Component.Parameter()
    # 11 | c1           | mass transfer constant (mode 1) or logical unit for data file (mode 2)
    c1 = Component.Parameter()
    # 12 | c2_raw       | mass transfer exponent minus 1 (mode 1) or number of data points (mode 2)
    c2_raw = Component.Parameter()
    # 13 | iptr         | print flag for data echo (mode 2 only)
    iptr = Component.Parameter()

    # Inputs
    # 1 | twi    | inlet water temperature (C)
    twi = Component.Input()
    # 2 | flwi   | inlet water flow rate (kg/hr)
    flwi = Component.Input()
    # 3 | tdb    | ambient dry bulb temperature (C)
    tdb = Component.Input()
    # 4 | twb    | ambient wet bulb temperature (C)
    twb = Component.Input()
    # 5 | tWant  | desired outlet water temperature for controlMode 1 (C)
    tWant = Component.Input()
    # 6 | tMain  | temperature of make-up water for sump (C)
    tMain = Component.Input()
    # 7 | qSump  | auxiliary energy added to the sump (kJ/h)
    qSump = Component.Input()
    # 8..7+maxcel | frac_i | relative fan speed for each cell (or cell availability flag in controlMode 1)
    # TODO-NEEDS CONVERSION REVIEW: variable number of inputs based on maxcel parameter; stored in list self._frac_inputs

    # Outputs
    # 1  | tsump   | sump exit water temperature (C)
    tsump = Component.Output()
    # 2  | flwi_out| water flow rate out (same as input) (kg/hr)
    flwi_out = Component.Output()
    # 3  | power   | total power consumption of tower fans (kJ/h)
    power = Component.Output()
    # 4  | qcells  | total heat rejection of tower cells (kJ/h)
    qcells = Component.Output()
    # 5  | two     | overall exit water temperature of cells (C)
    two = Component.Output()
    # 6  | flwl    | water loss (evaporation) rate (kg/hr)
    flwl = Component.Output()
    # 7  | taomix  | "mixing cup" dry bulb temperature of total exiting air (C)
    taomix = Component.Output()
    # 8  | twbmix  | "mixing cup" wet bulb temperature of total exiting air (C)
    twbmix = Component.Output()
    # 9  | waotot  | "mixing cup" humidity ratio of total exiting air
    waotot = Component.Output()
    # 10 | flatot  | total air flow rate (kg/hr)
    flatot = Component.Output()
    # 11 | deltau  | internal energy change of sump (kJ)
    deltau = Component.Output()
    # 12..11+maxcel | frac_out_i | relative fan speed of each cell
    # TODO-NEEDS CONVERSION REVIEW: variable number of outputs based on maxcel parameter; stored in list self._frac_outputs

    def presim_setup(self, **kwargs):
        super().presim_setup(**kwargs)

        mode_val = round(self.mode.v)
        igeom_val = round(self.igeom.v)
        maxcel_val = round(self.maxcel.v)
        ti = self.tstart.v

        if (mode_val != 1 and mode_val != 2):
            raise ValueError('The mode must be 1 or 2.')
        if (round(self.controlMode.v) != 0 and round(self.controlMode.v) != 1):
            raise ValueError('The control mode parameter must be 0 or 1.')
        if (igeom_val != 1 and igeom_val != 2):
            raise ValueError('The flow geometry parameter must be 1 or 2.')
        if (maxcel_val <= 0 or maxcel_val > 8):
            raise ValueError(f'The number of cells in the cooling tower must be between 1 and {nMaxCells}.')

        msump = self.vsump.v * rhow

        # Mode 2: determine ntu parameters from test data, curve-fit with linear regression
        if mode_val == 2:
            lu = round(self.c1.v)
            ndata = round(self.c2_raw.v)
            iptr_val = round(self.iptr.v)

            if ndata > ndmax:
                raise RuntimeError('Exceeded the maximum number of data points allowed')

            # TODO-NEEDS CONVERSION REVIEW: File I/O using Fortran logical unit; Python uses file paths directly
            # Open file and read performance data
            ok = False
            mstard = np.zeros(ndmax)
            epsd = np.zeros(ndmax)
            x = np.zeros((2, ndmax))
            y_arr = np.zeros(ndmax)
            coef = np.zeros(2)
            iwarn = 0

            try:
                with open(str(lu), 'r') as f:
                    for n in range(ndata):
                        line = f.readline()
                        vals = line.split()
                        fana, tdb_d, twb_d, flw1, tw1, tw2 = [float(v) for v in vals[:6]]

                        if tw1 < tw2:
                            raise ValueError('The cooling tower data contains a point at which the water outlet temperature greater than the water inlet temperature. The point needs to be removed in order to generate a curve fit for the data.')

                        psydat = humid_air(['T_wb','W','H','V'],T_db=tdb_d, T_wb=twb_d, P=patm)
                        twb_d = psydat['T_wb']
                        wa1 = psydat['W']
                        ha1 = psydat['H']
                        rhowa = 1./psydat['V']

                        if tw2 <= twb_d:
                            raise ValueError('The cooling tower data has a water outlet temperature less than or equal to the inlet air wet bulb temperature')

                        psydat1 = humid_air(['W','H'], T_db=tw1, T_wb=tw1, P=patm)
                        ww1 = psydat1['W']
                        hw1 = psydat1['H']
                        psydat2 = humid_air(['W','H'], T_db=tw2, T_wb=tw2, P=patm)
                        ww2 = psydat2['W']
                        hw2 = psydat2['H']

                        fla = rhowa * fana
                        ra = flw1 / fla
                        cs = (hw1 - hw2) / (tw1 - tw2)
                        mstar = cs / ra / cpw
                        qmax = fla * (hw1 - ha1)
                        eps = flw1 * cpw * (tw1 - tw2) / qmax
                        if eps >= 1.0:
                            raise ValueError('The cooling tower data gives an air-effectiveness that is greater than or equal to one')
                        eps = min(eps, (1. / mstar - 1.0e-6))

                        # Iteration is necessary to determine the air-side heat transfer effectiveness
                        for iter_n in range(imax):
                            epsold = eps
                            if igeom_val == 1:
                                c = (1. - eps) / (1 - mstar * eps)
                                ntu = -np.log(c) / (1. - mstar)
                            else:
                                ntu = -np.log(1 + np.log(1 - mstar * eps) / mstar)
                            ha2 = ha1 + eps * (hw1 - ha1)
                            hw = ha1 + (ha2 - ha1) / (1. - np.exp(-ntu))
                            tw = ts(hw)
                            hcheck = hw
                            if hcheck < 9.473 or hcheck > 355.137:
                                iwarn += 1
                                # note: inclusion of this warning counter can sometimes cause simulation errors
                                # (due to too many warnings). Not clear what the preferred fix would be.
                            psydat_tw = humid_air(['W','H'], T_db=tw, T_wb=tw, P=patm)
                            ww = psydat_tw['W']
                            hw = psydat_tw['H']
                            wa2 = ww + (wa1 - ww) * np.exp(-ntu)
                            flw2 = flw1 - fla * (wa2 - wa1)
                            qcell = flw1 * cpw * (tw1 - trefw) - flw2 * cpw * (tw2 - trefw)
                            eps = qcell / qmax
                            eps = min(eps, (1. / mstar - 1.0e-6), 0.9999999)
                            if abs(eps - epsold) <= tol:
                                break

                        # Store data for curve-fitting and later use
                        x[0, n] = 1.
                        x[1, n] = np.log(ra)
                        y_arr[n] = np.log(ntu)
                        mstard[n] = mstar
                        epsd[n] = eps

                        # Check for variation in the flowrate data
                        if not ok:
                            for j in range(n):
                                diff = abs(np.exp(x[1, j]) - np.exp(x[1, n]))
                                if diff > (0.1 * np.exp(x[1, n])):
                                    ok = True

            except Exception:
                raise RuntimeError('Error reading performance data file for Type 162 Cooling Tower')

            if not ok:
                raise ValueError('The ratio of the mass flowrate of water to the mass flowrate of dry air calculated from the performance data does not have enough variation for a good curve fit. New cooling tower data is required.')

            # Do curve-fit and store resulting parameters in parameter array
            # TODO-NEEDS LIBRARY: LinearRegression - use numpy least squares instead
            ndata_val = ndata
            A = x[:, :ndata_val].T  # shape (ndata, 2)
            coef, _, _, _ = np.linalg.lstsq(A, y_arr[:ndata_val], rcond=None)

            self.mode.v = 1.0
            self.c1.v = np.exp(coef[0])
            self.c2_raw.v = coef[1] - 1.

        # Coefficients of ntu correlation for either mode 1 or 2
        self._c1 = self.c1.v
        self._c2 = 1.0 + self.c2_raw.v

        # Sump mass
        self._msump = self.vsump.v * rhow

        # Stored sump temperature (dynamic array slot 1 in TRNSYS)
        self._ti = ti

        # Set the Initial Values of the Outputs
        self.tsump.v = self.twi.v
        self.flwi_out.v = self.flwi.v
        self.power.v = 0.
        self.qcells.v = 0.
        self.two.v = self.twi.v
        self.flwl.v = 0.
        self.taomix.v = self.tdb.v
        self.twbmix.v = self.twb.v
        self.waotot.v = 0.
        self.flatot.v = 0.
        self.deltau.v = 0.
        # TODO-NEEDS CONVERSION REVIEW: variable-length cell outputs (11+i) initialized to frac inputs (6+i)
        # These are managed via self._frac_out list in calculate()

    def calculate(self):
        super().calculate()

        maxcel_val = round(self.maxcel.v)

        # Re-read parameters (needed if another unit of this type was called last)
        mode_val = round(self.mode.v)
        controlMode_val = round(self.controlMode.v)
        igeom_val = round(self.igeom.v)
        c1_val = self.c1.v
        c2_val = 1.0 + self.c2_raw.v
        msump = self.vsump.v * rhow

        # Get the Current Inputs to the Model
        # TODO-NEEDS CONVERSION REVIEW: inputs 8..7+maxcel are variable; accessed by index in Fortran
        # In Python, frac_inputs must be fetched separately; placeholder list shown below
        # self._frac_inputs should be a list of input values for cells 1..maxcel
        # For now, this conversion maps them symbolically. Integration with Component input indexing needed.

        # Local working state (shared with inner TowerCellLoop)
        frac = np.zeros(nMaxCells)
        two = self.twi.v
        flwl = 0.
        qcells = 0.
        flatot = 0.
        airsum = 0.
        tNow = self.twi.v
        wa1 = 0.
        ha1 = 0.
        rhowa = 1.0
        iwarn = 0

        def get_cell_input(i):
            """
            Get the relative fan speed / availability input for cell i (1-indexed).
            TODO-NEEDS CONVERSION REVIEW: getInputValue(7+i) accesses the i-th cell input.
            This must be connected to the Variable-length cell inputs in the esclab model.
            """
            # TODO-NEEDS CONVERSION REVIEW: variable-length input access; requires model integration
            return 0.0

        def tower_cell_loop():
            """
            Subroutine TowerCellLoop: loop through all tower cells and compute heat/mass transfer.
            Returns (two, flwl, qcells, flatot, airsum, tNow, wa1, ha1, rhowa).
            """
            nonlocal iwarn
            two_loc = self.twi.v
            flwl_loc = 0.
            qcells_loc = 0.
            flatot_loc = 0.
            airsum_loc = 0.
            tNow_loc = self.twi.v
            wa1_loc = 0.
            ha1_loc = 0.
            rhowa_loc = 1.0
            twb_loc = self.twb.v

            ncell_local = int(sum(1 for i in range(1, maxcel_val + 1) if abs(frac[i - 1]) > tol))

            if ncell_local > 0 and self.flwi.v > 0.:
                tw1 = self.twi.v
                flw1 = self.flwi.v / ncell_local
                # TODO-NEEDS LIBRARY: humid_air for moist air property calculations
                psydat_amb = humid_air(['T_wb','W','H','V'], T_db=self.tdb.v, T_wb=self.twb.v, P=patm)
                twb_loc = psydat_amb['T_wb']
                wa1_loc = psydat_amb['W']
                ha1_loc = psydat_amb['H']
                rhowa_loc = 1./psydat_amb['V']

                # psydat(2): drybulb temperature, psydat(3): wetbulb temperature
                psydat_tw1 = humid_air(['W','H'], T_db=tw1, T_wb=tw1, P=patm)
                ww1 = psydat_tw1['W']
                hw1 = psydat_tw1['H']
                tw2 = twb_loc

                # Individual cell analysis
                for i in range(1, ncell_local + 1):
                    fla_i = rhowa_loc * self.fanmax.v * frac[i - 1]
                    if fla_i < 1.0e-6:
                        # water loss, tower heat rejection, and quantities used in calculating the air outlet
                        # state are unchanged since there is no air flow in this cell.
                        pass
                    else:
                        ra1 = flw1 / fla_i
                        ntu = c1_val * ra1**c2_val
                        for iter_n in range(imax):
                            told = tw2
                            psydat_tw2 = humid_air(['W','H'], T_db=tw2, T_wb=tw2, P=patm)
                            ww2 = psydat_tw2['W']
                            hw2 = psydat_tw2['H']
                            cs = (hw1 - hw2) / (tw1 - tw2)
                            mstar = cs / ra1 / cpw

                            if igeom_val == 1:
                                if mstar > 5.:
                                    eps = 1. / mstar  # the equation for c in the else section degrades to 1/mstar for values of mstar>5
                                else:
                                    c_val = np.exp(-ntu * (1. - mstar))
                                    eps = (1. - c_val) / (1. - mstar * c_val)
                                # if eps != eps then it has gotten set to the NaN condition.
                                if eps != eps:
                                    eps = 0.00001
                            else:
                                eps = (1 - np.exp(-mstar * (1. - np.exp(-ntu)))) / mstar

                            ha2 = ha1_loc + eps * (hw1 - ha1_loc)
                            hw = ha1_loc + (ha2 - ha1_loc) / (1. - np.exp(-ntu))
                            tw = ts(hw)
                            hcheck = hw

                            # REVISIT: this check can cause problems. Is there a better way we can do this?
                            if iter_n < 3 and (hcheck < 9.473 or hcheck > 355.137):
                                iwarn += 1
                                # The correlation for the saturation temperature used in the cooling tower model was
                                # used with an enthalpy value outside the allowable range.

                            psydat_tw_i = humid_air(['W','H'], T_db=tw, T_wb=tw, P=patm)
                            ww = psydat_tw_i['W']
                            hw = psydat_tw_i['H']
                            wa2 = ww + (wa1_loc - ww) * np.exp(-ntu)
                            ra2 = ra1 - (wa2 - wa1_loc)
                            tw2 = trefw + (ra1 * cpw * (tw1 - trefw) - (ha2 - ha1_loc)) / ra2 / cpw
                            if abs(tw2 - told) <= tol:
                                break

                        # Determine cumulative water loss, tower heat rejection, and quantities used
                        # in calculating the air outlet state.
                        flwl_loc = flwl_loc + fla_i * (wa2 - wa1_loc)
                        qcells_loc = qcells_loc + fla_i * (ha2 - ha1_loc)
                        flatot_loc = flatot_loc + fla_i
                        airsum_loc = airsum_loc + fla_i * ha2

            # Determine quantities required for "mixing cup" exit air conditions if there is no water flow
            elif ncell_local > 0 and self.flwi.v < 1.0e-6:
                psydat_amb2 = humid_air(T_db=self.tdb.v, T_wb=self.twb.v, P=patm)
                twb_loc = psydat_amb2['T_wb']
                wa1_loc = psydat_amb2['W']
                ha1_loc = psydat_amb2['h']
                ha2 = ha1_loc
                for i in range(1, ncell_local + 1):
                    fla_i = rhowa_loc * self.fanmax.v * frac[i - 1]
                    flatot_loc = flatot_loc + fla_i
                    airsum_loc = airsum_loc + fla_i * ha2
            # ************************* End of loop for cells *********************

            # Determine temperature of sump, transient of vsump > 0, otherwise steady-state analysis used
            flwo = self.flwi.v - flwl_loc
            if flwo > 0.:
                two_loc = trefw + (self.flwi.v * cpw * (self.twi.v - trefw) - qcells_loc) / flwo / cpw

            tNow_loc = two_loc
            return two_loc, flwl_loc, qcells_loc, flatot_loc, airsum_loc, tNow_loc, wa1_loc, ha1_loc, rhowa_loc

        if controlMode_val < 1:
            # the relative fan speeds are read from the user inputs.
            # Get relative fan flows and associated power requirements
            power_val = 0.
            ncell = 0
            for i in range(1, maxcel_val + 1):
                cell_input = get_cell_input(i)
                if abs(cell_input) > tol:
                    ncell += 1
                    if cell_input > tol:
                        frac[ncell - 1] = cell_input
                        # the fan power is computed based on the user specified relative speed.
                        power_val = power_val + self.pwrmax.v * (frac[ncell - 1])**3.
                        # the actual flow, however, is the maximum of the user-specified relative speed and the
                        # user specified natural convection flow
                        frac[ncell - 1] = max(frac[ncell - 1], (self.fanoff.v / self.fanmax.v))
                    else:
                        frac[ncell - 1] = self.fanoff.v / self.fanmax.v

            # call the subroutine that does most of the air side and all of the water side calculations.
            two, flwl, qcells, flatot, airsum, tNow, wa1, ha1, rhowa = tower_cell_loop()

        else:
            # need to guess fan speeds until one is found that satisfies the user input approach temperature.
            # if there is water flow, find an appropriate fan speed
            if self.flwi.v > 0.:

                # determine how many cells are available (based on the value of the inputs).
                # Note: can only get to this part of the code if controlMode = 1
                ncell = 0
                for i in range(1, maxcel_val + 1):
                    if abs(get_cell_input(i)) > tol:
                        ncell += 1

                # clear out the array that stores the relative fan speed of each cell and the variable
                # that accumulates the power.
                power_val = 0.
                frac[:maxcel_val] = 0.
                step = 0.1

                # first check full speed operation
                converged = False

                # start with the lowest possible fan speed
                fracGuess = self.fanoff.v / self.fanmax.v
                firstStep = float(round((fracGuess + 0.1) * (10**1)) / (10**1)) - fracGuess
                iterCount = 1
                errSign = 0.
                errSignLast = 0.
                error = 0.

                while not converged:
                    power_val = 0.  # 02.15.2017
                    for i in range(1, ncell + 1):
                        frac[i - 1] = fracGuess
                        power_val = power_val + self.pwrmax.v * (frac[i - 1])**3.

                    # call the subroutine that does most of the air side and all of the water side calculations.
                    two, flwl, qcells, flatot, airsum, tNow, wa1, ha1, rhowa = tower_cell_loop()

                    # update the arrays that contain the current guesses and the current results.
                    if iterCount == 1:
                        # the guess for the first iteration was the fan at minimum speed
                        if tNow < self.tWant.v:
                            # the fan can't go any slower so this is the best we can do.
                            converged = True
                        else:
                            error = tNow - self.tWant.v
                            errSign = np.sign(tNow - self.tWant.v)
                            errSignLast = np.sign(tNow - self.tWant.v)
                    else:
                        error = tNow - self.tWant.v
                        errSign = np.sign(tNow - self.tWant.v)

                    # exit if the error is small enough.
                    if abs(error) < 0.01:
                        converged = True

                    # exit if we've already iterated 1000 times.
                    if iterCount >= 1000:
                        # The cooling tower failed to converge on a fan speed that satisfies the desired
                        # water outlet temperature.
                        converged = True

                    # determine a new guessed fraction
                    if errSignLast != errSign:
                        fracGuess = fracGuess - step   # back up one step
                        step = step / 10.              # make the step size smaller
                        errSign = errSignLast          # reset the sign indicator

                    if iterCount == 1:
                        fracGuess = fracGuess + firstStep
                    else:
                        fracGuess = fracGuess + step
                        if fracGuess > 1.:
                            converged = True

                    iterCount += 1

            else:
                # if there is no water flow, turn the fans off.
                power_val = 0.
                ncell = 0
                for i in range(1, maxcel_val + 1):
                    if abs(get_cell_input(i)) > tol:
                        ncell += 1
                        frac[ncell - 1] = 0.

                two, flwl, qcells, flatot, airsum, tNow, wa1, ha1, rhowa = tower_cell_loop()

        # if the sump volume is positive, solve the differential equation to obtain its new temperature.
        if self.vsump.v > 0.:
            ti = self._ti
            a = -(cpw * self.flwi.v + self.uaSump.v) / (cpw * msump)
            flwo = self.flwi.v - flwl
            b_val = (cpw * (flwl * self.tMain.v + flwo * two) + self.qSump.v - self.uaSump.v * self.tdb.v) / (cpw * msump)
            # TODO-NEEDS LIBRARY: solveDiffEq for ODE integration (TRNSYS-specific function)
            # Call solveDiffEq(a, b_val, ti, tf, tsump)
            # Placeholder: Euler integration
            h = self.model.settings.timestep
            # TODO-NEEDS CONVERSION REVIEW: solveDiffEq solves dy/dt = a*y + b analytically in TRNSYS
            tf = ti + h * (a * ti + b_val)
            tsump_val = tf
            self._ti = tf
        else:
            if self.flwi.v < 1.0e-6:
                # the sump volume is negligible and there is no water flow.
                tsump_val = self.twi.v
            else:
                rloss = flwl / self.flwi.v
                tsump_val = (1. - rloss) * two + rloss * self.tMain.v
                qSumpLoss = 0.

        # Determine "mixing cup" exit air conditions and the change in internal energy of the sump
        # (relative to the beginning of the simulation).
        if flatot > 0.:
            haotot = airsum / flatot
            waotot_val = flwl / flatot + wa1
            psydat_mix = humid_air(['T_db', 'T_wb', 'W', 'H'], T_db=taomix_val, T_wb=twbmix_val, W=waotot_val, H=haotot, P=patm)
            taomix_val = psydat_mix['T_db']
            twbmix_val = psydat_mix['T_wb']
            waotot_val = psydat_mix['W']
            haotot = psydat_mix['H']
        else:
            waotot_val = wa1
            taomix_val = self.tdb.v
            twbmix_val = self.twb.v

        deltau_val = msump * cpw * (tsump_val - self.tstart.v)
        if msump <= 0.:
            deltau_val = 0.
        qeff = self.flwi.v * cpw * (self.twi.v - tsump_val)

        # Set the Outputs from this Model
        self.tsump.v = tsump_val
        self.flwi_out.v = self.flwi.v
        self.power.v = power_val
        self.qcells.v = qcells
        self.two.v = two
        self.flwl.v = flwl
        self.taomix.v = taomix_val
        self.twbmix.v = twbmix_val
        self.waotot.v = waotot_val
        self.flatot.v = flatot
        self.deltau.v = deltau_val
        # TODO-NEEDS CONVERSION REVIEW: variable-length cell outputs (11+i); frac(1) written for each active cell
        # self._frac_out[i] = frac[0] if abs(get_cell_input(i)) > tol else 0.
