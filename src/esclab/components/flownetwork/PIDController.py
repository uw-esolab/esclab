"""PID Controller component model with anti-windup (Type 23)."""

# TRNSYS Type 23: PID Controller with anti-windup
# ----------------------------------------------------------------------------------------------------------------------
#
# This routines implements a Proportional, Integral and Derivative (PID) controller with anti-windup.
#
# The numerical implementation is based on:
# Astrom, K.J. and Wittenmark, B. - Computer controlled Systems, 2nd Edition
# Prentice Hall, Englewood Cliffs, NJ - 1990
# ISBN 0-13-168600-3
#
#
# Inputs
# ----------------------------------------------------------------------------------------------------------------------
# Nb | Variable      | Description                                                    | Input  Units   | Internal Units
# ---|---------------|----------------------------------------------------------------|----------------|----------------
#  1 | ySet          | Setpoint                                                       | any            | any
#  2 | y             | Controlled variable                                            | any            | any
#  3 | onOff         | Controller ON/OFF signal (ON if abs(onOff)>1e-6)               | -              | -
#  4 | uMin          | Minimum value for the control signal                           | any            | any
#  5 | uMax          | Maximum value for the control signal                           | any            | any
#  6 | uThreshold    | Threshold for u. if abs(u)<uThreshold, u = 0                   | any            | any
#  7 | Kc            | Proportional gain                                              | any            | any
#  8 | Ti            | Integral action time constant (0 means no integral action)     | h              | h
#  9 | Td            | Derivative action time constant (0 means no derivative action) | h              | h
# 10 | Tt            | Tracking (anti-windup) time constant (0 means no anti-windup)  | h              | h
# 11 | b             | Fraction of ySet used in the proportional action               | -              | -
# 12 | g             | Fraction of ySet used in the derivative action                 | -              | -
# 13 | N             | High-frequency limit on derivative action                      | -              | -
#
#
# Parameters
# ----------------------------------------------------------------------------------------------------------------------
# Nb | Variable      | Description                                                    | Param. Units   | Internal Units
# ---|---------------|----------------------------------------------------------------|----------------|----------------
#  1 | mode          | 0 = Non-Iterative controller (perform calc. after converg.)    | -              | -
#    |               | 1 = Iterative controller                      |                | -              | -
#  2 | nStick        | Number of iterations after which the controller is "stuck"     | -              | -
#
#
# Outputs
# ----------------------------------------------------------------------------------------------------------------------
# Nb | Variable      | Description                                                    | Output  Units  | Internal Units
# ---|---------------|----------------------------------------------------------------|----------------|----------------
#  1 | u             | Control signal                                                 | any            | any
#  2 | e             | tracking error                                                 | any            | any
#  3 | v             | Unsaturated control signal                                     | any            | any
#  4 | vp            | Proportional action (unsaturated)                              | any            | any
#  5 | vi            | Integral action (unsaturated)                                  | any            | any
#  6 | vd            | Derivative action (unsaturated)                                | any            | any
#  7 | status        | Controller status                                              | -              | -
#
#
# Use of the Storage array
# ----------------------------------------------------------------------------------------------------------------------
#
# Not used. Output array is used for storage (previous call in mode 0, last call previous time step in mode 1)
# out( 8) = ySet  (stored: ySet at end of previous timestep)
# out( 9) = y     (stored: y at end of previous timestep)
# out(10) = u     (stored: u at end of previous timestep)
# out(11) = v     (stored: v at end of previous timestep)
# out(12) = vp    (stored: vp at end of previous timestep)
# out(13) = vi    (stored: vi at end of previous timestep)
# out(14) = vd    (stored: vd at end of previous timestep)
# out(15) = Number of time steps the controller was stuck
#
# ----------------------------------------------------------------------------------------------------------------------
# Copyright 2011 Solar Energy Laboratory, University of Wisconsin-Madison. All rights reserved.
#
# This type was originally written by Michael Kummert at the Fondation Universitaire Luxembourgeoise in 1998.
# It was re-engineered and cleaned in 2004 for TRNSYS 16
#
# Modifications:
# 2006-03-06 - JWT - Changed intent of out(*) to inout
# 2009-06-00 - TPM - conversion to version 17 conventions
# 2012-01-19 - MJD - Placed end of timestep section at end of source code
# 2013-05-31 - TPM - modified the version 17 type to work the same as the version 16 type
# 2015-03-03 - DEB - added SSR code. cleaned up an error message.

import numpy as np

from esclab.simulate import Component


class PIDController(Component):
    """
    PID Controller with anti-windup (TRNSYS Type 23).

    Implements a Proportional, Integral and Derivative (PID) controller with anti-windup
    based on Astrom & Wittenmark (1990).
    """
    trnsys_type = "23"

    # Parameters
    mode = Component.Parameter()    # 0 = Non-Iterative controller, 1 = Iterative controller
    nStick = Component.Parameter()  # Number of iterations after which the controller is "stuck"

    # Inputs
    ySet = Component.Input()        # Setpoint
    y = Component.Input()           # Controlled variable
    onOff = Component.Input()       # Controller ON/OFF signal (ON if abs(onOff)>1e-6)
    uMin = Component.Input()        # Minimum value for the control signal
    uMax = Component.Input()        # Maximum value for the control signal
    uThreshold = Component.Input()  # Threshold for u. if abs(u)<uThreshold, u = 0
    Kc = Component.Input()          # Proportional gain
    Ti = Component.Input()          # Integral action time constant (0 means no integral action) [h]
    Td = Component.Input()          # Derivative action time constant (0 means no derivative action) [h]
    Tt = Component.Input()          # Tracking (anti-windup) time constant (0 means no anti-windup) [h]
    b = Component.Input()           # Fraction of ySet used in the proportional action
    g = Component.Input()           # Fraction of ySet used in the derivative action
    N = Component.Input()           # High-frequency limit on derivative action

    # Outputs (1-7: computed; 8-15: stored state)
    u = Component.Output()          # Control signal
    e = Component.Output()          # Tracking error
    v = Component.Output()          # Unsaturated control signal
    vp = Component.Output()         # Proportional action (unsaturated)
    vi = Component.Output()         # Integral action (unsaturated)
    vd = Component.Output()         # Derivative action (unsaturated)
    status = Component.Output()     # Controller status
    # Stored state outputs (values at end of previous timestep)
    ySet_prev = Component.Output()  # out(8): ySet at end of previous timestep
    y_prev = Component.Output()     # out(9): y at end of previous timestep
    u_prev = Component.Output()     # out(10): u at end of previous timestep
    v_prev = Component.Output()     # out(11): v at end of previous timestep
    vp_prev = Component.Output()    # out(12): vp at end of previous timestep
    vi_prev = Component.Output()    # out(13): vi at end of previous timestep
    vd_prev = Component.Output()    # out(14): vd at end of previous timestep
    nStuckTimeSteps = Component.Output()  # out(15): number of time steps the controller was stuck

    def presim_setup(self, **kwargs):
        super().presim_setup(**kwargs)

        # Sampling rate (TRNSYS time step)
        h = self.model.settings.timestep

        # nStick < 1 is converted to the largest odd number less than (nMaxIterations-5)
        # TODO-NEEDS CONVERSION REVIEW: nMaxIterations not available in esclab; this block may need revisiting
        nStick_val = int(self.nStick.v + 0.01)
        if nStick_val < 1:
            # TODO-NEEDS CONVERSION REVIEW: getnMaxIterations()-4 is TRNSYS-specific
            pass

        # Set the Initial Values of the Outputs
        self.u.v = 0.0
        self.e.v = 0.0
        self.v.v = 0.0
        self.vp.v = 0.0
        self.vi.v = 0.0
        self.vd.v = 0.0
        self.status.v = 0.0
        self.ySet_prev.v = 0.0
        self.y_prev.v = 0.0
        self.u_prev.v = 0.0
        self.v_prev.v = 0.0
        self.vp_prev.v = 0.0
        self.vi_prev.v = 0.0
        self.vd_prev.v = 0.0
        self.nStuckTimeSteps.v = 0.0

    def calculate(self):
        super().calculate()

        # Sampling rate (TRNSYS time step)
        h = self.model.settings.timestep

        mode_val = int(self.mode.v + 0.01)
        nStick_val = int(self.nStick.v + 0.01)

        # nStick < 1 is converted to the largest odd number less than (nMaxIterations-5)
        if nStick_val < 1:
            # TODO-NEEDS CONVERSION REVIEW: getnMaxIterations()-4 is TRNSYS-specific; using a large default
            nStick_val = 999
            if nStick_val % 2 == 0:
                nStick_val = nStick_val - 1

        # --- Read inputs ---
        onOff_val = int(self.onOff.v)

        self.e.v = self.ySet.v - self.y.v

        # Default values for extra parameters
        Tt_val = self.Tt.v
        if Tt_val < 0.0:
            Tt_val = self.Ti.v
        b_val = self.b.v
        if b_val < 0.0:
            b_val = 1.0
        g_val = self.g.v
        if g_val < 0.0:
            g_val = 1.0
        N_val = self.N.v
        if N_val < 0.0:
            N_val = 10.0

        # --- Recall stored values (...1 means at the end of previous time step) ---
        ySet1 = self.ySet_prev.v
        y1 = self.y_prev.v
        u1 = self.u_prev.v
        v1 = self.v_prev.v
        vp1 = self.vp_prev.v
        vi1 = self.vi_prev.v
        vd1 = self.vd_prev.v
        e1 = ySet1 - y1

        # --- Unsaturated control signal ---

        # if the Max. number of iterations has been reached in this time step, keep values constant (not previous time step but just previous call)
        # TODO-NEEDS CONVERSION REVIEW: getTimestepIteration() is TRNSYS-specific; esclab iteration tracking may differ
        if (mode_val == 1) and (self.model.iteration >= nStick_val):
            self.v.v = self.v.v
            self.vp.v = self.vp.v
            self.vi.v = self.vi.v
            self.vd.v = self.vd.v
        else:
            # Proportional
            self.vp.v = self.Kc.v * (b_val * self.ySet.v - self.y.v)
            # Integral
            if self.Ti.v > 0.0:
                # If anti-windup, de-saturate stored value
                if Tt_val > 0.0:
                    self.vi.v = vi1 + h / Tt_val * (u1 - v1)
                # Otherwise start from previous value
                else:
                    self.vi.v = vi1
                # Integral action (backward difference approximation)
                self.vi.v = self.vi.v + h * self.Kc.v / self.Ti.v * self.e.v
            else:
                self.vi.v = 0.0
            # Derivative
            if self.Td.v > 0.0:
                self.vd.v = (self.Td.v / (self.Td.v + N_val * h) * vd1
                             - (self.Kc.v * N_val * self.Td.v) / (self.Td.v + N_val * h)
                             * ((self.y.v - g_val * self.ySet.v) - (y1 - g_val * ySet1)))
            else:
                self.vd.v = 0.0
            self.v.v = self.vp.v + self.vi.v + self.vd.v

        # --- Saturated control signal ---
        if onOff_val > 0:
            self.status.v = 1.0
            self.u.v = self.v.v
            # absolute value too low (e.g. minimum flowrate to start a pump, etc.)
            if abs(self.v.v) < self.uThreshold.v:
                self.u.v = 0.0
                self.status.v = self.status.v + 2.0
            # Special case if uMin > uMax (u is always 0)
            elif self.uMin.v > self.uMax.v:
                self.u.v = 0.0
                self.status.v = self.status.v + 4.0
                self.status.v = self.status.v + 8.0
            # Min and max
            elif self.v.v < self.uMin.v:
                self.u.v = self.uMin.v
                self.status.v = self.status.v + 4.0
            elif self.v.v > self.uMax.v:
                self.u.v = self.uMax.v
                self.status.v = self.status.v + 8.0
        else:
            # Controller is OFF
            self.u.v = 0.0
            self.v.v = 0.0
            self.vp.v = 0.0
            self.vi.v = 0.0
            self.vd.v = 0.0
            self.status.v = 0.0

    def converged(self):
        super().converged()

        # Perform Any "End of Timestep" Manipulations That May Be Required
        self.ySet_prev.v = self.ySet.v
        self.y_prev.v = self.y.v
        self.u_prev.v = self.u.v
        self.v_prev.v = self.v.v
        self.vp_prev.v = self.vp.v
        self.vi_prev.v = self.vi.v
        self.vd_prev.v = self.vd.v
        # Increment the "Stuck" counter
        nStuck = int(self.nStuckTimeSteps.v)
        nStick_val = int(self.nStick.v + 0.01)
        if nStick_val < 1:
            # TODO-NEEDS CONVERSION REVIEW: getnMaxIterations()-4 is TRNSYS-specific; using a large default
            nStick_val = 999
            if nStick_val % 2 == 0:
                nStick_val = nStick_val - 1
        # The max number of iterations in this time step has been reached
        # TODO-NEEDS CONVERSION REVIEW: getTimestepIteration() is TRNSYS-specific; esclab iteration tracking may differ
        if self.model.iteration >= nStick_val + 1:
            nStuck = nStuck + 1
            self.nStuckTimeSteps.v = float(nStuck)
