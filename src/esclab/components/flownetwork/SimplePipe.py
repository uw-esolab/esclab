"""Simple pipe component model (Type 4001)."""

import numpy as np
from eeslib import fluid_properties as fp
from esclab.components.esol_properties import Incompressible as Inc

from esclab.simulate import Component


def FricFactor_IC(Rough, Reynold, guess):
    """
    ***************************************************************************************************
    Friction factor (taken from Piping loss model)
    ***************************************************************************************************
    Uses an iterative method to solve the implicit friction factor function.
    For more on this method, refer to Fox, et al., 2006 Introduction to Fluid Mechanics.
    """
    # Rough is relative roughness [--]

    Acc = .00001

    if(Reynold < 2750.):
        result = 64./np.max([Reynold,1.])
        return

    X = 1/np.sqrt(guess)  # 1. / 0.03
    TestOld = X + 2. * np.log10(Rough / 3.7 + 2.51 * X / Reynold)
    Xold = X
    X = X*0.7 # 1. / (0.03 + 0.005)
    NumTries = 0

    while True:
        NumTries = NumTries + 1
        Test = X + 2 * np.log10(Rough / 3.7 + 2.51 * X / Reynold)
        if (abs(Test - TestOld) <= Acc):
            result = 1. / (X * X)
            break

        if (NumTries > 20):
            print(" Could not find friction factor solution")
            return

        Slope = (Test - TestOld) / (X - Xold)
        Xold = X
        TestOld = Test
        X = max((Slope * X - Test) / Slope,1.e-5)

    return result


class SimplePipe(Component):
    """
    Object: ESOL4001-SimplePipe
    Simulation Studio Model: ESOL4001-SimplePipe


    Author: Matt Tuman
    Editor: Mike Wagner
    Date:     September 19, 2022
    last modified: October 13, 2022
    Ported by: Mike Wagner, February 12, 2026
    """
    # Create types
    Pipe_ID = Component.Parameter()         #pipe inner diameter
    Length_Pipe = Component.Parameter()     #length of pipe
    Roughness = Component.Parameter()       #Pipe Roughness
    fluid = Component.Parameter()           #Fluid going through valve (Each ID is a different SAM fluid (21 = Therminol))
    DELTA_P = Component.Parameter(0.)
    guess = Component.Parameter(28.)        #initial guess for friction factor

    m_dot_in = Component.Input()
    T_in = Component.Input()
    P_in = Component.Input()

    m_dot_out = Component.Output()
    T_out = Component.Output()
    P_out = Component.Output()
    DELTA_P = Component.Output()
    ff = Component.Output()
    _props = Inc()

    def presim_setup(self, **kwargs):
        super().presim_setup(**kwargs)

        #Set the Initial Values of the Outputs (#,Value)
        self.m_dot_out.v = self.m_dot_in.v     # mass flow leaving pipe [kg/s]
        self.T_out.v = self.T_in.v             # temperature leaving the pipe
        self.P_out.v = self.P_in.v             # pressure leaving the pipe
        self.DELTA_P.v = self.DELTA_P.v        # pressure drop through pipe
        self.ff.v = self.guess.v               # Friction factor guess
        return

    def calculate(self):
        super().calculate()

        if(self.m_dot_in.v > 0.00): # okay to continue

            self.m_dot_out.v = self.m_dot_in.v
            self.guess.v =  max([self.ff.v,0.1]) # Friction Factor Guess from last iteration

            # Calculating Pressure Drop
            fluid_name = str(self.fluid.v) if self.fluid.v == self.fluid.v else "Nitrate Salt"
            try:
                rho = self._props.density(fluid_name, self.T_in.v, self.P_in.v)
            except Exception:
                rho = fp.density(fluid_name, T=self.T_in.v, P=self.P_in.v)
            try:
                visc = self._props.viscosity(fluid_name, self.T_in.v, self.P_in.v)
            except Exception:
                visc = fp.viscosity(fluid_name, T=self.T_in.v, P=self.P_in.v)
            vel = self.m_dot_in.v/rho/(3.14/4. *self.Pipe_ID.v**2.)
            Re = rho* vel * self.Pipe_ID.v/visc

            self.ff.v = FricFactor_IC(self.Roughness.v/self.Pipe_ID.v,Re,self.guess.v)
            K_T = (8.*self.ff.v*self.Length_Pipe.v)/((3.14**2.)*(self.Pipe_ID.v**5.)*rho)
            self.DELTA_P.v = K_T * self.m_dot_in.v**2.
            self.P_out.v = self.P_in.v - self.DELTA_P.v

            # Finding change in htf temperature
            self.T_out.v = self.T_in.v     # incompressible fluid

        else: # flow is not possible, calculations will fail

            # self.T_out.v = getOutputValue(2)
            # P_out = getOutputValue(3)

            self.DELTA_P.v = 0.
