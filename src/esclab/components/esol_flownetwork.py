from esclab.simulate import *
from eeslib import fluid_properties as fp
from eeslib.functions import convert
import numpy as np

# ***************************************************************************************************
#  Friction factor (taken from Piping loss model)
# ***************************************************************************************************
#  Uses an iterative method to solve the implicit friction factor function.
#  For more on this method, refer to Fox, et al., 2006 Introduction to Fluid Mechanics.
def FricFactor_IC(Rough, Reynold, guess):
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
# -----------------------------------------------


class SimplePipe(Component):
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

    def setup(self, **kwargs):
        super().setup(**kwargs)

        #Set the Initial Values of the Outputs (#,Value)
        self.m_dot_out.value = self.m_dot_in.value     # mass flow leaving pipe [kg/s]
        self.T_out.value = self.T_in.value             # temperature leaving the pipe
        self.P_out.value = self.P_in.value             # pressure leaving the pipe
        self.DELTA_P.value = self.DELTA_P.value        # pressure drop through pipe
        self.ff.value = self.guess.value               # Friction factor guess
        return 

    def calculate(self):
        super().calculate()

        if(self.m_dot_in.value > 0.00): # okay to continue
            
            self.m_dot_out.value = self.m_dot_in.value
            self.guess.value =  max([self.ff.value,0.1]) # Friction Factor Guess from last iteration
            
            # Calculating Pressure Drop
            rho =  fp.density(self.fluid.value, T = self.T_in.value, P = self.P_in.value)
            visc = fp.viscosity(self.fluid.value, T = self.T_in.value, P = self.P_in.value)
            vel = self.m_dot_in.value/rho/(3.14/4. *self.Pipe_ID.value**2.)
            Re = rho* vel * self.Pipe_ID.value/visc 
            import eeslib
            
            self.ff.value = FricFactor_IC(self.Roughness.value/self.Pipe_ID.value,Re,self.guess.value)
            K_T = (8.*self.ff.value*self.Length_Pipe.value)/((3.14**2.)*(self.Pipe_ID.value**5.)*rho)
            self.DELTA_P.value = K_T * self.m_dot_in.value**2.
            self.P_out.value = self.P_in.value - self.DELTA_P.value
            
            # Finding change in htf temperature
            self.T_out.value = self.T_in.value     # incompressible fluid
            
        else: # flow is not possible, calculations will fail

            # self.T_out.value = getOutputValue(2)
            # P_out = getOutputValue(3)

            self.DELTA_P.value = 0.
            # -----------------------------------------------------------------------------------------------------------------------
        

if __name__ == "__main__":
    
    P = SimplePipe()
    P.Pipe_ID.value = .2
    P.Length_Pipe.value = 100
    P.fluid.value = 'Air'
    P.Roughness.value = 1.e-5

    P.m_dot_in.value = 10
    P.P_in.value = convert('bar','Pa')
    P.T_in.value = 400

    P.setup()
    P.calculate()