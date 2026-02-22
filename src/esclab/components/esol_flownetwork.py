from esclab.simulate import *
from eeslib import fluid_properties as fp
from eeslib.functions import convert
import numpy as np
from esclab.components.esol_properties import Incompressible as Inc

# 23
# 57
# 65
# 162
# 603
# 4001

# 4006
# 4007
# 4008
# 6001
# 6003
# 6007
# 6011
# 6014
# 6016
# 6017
# 6019
# 6022
# 6027
# 6028

class TeeOut(Component):
    #    PARAMETERS
    Fluid_ID = Component.Parameter()           # Fluid ID: 40 for dowtherm A
    Solving_Method = Component.Parameter()     # Adjust solver

    #    INPUTS
    m_dot = Component.Parameter()              # Mass flow into tee [kg/s]
    Pressure = Component.Parameter()           # Pressure at inlet of the tee [Pa]
    Temperature = Component.Parameter()        # Temperature of fluid into tee [C]
    f = Component.Parameter()                  # Fraction of flow of tee [0-1]
    P_1_out = Component.Parameter()            # Pressure at the end of loop 1 [Pa]
    P_2_out = Component.Parameter()            # Pressure at the end of loop 2 [Pa]
    m_dot_prev = Component.Parameter()         # Previous mass flow rate into the tee [kg/s]
    P_prev = Component.Parameter()             # Previous pressure at the inlet of the tee [Pa]
    error = Component.Parameter()              # Error in pressure due do a component with specified pressure
    Mass_Counter = Component.Parameter()       # Total amount of mass in the system up to this point

    #    OUTPUTS
    m_dot_1 = Component.Output()
    P_1 = Component.Output()
    Temp_1 = Component.Output()
    m_dot_2 = Component.Output()
    P_2 = Component.Output()
    Temp_2 = Component.Output()
    f_out = Component.Output()
    m_dot_out = Component.Output()
    Pressure_out = Component.Output()
    Mass_Counter_out = Component.Output()
    

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
    """
    Object: ESOL4001-SimplePipe
    Simulation Studio Model: ESOL4001-SimplePipe
    

    Author: Matt Tuman
    Editor: Mike Wagner
    Date:	 September 19, 2022
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
            rho =  fp.density(self.fluid.v, T = self.T_in.v, P = self.P_in.v)
            visc = fp.viscosity(self.fluid.v, T = self.T_in.v, P = self.P_in.v)
            vel = self.m_dot_in.v/rho/(3.14/4. *self.Pipe_ID.v**2.)
            Re = rho* vel * self.Pipe_ID.v/visc 
            import eeslib
            
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
            # -----------------------------------------------------------------------------------------------------------------------
        

# -----------------------------------------------------------------------------------------------
class VarSpeedPump(Component):
    """
    Docstring for VarSpeedPump
    Subroutine Type4004
    Object: ESOL4004-VarPump
    
    Simulation Studio Model: ESOL4004-VarPump
    
    Author: Matt Tuman
    Editor: Mike Wagner
    Date:	 January 05, 2023
    last modified: January 05, 2023
    Converted by: Mike Wagner, February 12, 2026
    """

    #     PARAMETERS
    N_pumps_parallel = Component.Parameter()         #  Number of pumps that are in parallel
    Pump_a0 = Component.Parameter()                  #  head = Pump_a0 + Pump_a1*Q_dot + Pump_a2*Q_dot^2 [-]
    Pump_a1 = Component.Parameter()                  #  head = Pump_a0 + Pump_a1*Q_dot + Pump_a2*Q_dot^2 [1/m^3/s]
    Pump_a2 = Component.Parameter()                  #  head = Pump_a0 + Pump_a1*Q_dot + Pump_a2*Q_dot^2 [1/m^3/s]
    NPSH_a0 = Component.Parameter()                  #  NPSHr = Pump_a0 + Pump_a1*Q_dot + Pump_a2*Q_dot^2 [-]
    NPSH_a1 = Component.Parameter()                  #  NPSHr = Pump_a0 + Pump_a1*Q_dot + Pump_a2*Q_dot^2 [1/m^3/s]
    NPSH_a2 = Component.Parameter()                  #  NPSHr = Pump_a0 + Pump_a1*Q_dot + Pump_a2*Q_dot^2 [1/m^3/s]
    RPM_full = Component.Parameter()                 #  RPM of pump at 100% speed 
    LR = Component.Parameter()                       #  Learning rate: Close to 1 means that the mass flow can change significantly [0-1]
    D_outlet = Component.Parameter()                 #  Diameter of the pump outlet [-]
    Fluid_ID = Component.Parameter()                 #  Fluid ID
    Pump_Solver = Component.Parameter()              #  Determines if hydraulic solver is used to update mass flow

    #     INPUTS
    Mass_Flow = Component.Input()                #  Mass flow into the pump [kg/s]
    Pressure = Component.Input()                 #  Pressure of fluid at inlet of the pump [Pa]
    Temperature = Component.Input()              #  Temperature of fluid at inlet of pump [C]
    speed = Component.Input()                    #  Speed that the pump is operating at [0-1]
    error = Component.Input()                    #  Error accumulated within pressure drops throughout the system [Pa]
    mass_count = Component.Input()               #  Total mass counted before the pump [kg]

    #      Outputs
    m_dot_out = Component.Output() #  Mass Flow
    P_out = Component.Output() #  Pressure
    Temperature_out = Component.Output() #  Temperature
    mass_count_out = Component.Output()
    cavitation = Component.Output()
    g = 9.81 #  gravity

    # -----------------------------------------------------------------------------------------------------------------------
    def presim_setup(self, **kwargs):
        # Set the Initial Values of the Outputs (#,Value)
        self.m_dot_out.v = self.Mass_Flow.v
        self.P_out.v = self.Pressure.v
        self.Temperature_out.v = self.Temperature.v
        self.mass_count_out.v = self.mass_count.v
        self.cavitation.v = 0    #  Not cavitating initially
    # -----------------------------------------------------------------------------------------------------------------------

    # -----------------------------------------------------------------------------------------------------------------------
    def calculate(self):
        # -----------------------------------------------------------------------------------------------------------------------
        # Post-convergence
        if self.model.is_converged:
            # Once model has converged, check if pump is likely cavitating
            rho = Inc.density(self.Fluid_ID.v, self.Temperature.v, 0.0)
            Q_dot = self.Mass_Flow.v/rho/self.N_pumps_parallel.v
            NPSHr = self.speed.v**2 * (self.NPSH_a0.v + self.NPSH_a1.v*(Q_dot/self.speed.v) + self.NPSH_a2.v*(Q_dot/self.speed.v)**2) 
            NPSH_meas = self.Pressure.v/rho/self.g - 1034000.0/rho/self.g #  NPSH in the simulation is relative to the pressure in the expansion tank (reason for subtracting 150 psi)
            if(NPSHr>NPSH_meas):
                self.cavitation.v = 1.0
            else:
                self.cavitation.v = 0.0
            return


        # -----------------------------------------------------------------------------------------------------------------------
        if (self.model.iteration == 0) or (self.model.iteration == 0 and self.model.time == self.model.settings.timestep):
            #  Set output values to the computed values from the previous timestep (don't want any manipulation without feedback)
            #  Determine pressure increase in pump
            rho = Inc.density(self.Fluid_ID.v, self.Temperature.v, 0.0)
            delta_P = self.g*self.speed.v**2*(
                rho*self.Pump_a0.v + 
                self.Pump_a1.v/self.speed.v*(self.Mass_Flow.v/self.N_pumps_parallel.v) + 
                self.Pump_a2.v/rho/self.speed.v**2*(self.Mass_Flow.v/self.N_pumps_parallel.v)**2
                )

            # Set the Initial Values of the Outputs (#,Value)
            self.m_dot_out.v = self.Mass_Flow.v
            self.P_out.v = self.Pressure.v+delta_P
            self.Temperature_out.v = self.Temperature.v
            self.mass_count_out.v = self.mass_count.v 
            # self.cavitation.v #don't change

        else:
            if(self.Pump_Solver.v == 0.0):
                #  Compute Density of fluid
                rho = Inc.density(self.Fluid_ID.v, self.Temperature.v, 0.0)

                #  Compute head loss in system
                delta_P = self.P_out.v-(self.Pressure.v+self.error.v)
                H_L = delta_P/rho/9.81 * 3.28084 #  [ft]

                #  Compute new flow rate corresponding to head loss in system
                #  (solving a quadratic equation for the pump curve fit)
                A = self.Pump_a2.v
                B = self.Pump_a1.v*self.speed.v
                D = self.Pump_a0.v*(self.speed.v**2) - H_L

                discriminant = (B)**2 - 4*A*D
                if (discriminant>=0):
                    sol1 = (-B + np.sqrt(discriminant))/(2*A)
                    sol2 = (-B - np.sqrt(discriminant))/(2*A)
                else:
                    discriminant = 0
                    sol1 = (-B + np.sqrt(discriminant))/(2*A)
                    sol2 = (-B - np.sqrt(discriminant))/(2*A)

                Q_dot_new = max([sol1, sol2]) * 0.00006309019640343866 #  [m^3/s]
                m_dot_new = Q_dot_new*rho

                #  Update mass flow rate according to learning rate
                m_dot_adj = self.LR.v*m_dot_new + (1-self.LR.v)*self.m_dot_out.v/self.N_pumps_parallel.v

                P_out = (self.Pressure.v+delta_P)
                
                # Set the Outputs from this Model (#,Value)
                self.m_dot_out.v = m_dot_adj*self.N_pumps_parallel.v
                self.P_out.v = P_out
                self.Temperature_out.v = self.Temperature.v
                self.mass_count_out.v = self.mass_count.v
                # self.cavitation.v # don't change

            else:
                #  Determine pressure increase in pump
                rho = Inc.density(  self.Fluid_ID.v, self.Temperature.v, 0.0)
                delta_P = self.g*self.speed.v**2*(rho*self.Pump_a0.v + 
                    self.Pump_a1.v/self.speed.v*(self.Mass_Flow.v/self.N_pumps_parallel.v) + 
                    self.Pump_a2.v/rho/self.speed.v**2*(self.Mass_Flow.v/self.N_pumps_parallel.v)**2)

                # Set the Outputs from this Model (#,Value)
                self.m_dot_out.v = self.Mass_Flow.v
                self.P_out.v = self.Pressure.v + delta_P
                self.Temperature_out.v = self.Temperature.v
                self.mass_count_out.v = self.mass_count.v
                # self.cavitation.v # don't change
    # -----------------------------------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------






if __name__ == "__main__":
    
    P = SimplePipe()
    P.Pipe_ID.v = .2
    P.Length_Pipe.v = 100
    P.fluid.v = 'Air'
    P.Roughness.v = 1.e-5

    P.m_dot_in.v = 10
    P.P_in.v = convert('bar','Pa')
    P.T_in.v = 400

    P.presim_setup()
    P.calculate()