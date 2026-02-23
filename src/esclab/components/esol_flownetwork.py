"""
Provides:
-------------------------------------------------------------------------------------------
SimplePipe      | A simple pipe model that calculates pressure drop based on the 
Type4001        | Darcy-Weisbach equation and updates the outlet pressure and temperature
                | accordingly. The friction factor is calculated using an iterative 
                | method to solve the implicit function for the friction factor.
-------------------------------------------------------------------------------------------
VarSpeedPump    | A variable speed pump model that calculates the pressure increase 
Type4004        | based on a pump curve fit and updates the mass flow rate according to 
                | a learning rate. The model also checks for cavitation based on the NPSH
                | of the pump and the measured NPSH in the system.
-------------------------------------------------------------------------------------------
TeeOut          | A tee junction model that splits the mass flow into two branches based 
Type4006        | on a specified fraction, and updates the pressures and temperatures 
                | at the outlets accordingly. The model also accounts for pressure drops 
                | in the system and updates the mass flow rate based on feedback from 
                | the previous iteration.
-------------------------------------------------------------------------------------------
"""


from esclab.simulate import *
from eeslib import fluid_properties as fp
from eeslib.functions import convert
import numpy as np
from scipy.interpolate import RectBivariateSpline
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
    """
    # Object: 4006-Tee_Out
    # Simulation Studio Model: ESOL4006-Tee_Out
    # 

    # Author: Matt Tuman
    # Editor: 
    # Date:	 October 20, 2022
    # last modified: October 24, 2022
    """

    #    PARAMETERS
    Fluid_ID = Component.Parameter()           # Fluid ID: 40 for dowtherm A
    Solving_Method = Component.Parameter()     # Adjust solver

    #    INPUTS
    m_dot = Component.Input()              # Mass flow into tee [kg/s]
    Pressure = Component.Input()           # Pressure at inlet of the tee [Pa]
    Temperature = Component.Input()        # Temperature of fluid into tee [C]
    f = Component.Input()                  # Fraction of flow of tee [0-1]
    P_1_out = Component.Input()            # Pressure at the end of loop 1 [Pa]
    P_2_out = Component.Input()            # Pressure at the end of loop 2 [Pa]
    m_dot_prev = Component.Input()         # Previous mass flow rate into the tee [kg/s]
    P_prev = Component.Input()             # Previous pressure at the inlet of the tee [Pa]
    error = Component.Input()              # Error in pressure due do a component with specified pressure
    Mass_Counter = Component.Input()       # Total amount of mass in the system up to this point

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

    def presim_setup(self, **kwargs):
        pass 

    def calculate(self):

        # Initialize output values for the first time step
        if self.model.is_first_step:
            # Use initial guess of fraction of flow to compute mass flow rates out of tee
            self.m_dot_1.v = self.m_dot.v*self.f.v
            self.m_dot_2.v = self.m_dot.v - self.m_dot_1.v

            # Set the Initial Values of the Outputs (#,Value)
            self.P_1.v = self.Pressure.v                   # Pressure sent down branch 1
            self.Temp_1.v = self.Temperature.v             # temperature sent down branch 1
            self.P_2.v = self.Pressure.v                   # pressure sent down branch 2
            self.Temp_2.v = self.Temperature.v             # temperature sent down branch 2
            self.f_out.v = self.f.v                        # previous flow fraction guessed
            self.m_dot_out.v = self.m_dot.v                # mass flow rate that entered the tee
            self.Pressure_out.v = self.Pressure.v          # pressure that entered the tee 
            self.Mass_Counter_out.v = self.Mass_Counter.v  # mass counter output used for the expansion system type
            return 
        
        #-----------------------------------------------------------------------------------------------------------------------
        # -----------  iteration calculations ----------------------
        #-----------------------------------------------------------------------------------------------------------------------

        #-------------------- SOLVE USING LINEARIZED SOLUTION-------------------------------------------------------------
        #
        if (self.Solving_Method.v == 0.0):
            if (self.model.iteration > 0) or (self.model.time > self.model.settings.timestep):
                # Compute K values
                self.m_dot_1.v = self.f.v*self.m_dot_prev.v
                self.m_dot_2.v = (1-self.f.v)*self.m_dot_prev.v
                K_T1 = (self.P_prev.v - self.P_1_out.v)/(self.m_dot_1.v**2)
                K_T2 = (self.P_prev.v - self.P_2_out.v)/(self.m_dot_2.v**2)

                # Compute new f value
                f_new = K_T2*self.m_dot_2.v/(K_T2*self.m_dot_2.v + K_T1*self.m_dot_1.v)

                # Update f value 
                #f = (f_new + f)/2.0
                LR = 0.5
                self.f.v = f_new*LR + self.f.v*(1-LR)
            
            # Set Outputs
            m_dot_1_new = self.f.v*self.m_dot.v
            m_dot_2_new = (1-self.f.v)*self.m_dot.v
            f_new = self.f.v

        #-----------------------------------------------------------------------------------------------------------------------
        #
        #-------------------- SOLVE USING DIRECT QUADRATIC SOLUTION-------------------------------------------------------------
        #
        elif (self.Solving_Method.v == 1.0):
            # Set learning rate
            LR = .41

            ### Compute k_t values using output from the previous iteration
            error_term = 3.0
            if(error_term == 1.0):
                K_T1 = (self.P_prev.v - (self.P_1_out.v+self.error.v)) / (self.m_dot_1.v**2)
                K_T2 = (self.P_prev.v - self.P_2_out.v) / (self.m_dot_2.v**2)
            elif(error_term == 2.0):
                K_T1 = (self.P_prev.v - self.P_1_out.v) / (self.m_dot_1.v**2)
                K_T2 = (self.P_prev.v - (self.P_2_out.v+self.error.v)) / (self.m_dot_2.v**2)
            else:
                K_T1 = (self.P_prev.v - self.P_1_out.v) / (self.m_dot_1.v**2)
                K_T2 = (self.P_prev.v - self.P_2_out.v) / (self.m_dot_2.v**2)

        #-----------------------------------------------------------------------------------------------------------------------
        #
        #--------------------------- GRADIENT DESCENT APPROACH -----------------------------------------------------------------
        ############################ WOULD NOT RECCOMEND USING ###############################

        elif (self.Solving_Method.v == 2.0):
            Beta = 0.01
            # First iteration: allow to solve for new mass flow rates
            if self.model.is_first_iteration:
                self.m_dot_1.v = self.f.v*self.m_dot.v
                self.m_dot_2.v = (1-self.f.v)*self.m_dot.v
                f_new = self.f.v

            else:
                # Second iteration or more
                # Load in the fraction, inlet pressure, and mass flow rate from the first iteration
                m_dot_1_new = self.f.v*self.m_dot_prev.v
                m_dot_2_new = (1-self.f.v)*self.m_dot_prev.v
                
                # Compute k_t values using output from the first iteration
                if (m_dot_1_new <= 0) :
                    K_T1 = (self.P_prev.v - self.P_1_out.v) / (.00000001**2)
                else:
                    K_T1 = (self.P_prev.v - self.P_1_out.v) / (m_dot_1_new**2)
                
                if (m_dot_2_new <= 0) :
                    K_T2 = (self.P_prev.v - self.P_2_out.v) / (.000000001**2)
                else:
                    K_T2 = (self.P_prev.v - self.P_2_out.v) / (m_dot_2_new**2)

                # Compute the new fraction value: f_adj
                dPdf = ((K_T2 - K_T1)*self.f.v - K_T2) * ((K_T2-K_T1)*self.f.v**2 - 2*K_T2*self.f.v + K_T2)
                f_adj = self.f.v - Beta*dPdf

                if (self.P_1_out.v < 0) :
                    f_adj = .0000000001
                if (self.P_2_out.v < 0) :
                    f_adj = .9999999999
                
                f_new = LR*f_adj + (1-LR)*self.f_out.v
                m_dot_1_new = self.m_dot.v*f_new
                m_dot_2_new = self.m_dot.v*(1-f_new)
            

        # IF USING ESOL4050-Parallel-Flow-Solver
        elif (self.Solving_Method.v == 3.0 ):  
            # Read in flow value computed by ESOL4050
            f_new = self.f.v
            #f_new_adj = f  # Added by C.Volkwein

            m_dot_1_new = f_new*self.m_dot.v
            m_dot_2_new = (1-f_new)*self.m_dot.v

        #-----------------------------------------------------------------------------------------------------------------------
        #Set the Outputs from this Model (#,Value)

        # Set output values
        self.m_dot_1.v = m_dot_1_new
        self.P_1.v = self.Pressure.v
        self.Temp_1.v = self.Temperature.v
        self.m_dot_2.v = m_dot_2_new
        self.P_2.v = self.Pressure.v
        self.Temp_2.v = self.Temperature.v
        self.f_out.v = f_new
        self.m_dot_out.v = self.m_dot.v
        self.Pressure_out.v = self.Pressure.v
        self.Mass_Counter_out.v = self.Mass_Counter.v
        return 

# -----------------------------------------------------------------------------------------------------------------------

class Valve(Component):
    """
    # Object: ESOL4007-Valve
    # Simulation Studio Model: ESOL4007-Valve
    # 

    # Author: Matt Tuman
    # Editor: 
    # Date:	 October 27, 2022
    # last modified: October 27, 2022
    """

    #    PARAMETERS
    Diameter = Component.Parameter()
    Fluid_ID = Component.Parameter()
    Valve_Type = Component.Parameter()

    #    INPUTS
    m_dot = Component.Input()
    Pressure = Component.Input()
    Temperature = Component.Input()
    fraction_open = Component.Input()
    mass_counter = Component.Input()
    Cv = Component.Input()

    #    OUTPUTS
    m_dot_out = Component.Output()
    Pressure_out = Component.Output()
    Temperature_out = Component.Output()
    mass_counter_out = Component.Output()
    Cv_out = Component.Output()
    VP_output = Component.Output()

    # stored data
    fraction_open = float('nan')

    def calculate(self):

        if self.model.is_first_step:
            ## COMPUTE PRESSURE DROP
            ##########################
            # Compute Volumetric Flow Rate
            rho_fluid = Inc.density(self.Fluid_ID.v, self.Temperature.v, 0.0)
            Q = self.m_dot.v/rho_fluid
            # Convert flowrate to gpm
            Q = Q*15850.323140625002 
            # Compute specific gravity of fluid 
            SG = rho_fluid/1000.0
            # Compute Cv of valve
            Cv = CV_data(self.Valve_Type.v, self.Diameter.v, self.fraction_open)
            # Compute pressure drop [Pa]
            dP = SG*Q**2/(Cv**2) * 6894.76
            
            #Set the Initial Values of the Outputs (#,Value)
            self.m_dot_out.v = self.m_dot.v     #1
            self.Pressure_out.v = self.Pressure.v - dP      #2
            self.Temperature_out.v = self.Temperature.v     #3
            self.mass_counter_out.v = self.mass_counter.v       #4
            self.Cv_out.v = Cv      #5
            self.VP_output.v = self.fraction_open        #6

            # Set Initial Values of Dynamic Storage
            # Call SetDynamicArrayValueThisIteration(1, fraction_open) # Store valve position

            return


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


def CV_data(Valve_Type, D_in, Valve_position):
    """
    Compute C_v for various valve types using 2D interpolation.
    
    Parameters
    ----------
    Valve_Type : int
        1 = Concentric Butterfly Valve
        2 = Triple Offset Butterfly Valve
    D_in : float
        Valve diameter [m]
    Valve_position : float
        Valve opening position [0-1], where 1 is fully open
    
    Returns
    -------
    float
        Flow coefficient C_v
    
    References
    ----------
    https://www.valvesonline.com.au/references/flow-rates/butterfly-valves/
    """
    C_v_min = 0.0001  # minimum CV allowed
    
    # Return minimum flow to prevent hydraulic model crashes
    if Valve_position == 0:
        return C_v_min
    
    D = D_in * 39.3701  # Convert m -> inches
    
    # Define valve data: diameters, positions, and CV values
    valve_data = {
        1: {  # Concentric Butterfly
            'D': np.array([4.0, 8.0, 12.0, 16.0, 20.0, 24.0, 32.0]),
            'Pos': np.array([0.0, 1/9, 2/9, 3/9, 4/9, 5/9, 6/9, 7/9, 8/9, 1.0]),
            'CV': np.array([
                [0.0, 0.5, 17.0, 36.0, 78.0, 139.0, 230.0, 364.0, 546.0, 600.0],
                [0.0, 3.0, 89.0, 188.0, 408.0, 727.0, 1202.0, 1903.0, 2854.0, 3136.0],
                [0.0, 5.0, 234.0, 495.0, 1072.0, 1911.0, 3162.0, 5005.0, 7507.0, 8250.0],
                [0.0, 8.0, 464.0, 983.0, 2130.0, 3797.0, 6282.0, 9942.0, 14913.0, 16388.0],
                [0.0, 14.0, 791.0, 1674.0, 3628.0, 6465.0, 10698.0, 16931.0, 25396.0, 27908.0],
                [0.0, 22.0, 1222.0, 2587.0, 5605.0, 9989.0, 16528.0, 26157.0, 39236.0, 43116.0],
                [0.0, 45.0, 2387.0, 4791.0, 8736.0, 13788.0, 20613.0, 31395.0, 48117.0, 68250.0],
            ])
        },
        2: {  # Triple Offset Butterfly
            'D': np.array([4.0, 8.0, 12.0, 16.0, 20.0, 24.0, 32.0]),
            'Pos': np.array([0.0, 1/9, 2/9, 3/9, 4/9, 5/9, 6/9, 7/9, 8/9, 1.0]),
            'CV': np.array([
                [0.0, 8.4, 29.3, 58.6, 92.0, 140.0, 200.0, 330.0, 370.0, 420.0],
                [0.0, 38.2, 140.0, 270.0, 420.0, 640.0, 900.0, 1500.0, 1690.0, 1920.0],
                [0.0, 88.4, 310.0, 620.0, 980.0, 1460.0, 2080.0, 3450.0, 3890.0, 4420.0],
                [0.0, 150.0, 530.0, 1060.0, 1660.0, 2490.0, 3540.0, 5870.0, 6620.0, 7520.0],
                [0.0, 270.0, 930.0, 1850.0, 2900.0, 4350.0, 6190.0, 10300.0, 11600.0, 13200.0],
                [0.0, 420.0, 1450.0, 2890.0, 4530.0, 6800.0, 9680.0, 16100.0, 18100.0, 20600.0],
                [0.0, 810.0, 2820.0, 5630.0, 8840.0, 13300.0, 18900.0, 31400.0, 35400.0, 40200.0],
            ])
        }
    }
    
    if Valve_Type not in valve_data:
        raise ValueError(f"Valve_Type {Valve_Type} not supported (1=Concentric, 2=Triple Offset)")
    
    data = valve_data[Valve_Type]
    D_vals = data['D']
    Pos_vals = data['Pos']
    CV_vals = data['CV']
    
    # Create 2D bivariate spline interpolator (cubic by default)
    # RectBivariateSpline expects data with x as columns, y as rows
    spl = RectBivariateSpline(D_vals, Pos_vals, CV_vals, kx=1, ky=1)
    
    # Evaluate at the requested diameter and position
    C_v = float(spl(D, Valve_position)[0, 0])
    
    # Ensure minimum CV
    return max(C_v, C_v_min)



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