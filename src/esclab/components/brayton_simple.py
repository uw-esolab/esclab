from esclab.simulate import *
from eeslib import fluid_properties
import numpy as np

# -------------------------------------------
class Turbine(Component):

    eta_s =    Component.Parameter(0.9)
    mdot_in =  Component.Input()
    T_in    =  Component.Input()
    P_in    =  Component.Input()
    h_out =    Component.Output()
    T_out =    Component.Output()
    P_out =    Component.Output()
    mdot_out = Component.Output()
    W =        Component.Output()
    fluid = 'Air'

    def calculate(self):
        self.mdot_out = self.mdot_in
        h_in = fluid_properties.enthalpy(self.fluid, T=self.T_in, P=self.P_in)
        s_in = fluid_properties.entropy(self.fluid, T=self.T_in, P=self.P_in)
        self.P_out = self.P_in / self.model.design.PR
        h_out_s = fluid_properties.enthalpy(self.fluid,P=self.P_out, s=s_in)
        self.h_out = h_in - (h_in-h_out_s)*self.eta_s
        self.T_out = fluid_properties.temperature(self.fluid, P=self.P_out, h=self.h_out)
        self.W = (h_in - self.h_out)*self.mdot_in        
        return 
    
# -------------------------------------------
class Compressor(Component):
    eta_s = Component.Parameter(0.85)
    PR = Component.Parameter(3.)

    mdot_in = Component.Input()
    T_in    = Component.Input()
    P_in    = Component.Input()

    h_out = Component.Output()
    T_out = Component.Output()
    P_out = Component.Output()
    mdot_out = Component.Output()
    W = Component.Output()

    fluid = 'Air'
        
    def calculate(self):
        self.mdot_out = self.mdot_in
        h_in = fluid_properties.enthalpy(self.fluid, T=self.T_in, P=self.P_in)
        s_in = fluid_properties.entropy(self.fluid, T=self.T_in, P=self.P_in)
        self.P_out = self.P_in * self.PR
        h_out_s = fluid_properties.enthalpy(self.fluid,P=self.P_out, s=s_in)
        self.h_out = h_in + (h_out_s - h_in)/self.eta_s
        self.T_out = fluid_properties.temperature(self.fluid, P=self.P_out, h=self.h_out)
        self.W = (h_in - self.h_out)*self.mdot_in  
        return 
    
# -------------------------------------------
class Combustor(Component):
    mdot_in = Component.Input()
    T_in    = Component.Input()
    P_in    = Component.Input()

    h_out = Component.Output()
    T_out = Component.Output()
    P_out = Component.Output()
    mdot_out = Component.Output()
    qdot = Component.Output()
        
    def __init__(self):
        return 
    
    def presim_setup(self, **kwargs):
        return super().presim_setup(**kwargs)

    def calculate(self):
        super().calculate()
        self.mdot_out = self.mdot_in 
        self.P_out = self.P_in 
        self.T_out = self.model.design.T_turb_in
        h_in = fluid_properties.enthalpy(self.model.design.fluid, T=self.T_in, P=self.P_in)
        self.h_out = fluid_properties.enthalpy(self.model.design.fluid, T=self.T_out, P=self.P_out)
        self.qdot = (self.h_out - h_in)*self.mdot_in
        return 
    
# -------------------------------------------
class Cooler(Component):
    dT_approach = Component.Parameter(20)

    mdot_in    = Component.Input()
    T_in       = Component.Input()
    P_in       = Component.Input()
    T_amb      = Component.Input()

    h_out = Component.Output()
    T_out = Component.Output()
    P_out = Component.Output()
    mdot_out = Component.Output()
    qdot = Component.Output()
    
    def calculate(self):
        self.mdot_out = self.mdot_in
        self.P_out = self.P_in
        h_in = fluid_properties.enthalpy(self.model.design.fluid, T=self.T_in, P=self.P_in)
        self.T_out = self.T_amb + self.dT_approach
        self.h_out = fluid_properties.enthalpy(self.model.design.fluid, P=self.model.design.P_low, T=self.T_out)
        self.qdot = (self.h_out - h_in)*self.mdot_in
        return 

class Summary(Component):
    W_turbine = Component.Input()
    W_compressor = Component.Input()
    Q_combustor = Component.Input()
    eta_cycle = Component.Output()
    
    def calculate(self):
        self.eta_cycle = (self.W_turbine + self.W_compressor)/self.Q_combustor

class Weather(Component):
    def presim_setup(self):
        self.T_amb = Component.Output()
        return 

    def calculate(self):
        super().calculate()
        self.T_amb = self.model.design.T_amb + 15*np.sin(self.model.time/1e5/np.pi) + 15*np.sin(self.model.time/1e5/np.pi*2+.5)
