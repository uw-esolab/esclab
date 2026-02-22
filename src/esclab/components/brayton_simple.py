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
        self.mdot_out.v = self.mdot_in.v
        h_in = fluid_properties.enthalpy(self.fluid, T=self.T_in.v, P=self.P_in.v)
        s_in = fluid_properties.entropy(self.fluid, T=self.T_in.v, P=self.P_in.v)
        self.P_out.v = self.P_in.v / self.model.design.PR
        h_out_s = fluid_properties.enthalpy(self.fluid,P=self.P_out.v, s=s_in)
        self.h_out.v = h_in - (h_in-h_out_s)*self.eta_s.v
        self.T_out.v = fluid_properties.temperature(self.fluid, P=self.P_out.v, h=self.h_out.v)
        self.W.v = (h_in - self.h_out.v)*self.mdot_in.v        
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
        self.mdot_out.v = self.mdot_in.v
        h_in = fluid_properties.enthalpy(self.fluid, T=self.T_in.v, P=self.P_in.v)
        s_in = fluid_properties.entropy(self.fluid, T=self.T_in.v, P=self.P_in.v)
        self.P_out.v = self.P_in.v * self.PR.v
        h_out_s = fluid_properties.enthalpy(self.fluid,P=self.P_out.v, s=s_in)
        self.h_out.v = h_in + (h_out_s - h_in)/self.eta_s.v
        self.T_out.v = fluid_properties.temperature(self.fluid, P=self.P_out.v, h=self.h_out.v)
        self.W.v = (h_in - self.h_out.v)*self.mdot_in.v  
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
        self.mdot_out.v = self.mdot_in.v 
        self.P_out.v = self.P_in.v 
        self.T_out.v = self.model.design.T_turb_in
        h_in = fluid_properties.enthalpy(self.model.design.fluid, T=self.T_in.v, P=self.P_in.v)
        self.h_out.v = fluid_properties.enthalpy(self.model.design.fluid, T=self.T_out.v, P=self.P_out.v)
        self.qdot.v = (self.h_out.v - h_in)*self.mdot_in.v
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
        self.mdot_out.v = self.mdot_in.v 
        self.P_out.v = self.P_in.v
        h_in = fluid_properties.enthalpy(self.model.design.fluid, T=self.T_in.v, P=self.P_in.v)
        self.T_out.v = self.T_amb.v + self.dT_approach.v
        self.h_out.v = fluid_properties.enthalpy(self.model.design.fluid, P=self.model.design.P_low, T=self.T_out.v)
        self.qdot.v = (self.h_out.v - h_in)*self.mdot_in.v
        return 

class Summary(Component):
    W_turbine = Component.Input()
    W_compressor = Component.Input()
    Q_combustor = Component.Input()
    eta_cycle = Component.Output()
    
    def calculate(self):
        self.eta_cycle.v = (self.W_turbine.v + self.W_compressor.v)/self.Q_combustor.v

class Weather(Component):
    def presim_setup(self):
        self.T_amb = Component.Output()
        return 

    def calculate(self):
        super().calculate()
        self.T_amb.v = self.model.design.T_amb + 15*np.sin(self.model.time/5/np.pi)
