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
        self.mdot_out.value = self.mdot_in.value
        h_in = fluid_properties.enthalpy(self.fluid, T=self.T_in.value, P=self.P_in.value)
        s_in = fluid_properties.entropy(self.fluid, T=self.T_in.value, P=self.P_in.value)
        self.P_out.value = self.P_in.value / self.model.design.PR
        h_out_s = fluid_properties.enthalpy(self.fluid,P=self.P_out.value, s=s_in)
        self.h_out.value = h_in - (h_in-h_out_s)*self.eta_s.value
        self.T_out.value = fluid_properties.temperature(self.fluid, P=self.P_out.value, h=self.h_out.value)
        self.W.value = (h_in - self.h_out.value)*self.mdot_in.value        
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
        self.mdot_out.value = self.mdot_in.value
        h_in = fluid_properties.enthalpy(self.fluid, T=self.T_in.value, P=self.P_in.value)
        s_in = fluid_properties.entropy(self.fluid, T=self.T_in.value, P=self.P_in.value)
        self.P_out.value = self.P_in.value * self.PR.value
        h_out_s = fluid_properties.enthalpy(self.fluid,P=self.P_out.value, s=s_in)
        self.h_out.value = h_in + (h_out_s - h_in)/self.eta_s.value
        self.T_out.value = fluid_properties.temperature(self.fluid, P=self.P_out.value, h=self.h_out.value)
        self.W.value = (h_in - self.h_out.value)*self.mdot_in.value  
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
    
    def setup(self, **kwargs):
        return super().setup(**kwargs)

    def calculate(self):
        super().calculate()
        self.mdot_out.value = self.mdot_in.value 
        self.P_out.value = self.P_in.value 
        self.T_out.value = self.model.design.T_turb_in
        h_in = fluid_properties.enthalpy(self.model.design.fluid, T=self.T_in.value, P=self.P_in.value)
        self.h_out.value = fluid_properties.enthalpy(self.model.design.fluid, T=self.T_out.value, P=self.P_out.value)
        self.qdot.value = (self.h_out.value - h_in)*self.mdot_in.value
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
        self.mdot_out.value = self.mdot_in.value 
        self.P_out.value = self.P_in.value
        h_in = fluid_properties.enthalpy(self.model.design.fluid, T=self.T_in.value, P=self.P_in.value)
        self.T_out.value = self.T_amb.value + self.dT_approach.value
        self.h_out.value = fluid_properties.enthalpy(self.model.design.fluid, P=self.model.design.P_low, T=self.T_out.value)
        self.qdot.value = (self.h_out.value - h_in)*self.mdot_in.value
        return 

class Summary(Component):
    W_turbine = Component.Input()
    W_compressor = Component.Input()
    Q_combustor = Component.Input()
    eta_cycle = Component.Output()
    
    def calculate(self):
        self.eta_cycle.value = (self.W_turbine.value + self.W_compressor.value)/self.Q_combustor.value

    def converge(self):
        # print(f'Time: {self.model.time}, Eta={self.eta_cycle.value:.4f}, Power={(self.W_turbine.value+self.W_compressor.value)*convert("W","MW"):.2f} MW')
        pass

class Weather(Component):
    def setup(self):
        self.T_amb = Component.Output()
        return 

    def calculate(self):
        super().calculate()
        self.T_amb.value = self.model.design.T_amb + 15*np.sin(self.model.time/5/np.pi)
