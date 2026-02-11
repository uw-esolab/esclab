from esclab.simulate import *
from eeslib import fluid_properties
import numpy as np

# -------------------------------------------
class Turbine(Component):
    def __init__(self):
        self.eta_s = Component.Parameter(0.9)

        self.mdot_in = Component.Input()
        self.T_in    = Component.Input()
        self.P_in    = Component.Input()

        self.h_out = Component.Output()
        self.T_out = Component.Output()
        self.P_out = Component.Output()
        self.mdot_out = Component.Output()
        self.W = Component.Output()
        
        self.fluid = 'Air'
        return 
    
    def setup(self):
        pass
    
    def calculate(self):
        super().calculate()

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
    def __init__(self):
        self.eta_s = Component.Parameter(0.85)
        self.PR = Component.Parameter(3.)

        self.mdot_in = Component.Input()
        self.T_in    = Component.Input()
        self.P_in    = Component.Input()

        self.h_out = Component.Output()
        self.T_out = Component.Output()
        self.P_out = Component.Output()
        self.mdot_out = Component.Output()
        self.W = Component.Output()

        self.fluid = 'Air'
        
        return 
    
    def setup(self, **kwargs):
        return super().setup(**kwargs)

    def calculate(self):
        super().calculate()
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
    def __init__(self):
        self.mdot_in = Component.Input()
        self.T_in    = Component.Input()
        self.P_in    = Component.Input()

        self.h_out = Component.Output()
        self.T_out = Component.Output()
        self.P_out = Component.Output()
        self.mdot_out = Component.Output()
        self.qdot = Component.Output()
        
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
    def __init__(self):
        self.dT_approach = Component.Parameter(20)

        self.mdot_in    = Component.Input()
        self.T_in       = Component.Input()
        self.P_in       = Component.Input()
        self.T_amb      = Component.Input()

        self.h_out = Component.Output()
        self.T_out = Component.Output()
        self.P_out = Component.Output()
        self.mdot_out = Component.Output()
        self.qdot = Component.Output()
        super().__init__()

    def setup(self):
        return 

    def calculate(self):
        super().calculate()
        self.mdot_out.value = self.mdot_in.value 
        self.P_out.value = self.P_in.value
        h_in = fluid_properties.enthalpy(self.model.design.fluid, T=self.T_in.value, P=self.P_in.value)
        self.T_out.value = self.T_amb.value + self.dT_approach.value
        self.h_out.value = fluid_properties.enthalpy(self.model.design.fluid, P=self.model.design.P_low, T=self.T_out.value)
        self.qdot.value = (self.h_out.value - h_in)*self.mdot_in.value
        return 

class Summary(Component):
    def __init__(self):
        self.W_turbine = Component.Input()
        self.W_compressor = Component.Input()
        self.Q_combustor = Component.Input()
        self.eta_cycle = Component.Output()
        return 
    
    def setup(self, **kwargs):
        return super().setup(**kwargs)

    def calculate(self):
        super().calculate()
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
