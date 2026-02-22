from esclab.simulate import *
import numpy as np

class FF(Component):
    def __init__(self):
        self.signal = Component.Output()

    def presim_setup(self):
        pass

    def calculate(self):
        self.signal.v = 1 if self.model.time % 24 < 8 else 0
        return 
    
class Object(Component):
    def __init__(self):
        self.signal = Component.Input(0.)
        self.gain = Component.Parameter(100)
        self.tau = Component.Parameter(12)
        self.Tamb = Component.Parameter(200)
        self.T = Component.Output()
        self.T0 = 0  #iniital temp

    def calculate(self):
        # post-convergence calculations
        if self.model.is_converged:
            self.T0 = self.T.v
            return 
        
        # iteration calculations
        dTdt = -(self.T0 -(self.Tamb.v + self.signal.v*self.gain.v))/self.tau.v
        self.T.v = self.T0 + dTdt*self.model.settings.timestep
        return 
    
# ---------------------------------------------------
model = Model()
model.settings.timestep = .1
model.settings.stop_time = 24*7

model.ff = FF()
model.object = Object()

model.initialize()

model.connect(model.ff.signal,    model.object.signal)

model.add_plotter([model.object.T],[model.ff.signal], nmax_points=1000, update_every=25)


while model.time < model.settings.stop_time:
    model.step()
