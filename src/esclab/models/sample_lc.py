from esclab.simulate import *
import numpy as np

class FF(Component):
    def __init__(self):
        self.signal = Component.Output()
        return super().__init__()

    def presim_setup(self):
        pass

    def calculate(self):
        self.signal.v = 1 if self.model.time % 24 < 8 else 0
        return 
    
class Object(Component):
    signal = Component.Input(0.)
    gain = Component.Parameter(100)
    tau = Component.Parameter(12)
    Tamb = Component.Parameter(200)
    T = Component.Output()

    def presim_setup(self):
        self.T0 = 0  # initial temp

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

model.wait_for_plots()