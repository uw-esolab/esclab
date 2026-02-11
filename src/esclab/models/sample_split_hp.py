from esclab.simulate import *
import numpy as np
from eeslib.functions import convert
import time 


class SplitHP(Model):
    def __init__(self):
        super().__init__()

model = SplitHP()

# ------------------- low cycle ----------------
# Compressor
from esclab.components.brayton_simple import Compressor
# cascade HX
class CascadeHX(Component):
    def setup(self, **kwargs):
        self.approach     = Component.Parameter(5)

        self.m_dot_hi_in  = Component.Input()
        self.h_hi_in      = Component.Input()
        self.m_dot_low_in = Component.Input()
        self.h_low_in     = Component.Input()

        return super().setup(**kwargs)
# Valve

# Evaporator