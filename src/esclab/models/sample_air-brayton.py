from esclab.simulate import *
from eeslib.functions import convert, converttemp
from esclab.components.brayton_simple import *

class AirBraytonModel(Model):
    class DesignValues:
        def __init__(self):
            # self.PR = np.array([3.,4.,.5])
            self.PR = 4.
            self.P_low = 1*convert('bar','Pa')
            self.T_amb = converttemp('F','K', 90)
            self.T_turb_in = converttemp('C','K',1200)
            self.mdot = 100  #kg/s
            self.fluid = 'Air'
    
    def __init__(self):
        super().__init__()
        self.design = AirBraytonModel.DesignValues()

# ==============================================================

# Create model
model = AirBraytonModel()

# Create component instances
model.weather = Weather()
model.turbine = Turbine()
model.compressor = Compressor()
model.combustor = Combustor()
model.cooler = Cooler()
model.summary = Summary()

model.initialize()

# Set parameters and initial values 
# turbine
model.turbine.eta_s   = 0.9
model.turbine.mdot_in = model.design.mdot
model.turbine.T_in    = model.design.T_turb_in
model.turbine.P_in    = model.design.P_low*model.design.PR
model.turbine.fluid = 'Air'
# Compressor
model.compressor.eta_s   = 0.85
model.compressor.PR      = 4.
model.compressor.mdot_in = model.design.mdot
model.compressor.T_in    = model.design.T_amb+20
model.compressor.P_in    = model.design.P_low
model.compressor.fluid = model.turbine.fluid
# Cooler 
model.cooler.mdot_in    = model.design.mdot
model.cooler.T_in       = model.design.T_amb+100
model.cooler.P_in       = model.design.P_low
model.cooler.T_amb      = model.design.T_amb
# Combustor
# initial value
model.combustor.P_in = model.design.P_low*model.design.PR
model.combustor.T_in = model.design.T_turb_in - 300
model.combustor.mdot_in  = model.design.mdot

# ---------------------------------------------
# Connections
# ---------------------------------------------
model.connect(model.turbine.mdot_out,     model.cooler.mdot_in       , )
model.connect(model.turbine.P_out,        model.cooler.P_in          , )
model.connect(model.turbine.T_out,        model.cooler.T_in          , )
model.connect(model.weather.T_amb,        model.cooler.T_amb         , )
model.connect(model.cooler.mdot_out,      model.compressor.mdot_in   , )
model.connect(model.cooler.P_out,         model.compressor.P_in      , )
model.connect(model.cooler.T_out,         model.compressor.T_in      , )
model.connect(model.compressor.mdot_out,  model.combustor.mdot_in    , )
model.connect(model.compressor.P_out,     model.combustor.P_in       , )
model.connect(model.compressor.T_out,     model.combustor.T_in       , )
model.connect(model.combustor.mdot_out,   model.turbine.mdot_in      , )
model.connect(model.combustor.P_out,      model.turbine.P_in         , )
model.connect(model.combustor.T_out,      model.turbine.T_in         , )
model.connect(model.turbine.W,            model.summary.W_turbine    , )
model.connect(model.compressor.W,         model.summary.W_compressor , )
model.connect(model.combustor.qdot,       model.summary.Q_combustor  , )

model.add_plotter([model.compressor.T_out, model.weather.T_amb, model.turbine.T_in, model.cooler.T_in],[model.compressor.W], y1label='Temperature', y2label='Work', update_every=10, nmax_points=300)
model.add_plotter([model.summary.eta_cycle], [model.summary.Q_combustor], y1label='Efficiency', y2label='Heat', update_every=10, nmax_points=100)
# model.add_plotter([model.summary.eta_cycle], [model.summary.Q_combustor], y1label='Efficiency', y2label='Heat', update_every=100, nmax_points=100)

# Optional: add a topology tab and export images.
# model.add_network_graph(show_tab=True, save_png=True, path_base='air_brayton_topology')
model.add_network_graph(show_tab=True)




model.settings.stop_time = 8760*3600*1/12
model.settings.start_time = 0
model.time = model.settings.start_time
model.settings.timestep = 3600

while model.time < model.settings.stop_time:
    model.step()

model.wait_for_plots()