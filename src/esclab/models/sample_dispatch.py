from esclab.simulate import *
import numpy as np
import gurobipy as gp 
from gurobipy import GRB 
from eeslib.functions import convert
import time 

class Design:
    def __init__(self):
        self.producer_rating = 1000  # J/s
        self.storage_capacity = self.producer_rating * 5 * convert('hr','s')  # J
        self.storage_initial_charge = self.storage_capacity * 0.5  # J
        self.c_rate = 0.25  # 1/hr  Fraction of total storage capacity that can be charged or discharged per hour
        self.consumer_efficiency = 0.5  
        self.optimization_horizon = 48 * convert('hr','s')  # s 
        self.control_horizon = 24 * convert('hr','s')  # s 


class Storage(Component):
    def __init__(self):
        super().__init__()
        self.capacity           = Component.Parameter(1)
        self.capacity_init      = Component.Parameter(1)
        self.c_rate             = Component.Parameter(1)  #Max charge/discharge rate, as a fraction of present charge
        self.loss_rate          = Component.Parameter(1)  #Loss as a fraction of capacity per s

        self.flow_in            = Component.Input()  #J/s
        self.flow_out           = Component.Input()  #J/s

        self.flow_in_actual     = Component.Output()  #J/s
        self.flow_out_actual    = Component.Output()  #J/s
        self.charge             = Component.Output()  #J
        self.losses             = Component.Output()  #J/s
        self.last_charge = 0.  #state of charge

    def setup(self, **kwargs):
        self.last_charge = self.capacity_init.value
        pass 

    def calculate(self):
        # Reality checks
        assert self.flow_in.value >= 0
        assert self.flow_out.value >= 0

        # Losses from last time step
        losses = max(self.last_charge * self.loss_rate.value * self.model.settings.timestep , 0)
        self.losses.value = losses
        # Energy balance
        self.charge.value = self.last_charge - losses + (self.flow_in.value - self.flow_out.value)*self.model.settings.timestep 
        
        self.flow_out_actual.value = self.flow_out.value
        self.flow_in_actual.value = self.flow_in.value 
        # Correct flows if resulting charge exceeds min/max capacity
        if self.charge.value < 0.:
            self.flow_out_actual.value += self.charge.value/self.model.settings.timestep 
            self.charge.value = 0.
        if self.charge.value >= self.capacity.value:
            self.flow_in_actual.value -= (self.charge.value - self.capacity.value)/self.model.settings.timestep
            self.charge.value = self.capacity.value
        if self.flow_out.value > self.flow_out_actual.value:
            pass
        return

    def converge(self):
        self.last_charge = self.charge.value

class Producer(Component):
    def __init__(self):
        super().__init__()
        self.rated_power    = Component.Parameter(1.)
        self.horizon        = Component.Parameter(24)

        self.power          = Component.Output()
        self.power_forecast = Component.Output()

    def setup(self, **kwargs):
        # initialize the power forecast array
        fctimes = np.arange(self.model.settings.start_time, self.model.settings.start_time + self.horizon.value, self.model.settings.timestep)
        self.power_forecast.value = self.__generate_forecast(fctimes)

    def calculate(self):

        self.power.value = self.power_forecast.value[0] 
        # Only update the new forecast value on the first iteration
        if self.model.iteration == 0:
            self.new_fc_value = self.__generate_forecast(np.array([self.model.time + self.horizon.value]))

        return
    
    def converge(self):
        self.power_forecast.value = np.roll(self.power_forecast.value, -1)
        self.power_forecast.value[-1] = self.new_fc_value[0]

    def __generate_forecast(self, time_values):
        tm_adj = (time_values/self.model.settings.timestep % 24.)*(np.pi)/(24.)
        y = np.sin(tm_adj)*self.rated_power.value #* 3/2 - self.rated_power.value*1/2
        y[y<0] = 0.
        lmask = np.abs(np.random.uniform(0,2,len(time_values))) - 0.5
        lmask[lmask<0] = 0.
        lmask[lmask>1] = 1.
        y = lmask*y
        return y
        
class Consumer(Component):
    def __init__(self):
        super().__init__()
        self.horizon    = Component.Parameter(1)
        self.efficiency = Component.Parameter(1)
        
        self.flow_in    = Component.Input()

        self.price      = Component.Output()
        self.revenue    = Component.Output()
        self.price_forecast = Component.Output()


    def setup(self, **kwargs):
        
        fctimes = np.arange(self.model.settings.start_time, self.model.settings.start_time+self.horizon.value, self.model.settings.timestep)
        self.price_forecast.value = self.__generate_forecast(fctimes)

    def calculate(self):
        self.price.value = self.price_forecast.value[0]
        # Only update the new forecast value on the first iteration
        if self.model.iteration == 0:
            self.new_fc_value = self.__generate_forecast(np.array([self.model.time + self.horizon.value]))

        self.revenue.value = self.flow_in.value * self.model.settings.timestep * self.efficiency.value * self.price.value

    def converge(self):
        self.price_forecast.value = np.roll(self.price_forecast.value, -1)
        self.price_forecast.value[-1] = self.new_fc_value[0]

    def __generate_forecast(self, time_values):
        tm_adj = (time_values % 24.)*(2*np.pi)/(24.)
        y = (np.cos(tm_adj)*np.random.normal(0,.1,len(tm_adj))+1)*1
        return y
    
class Scheduler(Component):
    def __init__(self):
        super().__init__()

        self.optimization_horizon    = Component.Parameter(1.)
        self.control_horizon         = Component.Parameter(1.)
        self.storage_initial_charge  = Component.Parameter(1.)
        self.storage_capacity        = Component.Parameter(1.)
        self.c_rate                  = Component.Parameter(1.)
        self.consumer_efficiency     = Component.Parameter(1.)

        # initialize schedules with dummy arrays of the right length
        self.charge_schedule    = Component.Input(np.array([]))
        self.price_schedule     = Component.Input(np.array([]))
        self.storage_charge     = Component.Input()

        self.flow_to_consumer   = Component.Output()
        self.flow_from_producer = Component.Output()
        self.price_now          = Component.Output()
        self.charge_avail_now   = Component.Output()
        self.num_iter           = Component.Output()

    def setup(self, **kwargs):
        # initialize schedules with dummy arrays of the right length
        da = np.ones(int(self.optimization_horizon.value/self.model.settings.timestep))
        self.charge_schedule.value = da
        self.price_schedule.value = da

        # Initialize storage charge state
        self.last_charge_state = self.storage_initial_charge.value

    def calculate(self):
        t_rel = int((self.model.time % self.control_horizon.value)*convert('s','hr'))
        if t_rel == 0:
            # Check for missing input data
            if np.isnan(self.storage_charge.value):
                # provide temporary data
                self.active_schedule = [{'flow_from_producer':0, 'flow_to_consumer':0}]
            else:
                if self.model.iteration < 2:
                    self.active_schedule = self.__run_opt_model()
        self.price_now.value = self.price_schedule.value[t_rel]
        self.charge_avail_now.value = self.active_schedule[t_rel]['charge']
        self.flow_to_consumer.value = self.active_schedule[t_rel]['flow_to_consumer']
        self.flow_from_producer.value = self.active_schedule[t_rel]['flow_from_producer']

    def converge(self):
        # Only update the initial charge state for the next step after the current step has converged
        self.last_charge_state = self.storage_charge.value
        self.num_iter.value = self.model.iteration

    def __run_opt_model(self):
        # md = self.model.design

        nt = len(self.price_schedule.value)
        T = range(nt)
        qmax = self.storage_capacity.value * self.c_rate.value / convert('hr','s')  # W
        smax = self.storage_capacity.value
        S0 = self.last_charge_state
        
        # Optimization Model
        om = gp.Model('storage_optimization')
        om.setParam('OutputFlag',0)
        
        # Variables
        q_out = om.addVars(T, lb=0., ub=qmax, name='q_out')  # flow to consumer
        q_in =  om.addVars(T, lb=0., ub=qmax, name='q_in')   # flow in used
        s = om.addVars(T, lb=0., ub=smax, name='s')  # storage inventory

        # Objective function
        eta = self.consumer_efficiency.value
        om.setObjective(gp.quicksum(q_out[t]*eta*self.price_schedule.value[t]*model.settings.timestep for t in T), GRB.MAXIMIZE)
        
        # ------------ Constraints
        # flow in is no greater than available flow in 
        om.addConstrs((q_in[t] <= self.charge_schedule.value[t] for t in T), 'flow_utilized')

        # energy balance on storage based on power consumed
        om.addConstrs((s[t] == (s[t-1] if t>0 else S0) - (q_out[t] - q_in[t])*model.settings.timestep for t in T), 'store_balance')

        om.optimize()

        if om.status == GRB.OPTIMAL:
            res = []
            for t in T:
                res.append({
                    'flow_from_producer':q_in[t].X,
                    'flow_to_consumer':q_out[t].X,
                    'charge':s[t].X,
                })
            return res 
            # print("Optimal objective value:", model.objVal)
        else:
            raise RuntimeError("Optimization model did not converge")


if __name__ == "__main__":
    model = Model()
    model.settings.stop_time = 30*convert('day','s')
    model.settings.timestep = 3600

    model.producer = Producer()
    model.storage = Storage()
    model.consumer = Consumer()
    model.scheduler = Scheduler()

    # ---------------------- Set parameters and initial values

    design = Design()

    # storage
    model.storage.capacity.value       = design.storage_capacity
    model.storage.c_rate.value         = design.c_rate  #Max charge/discharge rate, as a fraction of present charge
    model.storage.capacity_init.value  = design.storage_initial_charge
    model.storage.loss_rate.value      = 0.001 / convert('h','s')  #Loss as a fraction of capacity per s
    model.storage.flow_in.value        = 0  #J/s
    model.storage.flow_out.value       = 0  #J/s

    # Producer
    model.producer.rated_power.value    = design.producer_rating
    model.producer.horizon.value        = design.optimization_horizon

    # Consumer
    model.consumer.horizon.value    = design.optimization_horizon
    model.consumer.efficiency.value = design.consumer_efficiency
    model.consumer.flow_in.value    = 0.

    # Scheduler
    model.scheduler.optimization_horizon.value    = design.optimization_horizon
    model.scheduler.control_horizon.value         = design.control_horizon
    model.scheduler.storage_initial_charge.value  = design.storage_initial_charge
    model.scheduler.storage_capacity.value        = design.storage_capacity
    model.scheduler.c_rate.value                  = design.c_rate
    model.scheduler.consumer_efficiency.value     = design.consumer_efficiency


    model.initialize()
    
    # ------------- from ----------------------------- to ----------------------
    model.connect(model.scheduler.flow_from_producer    , model.storage.flow_in           ) 
    model.connect(model.scheduler.flow_to_consumer      , model.storage.flow_out          )
    model.connect(model.storage.flow_out_actual         , model.consumer.flow_in          ) #, log_n_iter=model.settings.max_iterations)
    model.connect(model.producer.power_forecast         , model.scheduler.charge_schedule )
    model.connect(model.consumer.price_forecast         , model.scheduler.price_schedule  )
    model.connect(model.storage.charge                  , model.scheduler.storage_charge  )

    model.add_plotter([model.storage.flow_in, model.storage.flow_out], [model.storage.charge, model.scheduler.charge_avail_now], update_every=24, nmax_points=24*7)
    model.add_plotter([model.consumer.revenue], [model.consumer.price], update_every=24, nmax_points=24*7)
    model.add_plotter([model.scheduler.num_iter], update_every=24)

    start = time.time()
    while model.time < model.settings.stop_time:
        model.step()
    print(f'time: {time.time()-start}')

    plt.show(block=True)