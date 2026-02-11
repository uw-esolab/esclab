# import numpy as __np
import matplotlib.pyplot as plt
import numpy as np
# import pandas as pd
import sys

# ------------------------------------------------------------------------
class Connection:
    def __init__(self, source, tol_rel, tol_abs, log_n_iter, learn_rate):

        self.tol_rel = tol_rel
        self.tol_abs = tol_abs 
        self.learn_rate = learn_rate
        
        BIG = 9e19
        self.err_rel = BIG
        self.err_abs = BIG

        self.source = source
        self.source_last_value = float('nan')

        self.log_n_iter = log_n_iter

        self.reset_for_step()
        return
    
    def reset_for_step(self):
        self.is_converged = False 
        self.n_iter = 0
        # [[value, abs err, rel err],]
        self.iter_log = np.zeros((max(self.log_n_iter,1),3))

    def compute_new_value(self):
        if not np.any(np.isnan(self.source_last_value)):
            return self.source_last_value + (self.source.value - self.source_last_value)*self.learn_rate
        else: 
            # nan case
            return self.source.value

    def check_convergence(self, new_value):
        # new_value = self.source.value 
        old_value = self.source_last_value
            
        self.err_abs = np.abs(new_value-old_value)
        self.err_rel = self.err_abs/np.maximum(np.abs(old_value),1.e-19)
        self.is_converged = np.all(self.err_abs < self.tol_abs) and np.all(self.err_rel < self.tol_rel)

        if self.log_n_iter > 0:
            la = np.array([self.source.value, self.err_abs, self.err_rel])
            if self.n_iter < self.log_n_iter: 
                self.iter_log[self.n_iter,:] = la
            else:
                self.iter_log = np.roll(self.iter_log, -1)
                self.iter_log[-1,:] = la

        self.n_iter += 1        
        return self.is_converged
    
# =========================================================================================   

class Component:
    class __io_base:
        def __init__(self):
            self.name = ''
            self.units = ''
            self.value = float('nan')
            self.is_connected = False
            pass
    class Input(__io_base):
        def __init__(self, initial_value=1.):
            super().__init__()
            self.connection = None  #instance of class Connection()
            self.value = initial_value

        def update_from_connection(self):
            if not self.connection == None:
                #check for nan
                if not np.any(np.isnan(self.connection.source.value)):
                    new_value = self.connection.compute_new_value()
                    converged = self.connection.check_convergence(new_value)
                    self.connection.source_last_value = self.connection.source.value
                    self.value = new_value
                    return converged
                else:
                    return False

    class Output(__io_base):
        def __init__(self):
            super().__init__()

    class Parameter:
        def __init__(self, value=float('nan'), units='', std_dev = None):
            self.value = value 
            self.units = units
            self.std_dev = std_dev
    # --------------------------------------------------------
    
    def __init__(self):
        """
        Define inputs, outputs, and parameters. 

        my_param_1 = Component.Parameter(value_1)

        my_input_1 = Component.Input(initial_value_1)
        my_input_2 = Component.Input(initial_value_2)
        
        my_output_1 = Component.Output()
        """
        pass 
    
    def __get_io_items(self, item_type, connected_only=False):
        io_list = []
        for item_name in dir(self):
            item = getattr(self,item_name)
            if isinstance(item, item_type):
                if connected_only:
                    if item.is_connected:
                        io_list.append(item)
                else:
                    io_list.append(item)
        return io_list

    def get_inputs(self, connected_only=False):
        return self.__get_io_items(Component.Input, connected_only)

    def get_outputs(self, connected_only=False):
        return self.__get_io_items(Component.Output, connected_only)
    
    def auto_assign_names(self):
        # Make sure all inputs/outputs have a name assigned
        allnames = []
        for member in dir(self):
            try:
                mo = getattr(self, member)
                if isinstance(mo, Component.Input) or isinstance(mo, Component.Output):
                    if self.name == '':
                        cname = type(self).__name__ 
                    else: 
                        cname = self.name
                    mo.name = cname + '.' + member
                    if isinstance(mo, Component.Output):
                        allnames.append(mo.name)
            except:
                    pass
        return allnames

    def setup(self, **kwargs):
        """
        Do initial calculations.
        Pre-simulation calculations here
        """
        pass
    
    def calculate(self):
        pass 
    
    def converge(self):
        pass 


# =========================================================================================        
class Model:
    class Settings:
        def __init__(self):
            self.timestep = 1  #sec
            self.start_time = 0  #sec
            self.stop_time = 24*3600 #sec
            self.tol_rel_global = 1.e-6
            self.tol_abs_global = 1.e-6
            self.max_iterations = 50 
    
    class OnlinePlotter:
        """
        y1 | [a, b, c, ...] component input or output values to be plotted on axis 1
        y2 | [d, e, f, ...] component input or output values to be plotted on axis 2
        """
        def __init__(self, y1, y2, y1lim, y2lim, y1label, y2label, nmax_points, update_every):
            assert isinstance(y1, type([]))
            assert isinstance(y2, type([])) or y2 == None

            self.current_step = -1

            self.nmax_points = nmax_points
            self.update_every = update_every
            self.y1_items = y1 
            self.y2_items = y2 
            colors = [
                '#377eb8', 
                '#ff7f00', 
                '#4daf4a',
                '#f781bf', 
                '#a65628', 
                '#984ea3',
                '#999999', 
                '#e41a1c', 
                '#dede00',
                "#00ffe5",
                "#5d0ab6",
                ]
            # Enable interactive mode
            plt.ion()

            # Create the figure and axis
            self.fig, self.ax1 = plt.subplots()
            self.ax2 = plt.twinx() if y2 != None else None

            self.y1_lines = []
            self.y2_lines = []

            c = 0
            for i in range(len(y1)):
                self.y1_lines.append( self.ax1.plot([], [], label=y1[i].name, color=colors[c%len(colors)])[0])  # Initialize an empty line
                c+=1
            if y2 != None:
                for i in range(len(y2)):
                    self.y2_lines.append( self.ax2.plot([], [], label=y2[i].name, color=colors[c%len(colors)])[0])  # Initialize an empty line
                    c+=1
            
            
            lines = self.y1_lines + self.y2_lines
            labels = [line.get_label() for line in lines]
            self.ax1.legend(lines, labels, loc='upper center', bbox_to_anchor=(0.5, 1.15), ncol=min(len(lines),3))

            self.ax1.set_xlabel('Time')

            # Set initial plot limits
            if y1lim != None:
                self.ax1.set_ylim(*y1lim)
            self.ax1.set_ylabel(y1label)
            if y2 != None:
                if y2lim != None:
                    self.ax2.set_ylim(*y2lim)
                self.ax2.set_ylabel(y2label)

            # Initialize data containers
            self.x_data =  np.zeros((self.nmax_points))
            self.y1_data = np.zeros((len(y1), self.nmax_points))
            if y2 != None:
                y2len = len(y2)
            else:
                y2len = 1
            self.y2_data = np.zeros((y2len, self.nmax_points))
            self.fig.tight_layout()
            return 
        
        def log_step(self, time):
            
            self.current_step += 1

            self.x_data = np.roll(self.x_data, -1)
            self.x_data[-1] = time

            y1vals = np.zeros(len(self.y1_items))
            if self.y2_items != None:
                y2vals = np.zeros(len(self.y2_items))
            # Append new data
            for j,yval in enumerate(self.y1_items):
                y1vals[j] = yval.value
            self.y1_data = np.roll(self.y1_data, -1, axis=1)
            self.y1_data[:,-1] = y1vals[:]
            
            if self.y2_items != None:
                for j,yval in enumerate(self.y2_items):
                    y2vals[j] = yval.value
                self.y2_data = np.roll(self.y2_data, -1, axis=1)
                self.y2_data[:,-1] = y2vals[:]

            if self.current_step % self.update_every == 0:
                self.__refresh_plot()

        def __refresh_plot(self):
            """
            x,y are numpy arrays
            """
            # update y1
            # Update the line data
            for i in range(len(self.y1_lines)):
                self.y1_lines[i].set_xdata(self.x_data)
                self.y1_lines[i].set_ydata(self.y1_data[i,:])
            for i in range(len(self.y2_lines)):
                self.y2_lines[i].set_xdata(self.x_data)
                self.y2_lines[i].set_ydata(self.y2_data[i,:])

            # Adjust the view
            self.ax1.relim()
            self.ax1.autoscale_view()
            if self.y2_items != None:
                self.ax2.relim()
                self.ax2.autoscale_view()

            # Redraw the plot
            self.fig.canvas.draw()
            self.fig.canvas.flush_events()
        
            # # Turn off interactive mode and show the final plot
            # plt.ioff()
            # plt.show()
    # ------- end OnlinePlotter ----------------------------------------
    def __init__(self):
        self.__components = []
        self.settings = Model.Settings()
        self.is_initialized = False
        self.time = 0
        self.plotters = []
        return
    
    def add_plotter(self, y1, y2=None, y1lim=None, y2lim=None, y1label='', y2label='', nmax_points = 1000, update_every=1):
        if not isinstance(y1, type([])):
            y1t = [y1]
        else:
            y1t = y1
        y2t = y2
        if y2 != None:
            if not isinstance(y2, type([])):
                y2t = [y2]

        self.plotters.append(Model.OnlinePlotter(y1t, y2t, y1lim, y2lim, y1label, y2label, nmax_points, update_every))

    def connect(self, source, destination, tol_rel = 1.e-6, tol_abs=1.e-6, log_n_iter = 0, learn_rate = 1.):
        
        if not isinstance(source, Component.Output):
            raise RuntimeError(f"Source connection object must be of type 'Component.Output'")
        if not isinstance(destination, Component.Input):
            raise RuntimeError(f"Destination connection object must be of type 'Component.Input'")

        destination.connection = Connection(source, tol_rel, tol_abs, log_n_iter, learn_rate)
        destination.is_connected = True
        source.is_connected = True
        return

    def initialize(self):
        self.time = self.settings.start_time
        self.iteration = -1

        output_names = ['time','iterations']
        for item in dir(self):
            itemobj = getattr(self, item)
            if isinstance(itemobj, Component):
                # Handle component setup here
                itemobj.model = self
                itemobj.name = ''
                itemobj.setup()
                cnames = itemobj.auto_assign_names()
                self.__components.append(itemobj)

                # Record all output data 
                output_names += cnames
        # Construct the historian database
        nstep = int((self.settings.stop_time - self.settings.start_time)/self.settings.timestep)
        self.historian = dict([[n,np.ones(nstep)*float('nan')] for n in output_names])
        
        self.is_initialized = True
        return 
    
    def step(self, n_steps = 1):
        # Check overall initialization for the model
        if not self.is_initialized:
            self.initialize()
        
        self.iteration = -1  # reset the current iteration
        err_rel_history = []
        err_abs_history = []
        # main loop 
        for i in range(self.settings.max_iterations):
            all_converged = True 
            self.iteration += 1

            # Run through list of components, gathering and updating inputs 
            max_abs_err = 0.
            max_rel_err = 0.
            for component in self.__components:
                # Update connections first
                for input in component.get_inputs(connected_only=True):
                    all_converged = all_converged & input.update_from_connection()
                    max_abs_err = max(max_abs_err, np.max(input.connection.err_abs))
                    max_rel_err = max(max_rel_err, np.max(input.connection.err_rel))

                # Calculate
                component.calculate()
            
            if all_converged:
                # print(f'Iterations: {i}')
                break
            if i == self.settings.max_iterations-1:
                for component in self.__components:
                    for input in component.get_inputs(connected_only=True):
                        if not input.connection.is_converged:
                            ps = f'{self.time} | {input.connection.source.name} not converged after {input.connection.n_iter} iterations.'
                            if input.connection.log_n_iter > 0:
                                ps += ' Iter log: ' + ' '.join(input.connection.iter_log[:,0].astype(str))
                            print(ps)

            err_rel_history.append(max_rel_err)
            err_abs_history.append(max_abs_err)

        # Convergence
        current_step = int( (self.time - self.settings.start_time)/self.settings.timestep )
        for component in self.__components:
            component.converge()
            for input in component.get_inputs(connected_only=True):
                input.connection.reset_for_step()

            self.historian['time'][current_step] = self.time
            self.historian['iterations'][current_step] = self.iteration
            for output in component.get_outputs():
                try:
                    if not isinstance(output.value, np.ndarray):
                        self.historian[output.name][current_step] = output.value
                except ValueError:
                    pass
        
        # Plotters 
        for plotter in self.plotters:
            plotter.log_step(self.time)

        self.time += self.settings.timestep

        # Terminal update
        percent = (self.time / (self.settings.stop_time - self.settings.start_time)) * 100
        if int(percent*1000) % 10 == 0:
            bar = '█' * int(percent // 2) + '-' * (50 - int(percent // 2))
            sys.stdout.write(f'\r|{bar}| {percent:.1f}% ({self.time:.2f} sec)')
            sys.stdout.flush()

        return 

