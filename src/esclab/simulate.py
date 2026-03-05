import pyqtgraph as qtg
from pyqtgraph.Qt import QtCore, QtWidgets, QtGui
import numpy as np
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
            return self.source_last_value + (self.source.v - self.source_last_value)*self.learn_rate
        else: 
            # nan case
            return self.source.v

    def check_convergence(self, new_value):
        # new_value = self.source.v 
        old_value = self.source_last_value
            
        self.err_abs = np.abs(new_value-old_value)
        self.err_rel = self.err_abs/np.maximum(np.abs(old_value),1.e-19)
        self.is_converged = np.all(self.err_abs < self.tol_abs) and np.all(self.err_rel < self.tol_rel)

        if self.log_n_iter > 0:
            la = np.array([self.source.v, self.err_abs, self.err_rel])
            if self.n_iter < self.log_n_iter: 
                self.iter_log[self.n_iter,:] = la
            else:
                self.iter_log = np.roll(self.iter_log, -1)
                self.iter_log[-1,:] = la

        self.n_iter += 1        
        return self.is_converged
    
# =========================================================================================   

class Component:
    # ------------------------------------------------------------------------
    class __io_base:
        """
        Base for Input, Parameter, and Output classes.

        Contains the value (v), units, name, and connection status.
        """
        def __init__(self):
            self.name = ''
            self.units = ''
            self.v = float('nan')
            self.is_connected = False
            return
        def set(self, value=None, units=None, name=None):
            if value is not None: 
                self.v = value
            if units is not None:
                self.units = units
            if name is not None:
                self.name = name
    # ------------ end class __io_base ----------------------------------------
    # -------------------------------------------------------------------------
    class Input(__io_base):
        def __init__(self, initial_value=1.):
            super().__init__()
            self.connection = None  #instance of class Connection()
            self.v = initial_value

        def update_from_connection(self):
            if not self.connection == None:
                #check for nan
                if not np.any(np.isnan(self.connection.source.v)):
                    new_value = self.connection.compute_new_value()
                    converged = self.connection.check_convergence(new_value)
                    self.connection.source_last_value = self.connection.source.v
                    self.v = new_value
                    return converged
                else:
                    return False
    # ------------ end class Input -------------------------------------------
    # ------------------------------------------------------------------------
    class Output(__io_base):
        def __init__(self):
            super().__init__()
    # ----------- end class Output ------------------------------------------
    # ------------------------------------------------------------------------
    class Parameter:
        def __init__(self, value=float('nan'), units='', std_dev = None):
            self.v = value 
            self.units = units
            self.std_dev = std_dev
    # ----------- end class Parameter -----------------------------------------
    
    # ------------------------------------------------------------------------
    #  Component class methods
    # ------------------------------------------------------------------------
    def __init__(self):
        """
        Define inputs, outputs, and parameters. 
        """
        return 
    
    def __get_io_items(self, item_type, connected_only=False):
        """
        Collect all inputs or outputs from this component, depending on the item_type argument.

        Parameters
        ----------
        item_type : type
            Component.Input or Component.Output
        connected_only : bool
            If True, only return inputs/outputs that are connected to another component. 
            If False, return all inputs/outputs regardless of connection status.
        """
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
        """
        Helper function to collect all inputs from this component.

        Parameters
        ----------
        connected_only : bool
            If True, only return inputs that are connected to another component. 
            If False, return all inputs regardless of connection status.
        """
        return self.__get_io_items(Component.Input, connected_only)

    def get_outputs(self, connected_only=False):
        """
        Helper function to collect all outputs from this component.
        
        Parameters
        ----------
        connected_only : bool
            If True, only return outputs that are connected to another component. 
            If False, return all outputs regardless of connection status.
        """
        return self.__get_io_items(Component.Output, connected_only)
    
    def set_values_from_dict(self, values_dict):
        """
        Helper function to set values of inputs and parameters from a dictionary. 
        The keys of the dictionary should match the attribute names of the inputs and parameters.

        Parameters
        ----------
        values_dict : dict
            Dictionary containing values to set for inputs and parameters. Keys should match attribute names.
        """
        for key, value in values_dict.items():
            if hasattr(self, key):
                attr = getattr(self, key)
                if isinstance(attr, Component.Input) or isinstance(attr, Component.Parameter):
                    attr.set(value=value)
                else:
                    raise RuntimeError(f"Attribute '{key}' is not an input or parameter.")
            else:
                raise RuntimeError(f"Component does not have an attribute named '{key}'.")
    
    def auto_assign_names(self):
        """
        Automatically assign names to all inputs and outputs of this component based on the 
        component's name and the attribute name of the input/output.
        Make sure all inputs/outputs have a name assigned
        """
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

    def presim_setup(self, **kwargs):
        """
        Do initial calculations.
        Pre-simulation calculations here
        """
        pass
    
    def calculate(self):
        """
        Do calculations for this component. This is where the main logic of the component goes, 
        and where outputs are calculated from inputs and parameters.
        
        Within the calculate step, the following flags are available to control the flow of calculations:
        ---------------------------------------------------------------------------------------------------------
        is_first_step       | True only on the first step() call of the simulation. Use for any calculations that 
                            | should only be done once at the beginning of the simulation.
        ---------------------------------------------------------------------------------------------------------
        is_first_iteration  | True only on the first iteration of each step() call. Use for any calculations that 
                            | should only be done once per time step, such as calculations that should be done 
                            | in order of components.
        ---------------------------------------------------------------------------------------------------------
        is_converged        | True only on the final iteration of each step() call, after convergence is reached. 
                            | Use for any calculations that should only be done once per time step, and that 
                            | require convergence to be reached first.
        ---------------------------------------------------------------------------------------------------------
        """
        pass 
    
    def postsim_calcs(self):
        pass


# =========================================================================================        
class Model:
    """
    Model class to hold components and run the simulation.

    The model workflow is as follows:
    * Use a script to create an instance of the Model class, and add components as attributes of the model.
    * Use the connect() method to connect component outputs to other component inputs.
    * Create any plotters using the add_plotter() method, which takes lists of component inputs and outputs 
      to be plotted on the y1 and y2 axes, respectively.
    * Optionally, call initialize() to initialize the model before running the simulation. This will automatically call the presim_setup() method of each component, and assign names to all inputs and outputs.
    * In a loop for all time steps, call the step() method to run the simulation. This will run through all components and update inputs from connections, calculate outputs, check for convergence, and log data for plotting and historian.

    Substeps in the time series workflow are as follows:

    Function        | Flags active       | Description
    ----------------|--------------------|------------------------------------------------------
    presim_setup()  |                    | Calculations before simulation begins. Data from other components is not available.
    step()          | is_first_step      | Calculations done only on the first timestep call 
    step()          | is_first_iteration | Calculations done only on the first iteration of each timestep. Components called in order. 
    step()          |                    | Subsequent iteration calculations requiring convergence.
    step()          | is_converged       | Final pass of calculations done after convergence is reached.
    postsim_calcs() |                    | Calculations done after the final time step is reached. Data from all time steps is available.
    """

    # ------------------------------------------------------------------------
    class Settings:
        timestep = 1  #sec
        start_time = 0  #sec
        stop_time = 24*3600 #sec
        tol_rel_global = 1.e-6
        tol_abs_global = 1.e-6
        max_iterations = 50 
        
        def __init__(self):
            pass
    # End class Settings -----------------------------------------------------
    
    # ------------------------------------------------------------------------
    class OnlinePlotter:
        app = None
        main_window = None
        tab_widget = None
        instances = []
        n_plotters = 0
        font_size_pt = 10
        min_font_size_pt = 6
        max_font_size_pt = 30

        """
        y1 | [a, b, c, ...] component input or output values to be plotted on axis 1
        y2 | [d, e, f, ...] component input or output values to be plotted on axis 2
        """
        def __init__(self, y1, y2, y1lim, y2lim, y1label, y2label, nmax_points, update_every, plotter_size=(.9,.9), tab_title=None):
            assert isinstance(y1, type([]))
            assert isinstance(y2, type([])) or y2 == None

            self.current_step = -1
            self.y1label = y1label
            self.y2label = y2label

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
            
            # Create shared pyqtgraph application and tabbed window
            self.app = Model.OnlinePlotter.app
            if self.app is None:
                self.app = qtg.mkQApp()
                Model.OnlinePlotter.app = self.app

            if Model.OnlinePlotter.main_window is None:
                main_window = QtWidgets.QMainWindow()
                main_window.setWindowTitle("Simulation Plot")
                tab_widget = QtWidgets.QTabWidget()

                controls_widget = QtWidgets.QWidget()
                controls_layout = QtWidgets.QHBoxLayout(controls_widget)
                controls_layout.setContentsMargins(6, 6, 6, 0)
                controls_layout.setSpacing(6)

                font_up_button = QtWidgets.QPushButton("🗚")
                font_down_button = QtWidgets.QPushButton("🗛")
                font_up_button.setFixedSize(30, 30)
                font_down_button.setFixedSize(30, 30)
                controls_layout.addWidget(font_up_button)
                controls_layout.addWidget(font_down_button)
                controls_layout.addStretch()

                container = QtWidgets.QWidget()
                container_layout = QtWidgets.QVBoxLayout(container)
                container_layout.setContentsMargins(0, 0, 0, 0)
                container_layout.setSpacing(0)
                container_layout.addWidget(controls_widget)
                container_layout.addWidget(tab_widget)
                main_window.setCentralWidget(container)

                font_down_button.clicked.connect(lambda: Model.OnlinePlotter.adjust_font_size(-1))
                font_up_button.clicked.connect(lambda: Model.OnlinePlotter.adjust_font_size(1))

                # Handle plotter size. If fractions are given, resize based on screen size. If absolute values are
                # given, use those. If None, use default size.
                screen = QtWidgets.QApplication.primaryScreen()
                if screen is not None:
                    screen_rect = screen.availableGeometry()
                    if plotter_size is not None:
                        if plotter_size[0] <= 1. and plotter_size[1] <= 1.:
                            main_window.resize(int(screen_rect.width() * plotter_size[0]), int(screen_rect.height() * plotter_size[1]))
                        else:
                            main_window.resize(plotter_size[0], plotter_size[1])
                    frame = main_window.frameGeometry()
                    frame.moveCenter(screen_rect.center())
                    main_window.move(frame.topLeft())
                elif plotter_size is not None:
                    if plotter_size[0] <= 1. and plotter_size[1] <= 1.:
                        main_window.resize(1000, 600)
                    else:
                        main_window.resize(plotter_size[0], plotter_size[1])

                main_window.show()
                Model.OnlinePlotter.main_window = main_window
                Model.OnlinePlotter.tab_widget = tab_widget

            self.win = qtg.GraphicsLayoutWidget()
            Model.OnlinePlotter.n_plotters += 1
            tab_label = tab_title if tab_title not in [None, ''] else y1label if y1label not in [None, ''] else f"Plot {Model.OnlinePlotter.n_plotters}"
            Model.OnlinePlotter.tab_widget.addTab(self.win, tab_label)

            # Create primary plot
            self.ax1 = self.win.addPlot()
            self.ax1.setLabel('bottom', 'Time')
            self.ax1.setLabel('left', y1label)
            self.legend_y1 = self.ax1.addLegend(offset=(10, 10))
            self.legend_y1.setBrush(QtGui.QBrush(QtGui.QColor(255, 255, 255, 40)))
            self.legend_y2 = None
            self.ax1.showGrid(x=True, y=True, alpha=0.3)
            
            if y1lim != None:
                self.ax1.setYRange(*y1lim)

            # Create lines for y1 axis
            self.y1_lines = []
            c = 0
            for i in range(len(y1)):
                pen = qtg.mkPen(color=colors[c%len(colors)], width=2)
                line = self.ax1.plot([], [], pen=pen, name=y1[i].name)
                self.y1_lines.append(line)
                c += 1
            
            # Create secondary y-axis if needed
            self.y2_lines = []
            if y2 != None:
                self.ax2 = qtg.ViewBox()
                self.ax1.showAxis('right')
                self.ax1.scene().addItem(self.ax2)
                self.ax1.getAxis('right').linkToView(self.ax2)
                self.ax2.setXLink(self.ax1)
                self.ax1.getAxis('right').setLabel(y2label, color='#c0c0c0')
                self.legend_y2 = qtg.LegendItem(offset=(-10, 10))
                self.legend_y2.setBrush(QtGui.QBrush(QtGui.QColor(255, 255, 255, 40)))
                self.legend_y2.setParentItem(self.ax1.vb)
                
                if y2lim != None:
                    self.ax2.setYRange(*y2lim)
                
                # Update views when resized
                def updateViews():
                    self.ax2.setGeometry(self.ax1.vb.sceneBoundingRect())
                    self.ax2.linkedViewChanged(self.ax1.vb, self.ax2.XAxis)
                
                updateViews()
                self.ax1.vb.sigResized.connect(updateViews)
                
                # Create lines for y2 axis
                for i in range(len(y2)):
                    pen = qtg.mkPen(color=colors[c%len(colors)], width=2, style=QtCore.Qt.DashLine)
                    line = qtg.PlotDataItem([], [], pen=pen, name=y2[i].name)
                    self.ax2.addItem(line)
                    self.legend_y2.addItem(line, y2[i].name)
                    self.y2_lines.append(line)
                    c += 1

            # Initialize data containers
            self.x_data =  np.zeros((self.nmax_points))
            self.y1_data = np.zeros((len(y1), self.nmax_points))
            if y2 != None:
                y2len = len(y2)
            else:
                y2len = 1
            self.y2_data = np.zeros((y2len, self.nmax_points))

            Model.OnlinePlotter.instances.append(self)
            self.apply_font_size()
            return 

        @classmethod
        def adjust_font_size(cls, delta):
            new_font_size = max(cls.min_font_size_pt, min(cls.max_font_size_pt, cls.font_size_pt + delta))
            if new_font_size == cls.font_size_pt:
                return

            cls.font_size_pt = new_font_size
            for plotter in cls.instances:
                plotter.apply_font_size()

        def apply_font_size(self):
            tick_font = QtGui.QFont()
            tick_font.setPointSize(Model.OnlinePlotter.font_size_pt)
            label_style = {'font-size': f"{Model.OnlinePlotter.font_size_pt + 2}pt"}

            left_axis = self.ax1.getAxis('left')
            bottom_axis = self.ax1.getAxis('bottom')
            left_axis.setStyle(tickFont=tick_font)
            bottom_axis.setStyle(tickFont=tick_font)

            self.ax1.setLabel('left', self.y1label, **label_style)
            self.ax1.setLabel('bottom', 'Time', **label_style)

            if self.legend_y1 is not None:
                label_font = self.legend_y1.font()
                label_font.setPointSize(Model.OnlinePlotter.font_size_pt)
                self.legend_y1.setFont(label_font)

                for _, label_item in self.legend_y1.items:
                    label_item.setText(label_item.text, size=f'{Model.OnlinePlotter.font_size_pt}pt')

            if self.y2_items is not None:
                right_axis = self.ax1.getAxis('right')
                right_axis.setStyle(tickFont=tick_font)
                right_axis.setLabel(self.y2label, color='#c0c0c0', **label_style)

                if self.legend_y2 is not None:
                    for _, label_item in self.legend_y2.items:
                        label_item.setText(label_item.text, size=f'{Model.OnlinePlotter.font_size_pt}pt')
        
        def log_step(self, time):
            
            self.current_step += 1

            self.x_data = np.roll(self.x_data, -1)
            self.x_data[-1] = time

            y1vals = np.zeros(len(self.y1_items))
            if self.y2_items != None:
                y2vals = np.zeros(len(self.y2_items))
            # Append new data
            for j,yval in enumerate(self.y1_items):
                y1vals[j] = yval.v
            self.y1_data = np.roll(self.y1_data, -1, axis=1)
            self.y1_data[:,-1] = y1vals[:]
            
            if self.y2_items != None:
                for j,yval in enumerate(self.y2_items):
                    y2vals[j] = yval.v
                self.y2_data = np.roll(self.y2_data, -1, axis=1)
                self.y2_data[:,-1] = y2vals[:]

            if self.current_step % self.update_every == 0:
                self.__refresh_plot()

        def __refresh_plot(self):
            """
            Fast update using pyqtgraph setData
            """
            # Update y1 lines
            for i in range(len(self.y1_lines)):
                self.y1_lines[i].setData(self.x_data, self.y1_data[i,:])
            
            # Update y2 lines
            for i in range(len(self.y2_lines)):
                self.y2_lines[i].setData(self.x_data, self.y2_data[i,:])

            # Process GUI events (much faster than matplotlib canvas.draw)
            self.app.processEvents()
    # ------- end OnlinePlotter --------------------------------------------------

    # ----------------------------------------------------------------------------
    # Model class methods --------------------------------------------------------
    # ----------------------------------------------------------------------------
    def __init__(self):
        self.__components = []
        self.settings = Model.Settings()
        self.is_initialized = False
        self.is_first_step = True
        self.is_first_iteration = True
        self.is_converged = False
        self.time = 0
        self.plotters = []
        return
    
    def add_plotter(self, y1, y2=None, y1lim=None, y2lim=None, y1label='', y2label='', nmax_points = 1000, update_every=1, tab_title=None):
        if not isinstance(y1, type([])):
            y1t = [y1]
        else:
            y1t = y1
        y2t = y2
        if y2 != None:
            if not isinstance(y2, type([])):
                y2t = [y2]

        self.plotters.append(Model.OnlinePlotter(y1t, y2t, y1lim, y2lim, y1label, y2label, nmax_points, update_every, tab_title=tab_title))

    def wait_for_plots(self):
        """
        Call to keep the plotter window open after simulation is complete.
        """
        app = Model.OnlinePlotter.app
        main_window = Model.OnlinePlotter.main_window
        # Check whether plotting was ever initialized
        if app is None or main_window is None:
            return
        # bring the plot window to the front if it is minimized or behind other windows
        if main_window.isMinimized():
            main_window.showNormal()
        # raise and activate the window to bring it to the front
        main_window.raise_()
        main_window.activateWindow()

        # start the Qt event loop to display the plot window and allow interaction
        app.exec()

    def connect(self, source, destination, tol_rel = 1.e-6, tol_abs=1.e-6, log_n_iter = 0, learn_rate = 1.):
        
        # In order for connections to be established correctly, the model 
        # must be initialized first to assign names to all inputs and outputs. 
        # If the model is not initialized, initialize it now.
        if not self.is_initialized:
            self.initialize()

        if not isinstance(source, Component.Output):
            raise RuntimeError(f"Source connection object must be of type 'Component.Output'")
        if not isinstance(destination, Component.Input):
            raise RuntimeError(f"Destination connection object must be of type 'Component.Input'")

        destination.connection = Connection(source, tol_rel, tol_abs, log_n_iter, learn_rate)
        destination.is_connected = True
        source.is_connected = True
        return

    def initialize(self):
        """
        Initialize the model before running the simulation. This will automatically call the
        presim_setup() method of each component, and assign names to all inputs and outputs.

        The components are given access to the model under the attributed <component>.model, 
        which can be used to access settings and model flags.
        """
        if self.is_initialized:
            print('Model is already initialized. No further action taken. Calling function is:' + self.__class__.__name__   )
            return

        # Set the initial simulation time and iteration count
        self.time = self.settings.start_time
        self.iteration = -1
        self.is_first_step = True
        

        # Initialize the list of outputs. This is extended based on the component settings
        output_names = ['time','iterations']

        # ------------------------------------------------------------
        # Loop through all attributes of the model to find components, assign names, and call presim_setup
        for item in dir(self):

            itemobj = getattr(self, item)

            # If the object is a component...            
            if isinstance(itemobj, Component):
                # Give all of the component instances access to the model
                itemobj.model = self
                itemobj.name = ''

                # Handle component presim_setup here
                itemobj.presim_setup()
                # Automatically assign names to all inputs and outputs of this component 
                cnames = itemobj.auto_assign_names()
                # Add this component to the model's list of components
                self.__components.append(itemobj)
                # Record all output data names for the historian
                output_names += cnames

        # Construct the historian database
        nstep = int((self.settings.stop_time - self.settings.start_time)/self.settings.timestep)
        self.historian = dict([[n,np.ones(nstep)*float('nan')] for n in output_names])
        
        # Mark the model as initialized
        self.is_initialized = True
        return 

    # ----------------------------------------------------------------------- 
    def step(self, n_steps = 1):
        # Check overall initialization for the model
        if not self.is_initialized:
            self.initialize()
        
        self.iteration = -1  # reset the current iteration
        self.is_first_iteration = True
        self.is_converged = False

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
            
            if all_converged and not self.is_first_iteration:
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
            self.is_first_iteration = False

        # Convergence complete for this step, now do post-convergence calculations and logging
        self.is_converged = True
        current_step = int( (self.time - self.settings.start_time)/self.settings.timestep )
        for component in self.__components:
            component.calculate()
            for input in component.get_inputs(connected_only=True):
                input.connection.reset_for_step()

            self.historian['time'][current_step] = self.time
            self.historian['iterations'][current_step] = self.iteration
            for output in component.get_outputs():
                try:
                    if not isinstance(output.v, np.ndarray):
                        self.historian[output.name][current_step] = output.v
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

