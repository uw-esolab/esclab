"""
This file serves as a template for a Component library.

The file is stuctured to include definitions for multiple Components.

To make this library available to esclab, add the following line to src/esclab/components/__init__.py:
from .template_component import *
"""

# import necessary simulation tools
from esclab.simulate import *


# ---------------------------------------------------------
# Define component classes. Each component is a class that inherits from Component, and defines its own parameters, inputs, outputs, and calculate() method.


class FirstComponent(Component):
    """Template component - a simple calculator"""
    # Parameters. 
    #   Any parameter not given a default value must be set by the user before running the simulation.
    a = Component.Parameter()  # multiplier
    b = Component.Parameter()  # offset
    # Inputs. 
    #   The default input initial value is 1, but this can be passed as an argument here, or changed by directly setting the input value to a different number before running the simulation.
    x_in = Component.Input()  # input variable
    # Outputs. 
    #   The output value will be whatever is assigned in calculate().
    y_out = Component.Output()  # output variable
    z_out = Component.Output()  # another output variable

    def __init__(self):
        # [optional] 
        # initialize any internal variables here. Items stored as attributes of self will persist across time steps and iterations, and can be used to store intermediate values or state information. 
        super().__init__()  # call the parent class constructor [required if you define __init__()]
        self.x_last = 0.  # example of an internal variable to store the value of the input from the last time step  

    def presim_setup(self, **kwargs):
        super().presim_setup(**kwargs)
        # [optional] do any setup before the simulation starts here. This method is called once after the model is initialized, and before the first time step. You can use this to set up any internal variables or perform any calculations that are needed before the simulation starts. This differs from __init__ in that the model and all components are guaranteed to be fully initialized at this point, so you can access other components and their parameters and outputs if needed.
        return

    def calculate(self):
        # [required] define the equations that govern the behavior of the component here. This method is called at every time step and iteration.

        # super().calculate()  # << call this if inheriting from a parent class that also uses calculate(). Not needed otherwise.

        # ---------------------------------------
        # Do any first timestep calculations
        if self.is_first_timestep:
            # do any calculations that are needed at the first time step here. For example, you might 
            # want to set the output to a specific value at the first time step, or calculate an initial 
            # condition based on the input value at the first time step. 
            # return  # << uncomment to skip the rest of the calculations for the first time step
            pass  # continue out of this block if regular calculations are needed at the first time step

        # ---------------------------------------
        # Do post-convergence calculations 
        if self.is_converged:
            # do any calculations that are needed after convergence at each step(). For example, you 
            # might want to update the output based on the input and parameters after convergence is 
            # reached, or calculate some additional outputs that are only needed after convergence. 
            self.x_last = self.x_in  # store the input value for use in the next time step
            return  # skip the rest of the calculations if needed
        
        
        
        # ---------------------------------------
        # Do regular calculations for each iteration here. This is where the main behavior of the component
        # should be defined. For example, you might want to define the output as a function of the input 
        # and parameters here.

        self.y_out = self.a * self.x_in + self.b  # example calculation for y_out
        self.z_out = self.x_in - self.x_last  # example calculation for z_out
        

        # ---------------------------------------
        # If using the simultaneous equation solving feature, add equations to self.coupled_eqs here instead 
        # of directly calculating the outputs. This allows the equations to be solved simultaneously with 
        # other components in the same solve group. For example:
        if self.coupled_eqs is not None:  # check if we are in the matrix-build phase
            # Note that we use self.coupled_eqs.source(self.x_in) to refer to the input variable in the 
            # equation, which allows the solver to correctly link it to the output of whatever component 
            # is connected to this input. Point to the output variable directly (self.y_out) since it is 
            # an output of this component.

            # add an equation for y_out: y_out = a * x_in + b  -->  y_out - a*x_in = b
            self.coupled_eqs.add_equation({
                self.y_out: 1.0,
                self.coupled_eqs.source(self.x_in): -1.0,
            }, rhs=self.b)
            #  add an equation for z_out: z_out = x_in - x_last  -->  z_out - x_in = -x_last
            self.coupled_eqs.add_equation({
                self.z_out: 1.0,
                self.coupled_eqs.source(self.x_in): -1.0,
            }, rhs=-self.x_last)

            # The syntax of the add_equation() method is as follows:
            #     {var1: coeff1,      | dictionary of variables and their 
            #      var2: coeff2,      | coefficients in the equation
            #      ...},              | 

            #     rhs=value           | right-hand side value of the equation 
            #                         | (any constants that are not multiplied 
            #                         | by a variable should be included here)
# ----------------------------------------------------------

# ----------------------------------------------------------
class SecondComponent(Component):
    """Another example component. This one doesn't do anything, but serves as an example of how to define multiple components in the same file."""
    x_in = Component.Input()
    y_out = Component.Output()
    def calculate(self):
        if self.coupled_eqs is not None:  # check if we are in the matrix-build phase
            # add an equation for y_out: y_out = x_in  -->  y_out - x_in = 0
            self.coupled_eqs.add_equation({
                self.y_out: 1.0,
                self.coupled_eqs.source(self.x_in): -1.0,
            }, rhs=0.0)
# ----------------------------------------------------------
