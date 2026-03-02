---
description: Load when converting TRNSYS Fortran types to Python esclab types. 
---

# Your role
You are a software engineer with a mechanical engineering (thermo-fluids) background and expertise in both Fortran and Python, particularly in the context of engineering simulations. Your task is to convert TRNSYS Fortran code into Python code for the esclab project.  


# Background
TRNSYS is a widely used software for simulating the behavior of transient systems, particularly in the field of energy and building performance. It is written in Fortran 90. TRNSYS is organized into: 
* TYPES: These are the building blocks of TRNSYS. They represent components such as solar collectors, heat exchangers, and storage tanks. Each TYPE has its own set of equations and parameters that define its behavior. 
* Decks: These are the input files that define the system being simulated. They specify which TYPES are used, how they are connected, and the parameters for each TYPE.

A major weakness of TRNSYS is that functionality is limited by Fortran and the limited set of libraries available, requiring basically everything to be written explicitly. TRNSYS also uses a lot of global variables, which can make it difficult to understand and maintain the code. Finally, TRNSYS lacks modern data structures, including classes, which requires special memory handling for TYPES that have multiple instances and are called multiple times in the same simulation.

## esclab structure

esclab structures are defined in the 'simulate.py' file, which is the main entry point for running simulations. The 'components' folder contains the code for the various components that can be used in simulations. The 'models' folder contains assembled components that are connected together to create a modeled system. TYPES-->components, Decks-->models. 

The eeslib library is a Python library that provides a collection of functions and classes for energy systems modeling and simulation. It is designed to be used in conjunction with esclab to provide additional functionality and support for more complex simulations. The eeslib library provides correlations, property calculations, unit conversions, and other utilities that are commonly needed in energy systems modeling. The eeslib library mimics some of the functionality provided by Engineering Equation Solver (EES), which is a popular software for solving engineering equations. In cases where converted code could be better suited for implementation in eeslib instead of esclab, the code should be implemented in eeslib and then called from esclab. 

esclab components contain:
* Inputs: These are the variables that are passed into the component from other components or from initial values provided in the model file. These are generally updated many times throughout the simulation.
* Outputs: These are the variables that the component calculates and passes back to other components or to the model file. They represent the results of the component's calculations.
* Parameters: These are the variables that define the behavior of the component. They can be set in the model file and can be used to customize the component's behavior for different simulations. They generally should not change during the simulation.

All Inputs/Outputs/Parameters use the Component.__io_base class which stores the data value in the member 'v'. It is generally necessary to reference any instance like this: self.input_name.v, self.output_name.v, self.parameter_name.v.

## Steps in the simulation process
1. The model file imports the necessary components and creates instances depending on the system being modeled. 
2. Connections between components are defined in the model file using the 'connect' function. This function specifies which outputs from one component are connected to which inputs of another component.
3. Data processing and plotting tools are declared in the model file.
4. Simulation settings are configured in the model file, including the time step, total simulation time, and any other relevant parameters defined in Model.settings.
5. The simulation is run by calling the 'step()' member of the model for each time step over the time horizon. During each time step, the following occurs:
    a. Each component is checked for initialization and will have its 'initialize()' method called if necessary.
    b. In a loop over each component instance:
      i. The connections between inputs and outputs of components are updated. 
      ii. The component's 'calculate()' method is called
        - calculate() is called at each iteration. Calculations can be differentiated based on the first timestep, the first iteration of each timestep, or any other iteration condition.
      iii. The component connections are checked for convervence. 
      iv. If all connections have converged, the 'converged()' method of each component is called one final time. Final computed values are stored in output data members.


# Objectives of esclab
esclab ("engineering simulation and controls lab") is a Python-based teaching tool that I'm writing to replace TRNSYS in my courses. The main objectives of esclab are to:
1. Provide a more modern and user-friendly interface for students to learn about energy systems modeling and simulation.
2. Leverage the power of Python and its libraries to enable more complex and realistic simulations. The use of speed-enhancing libraries such as NumPy is of particular importance and should be a major focus of the conversion process.
3. Encourage good coding practices and software engineering principles in the context of energy systems modeling.
4. Provide accessible and well-documented code that can be easily modified and extended by students.



# Mapping

## RUN TYPE 1

Please refer to the directory ./src/esclab/components/flownetwork/fortran-sources for the Fortran source files.

Source files are identified by the Type number. Generally, the file name contains the number, but sometimes this can instead by found on the first line of the file with the Subroutine definition. 

We are focusing on converting the following types:
4001
4004
4006
4007
4008
4012
4015
4016
4034
4035
4050
4097
4100
4101
4102
6001
6003
6007
6011
6014
6016
6017
6019
6022
6027
6028
6030
6031
6032
6034

You must follow the instructions below when converting the Fortran code to Python, in this order. Do not execute the code and test for output values until explicitly stated.

1. Using subagents assigned each to one file, starting with Type 4001, create a new python file in flownetwork corresponding to the type. Name the file with a descriptive name that reflects the functionality of the component. For example, Type 4001 is a simple pipe flow component, so it could be named 'PipeFlow.py'. The file name and component class should be named using PascalCase and should also reflect the functionality of the component, such as 'PipeFlow'. Modify the __init__.py file in the components folder to import the new component class.
2. Proceed line by line, making only direct conversions from Fortran to Python. Do not change the structure of the code. The goal is direct mapping and the code would fail if executed. In this step, you should handle the following:
  a. Map parameters and inputs declared at the top of the Fortran file to Component members. the getParameterValue(<number>) function should map to the Component.Parameter() member that corresponds to the parameter number in the Fortran code. The same applies for getInputValue(<number>). Do not use parameter_<num> or input_<num> variable names. Instead, use the corresponding variable name directly as the parameter name. 
  b. Declare outputs as Component.Output() members. You will need to identify appropriate output names based on the SetOutputValue() function calls near the end of the fortran file.
  c. Replace syntax like exponentiation with **, array indexing with [], and function definitions with def. 
  d. Replace if statements and loops with python syntax. Mark any conversions that are not straightforward with a comment like "# TODO-NEEDS CONVERSION REVIEW: " and a brief description of the issue.
  e. Never introduce value clamping, min/max functions, safety fallbacks, or other modifications that change the behavior of the code. The goal is a direct mapping in this step.
  f. Leave missing libraries and functions alone and mark with a comment like "# TODO-NEEDS LIBRARY: " and a brief description of the library or function needed. For example, if there is a call to a function in the FIT library for fluid properties, mark it with a comment like "# TODO-NEEDS LIBRARY: FIT library for fluid properties".
  g. Port all comments from the original Fortran code. If a line is <fortran code> !comment, the comment should be ported as a python comment on the line above the code. If there are block comments, these should be ported as block comments in Python. Do not add conversion comments like "Same structure as fortran" aside from explicitly stated TODO's in these instructions, but add comments provided in the original fortran source that explain the purpose of the code blocks.
3. Move code blocks that are inside the 'getIsStartTime()', 'getIsEndOfTimestep()', 'getIsFirstCallofTimestep()', 'getTimestepIteration()==0', and 'getCurrentTime()' blocks to if statements that check for the appropriate model flags or members. For example, code inside the 'getIsStartTime()' block should be moved inside an if statement checking for the model.is_first_step flag.
4. Convert property library calls. HTF properties map to esol_properties.Incompressible. Water properties map to eeslib.fluid_properties. Function arguments must be specified with the parameter name (e.g., T=Teval, P=Peval, etc.). Do not use surrogates for properties.  The fluid_id in fortran is a number, but in the new implementation, assume fluid_id is a string that is directly passed on. Do not use "fallback" fluid names like "Nitrate salt". Assume the property functions check for validity and return float values, and do not attempt to convert the return type or clamp/clip/limit input arguments. Flag any suspected units mismatches (e.g., J->kJ) with a comment like "# TODO-NEEDS UNITS CHECK: " and a brief description of the issue.
5. After direct conversion, note that there are some redundancies in setting and using parameters, inputs, and outputs. Find instances where a local scope variable is assigned the value of the Component.<class> member, and prefer instead to directly use the Component.<class> member in the code. For example, if there is a line like "param1 = self.myparameter.v", and then param1 is used in the code, it would be better to directly use "self.myparameter.v" instead of creating a new variable. Make these changes throughout the code, but do not change any of the underlying logic or structure of the code except to remove these redundancies. 

Do not proceed to Run Type 2 or 3 unless I explicitly ask you to do so.

## RUN TYPE 2

This step applies only to types that have already been converted and not to skipped (too complex) types.

In this step, follow the same steps as in RUN TYPE 1, except do not pursue type conversion but instead focus on the missing helper functions called in fortran-source/solar_field_modules.f90 and in other fortran files. Note property calculations are implemented multiple times and should not be converted. Instead, revert to eeslib or esol_properties as appropriate. 

Do not change anything in the existing converted TYPE files EXCEPT to add calls to the new helper functions where appropriate. Be sure to go back to the original Fortran context to understand how to properly call the helper functions in the context of the existing code.

Do not proceed to Run Type 3 unless I explicitly ask you to do so. Do not run tests to check output values.


## RUN TYPE 3

Common TODO flags across files to review:

Resolve the following todo's:
* TODO-NEEDS UNITS CHECK — kPa↔Pa and kJ/kg↔J/kg at eeslib call sites
* TODO-NEEDS LIBRARY — Not sure which of these remain
* TODO-NEEDS CONVERSION REVIEW — dynamic array storage patterns, variable input/output count (Types 4050, 6027)

Review units noting that esclab will always use base SI (K, Pa, J, kg, s). Where clues exist that original fortran units may have been different, correct the units with an appropriate convert() call:
* non-temperature: *convert('<non-SI>', '<SI>')
* temperature: converttemp('<non-SI>', '<SI>', <temp_value>)
Mark implemented units conversions with a comment like "# AUTO UNITS CONVERSION IMPLEMENTED: <description of conversion>" for clarity.