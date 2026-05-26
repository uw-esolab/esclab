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

You can find the TRNSYS developers guide specifying TRNSYS functions and documentation in src/esclab/components/flownetwork/fortran-source/07-ProgrammersGuide.pdf.

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
23
162

Not to be converted but are used in the deck:
65
57


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
3. Move code blocks that are inside the 'getIsStartTime()', 'getIsEndOfTimestep()', 'getIsFirstCallofTimestep()', 'getTimestepIteration()==0', and 'getCurrentTime()' blocks to if statements that check for the appropriate model flags or members. For example, code inside the 'getIsStartTime()' block should be moved inside an if statement checking for the model.is_first_step flag. Ignore code inside 'getIsIncludedInSSR()' blocks.
4. Convert property library calls. HTF properties map to esol_properties.Incompressible. Water properties map to eeslib.fluid_properties. Function arguments must be specified with the parameter name (e.g., T=Teval, P=Peval, etc.). Do not use surrogates for properties.  The fluid_id in fortran is a number, but in the new implementation, assume fluid_id is a string that is directly passed on. Do not use "fallback" fluid names like "Nitrate salt". Assume the property functions check for validity and return float values, and do not attempt to convert the return type or clamp/clip/limit input arguments. Flag any suspected units mismatches (e.g., J->kJ) with a comment like "# TODO-NEEDS UNITS CHECK: " and a brief description of the issue. Humid air properties in TRNSYS are calculated with the MoistAirProperties function. eeslib contains a call humid_air returning a dictionary of desired outputs (e.g, T_wb, T_db, RH) that should be used instead.
5. After direct conversion, note that there are some redundancies in setting and using parameters, inputs, and outputs. Find instances where a local scope variable is assigned the value of the Component.<class> member, and prefer instead to directly use the Component.<class> member in the code. For example, if there is a line like "param1 = self.myparameter.v", and then param1 is used in the code, it would be better to directly use "self.myparameter.v" instead of creating a new variable. Make these changes throughout the code, but do not change any of the underlying logic or structure of the code except to remove these redundancies. 

Do not proceed to other Run Types unless I explicitly ask you to do so.

## RUN TYPE 2

This step applies only to types that have already been converted and not to skipped (too complex) types.

In this step, follow the same steps as in RUN TYPE 1, except do not pursue type conversion but instead focus on the missing helper functions called in fortran-source/solar_field_modules.f90 and in other fortran files. Note property calculations are implemented multiple times and should not be converted. Instead, revert to eeslib or esol_properties as appropriate. 

Do not change anything in the existing converted TYPE files EXCEPT to add calls to the new helper functions where appropriate. Be sure to go back to the original Fortran context to understand how to properly call the helper functions in the context of the existing code.

Do not proceed to other Run Types unless I explicitly ask you to do so. Do not run tests to check output values.

## RUN TYPE 3

Deck file conversion

Refer to src/esclab/models/sample_air-brayton.py for an example of a esclab model file.

The TRNSYS .dck (deck) files specify instances of each of the TYPES and assign a unique UNIT number to each instance. The deck files also specify the connections between the inputs and outputs of the different instances, as well as the parameters and initial input values for each instance. Deck files will be converted to Python model files that specify the same information but in the context of the esclab structure. Convert using the following information:
1. New model files go in the 'models' folder with a descriptive name that reflects the system being modeled. For example, if the deck file is modeling a solar thermal system, the model file could be named 'SolarThermalSystem.py'. The model class should also be named using PascalCase and should reflect the system being modeled, such as 'SolarThermalSystem'. Modify the __init__.py file in the models folder to import the new model class.
2. Whole-line comments start with '*' characters. End-of-line comments start with '!' characters. Port all comments to the new Python file.
3. The Control Cards section specifies overall simulation configuration settings. Map settings to the Model.settings member variables as appropriate. Some settings may not map and can be ignored, but should be marked in comments with "# TODO-NEEDS CONVERSION REVIEW: " and a brief description of the setting and why it may not map.
4. Each UNIT instance contains the following subitems. Delegate the mapping of UNITS/TYPES to subagents assigned to each type.
  a. A declaration of the TYPE and UNIT number. The TYPE number should be used to identify which component class to instantiate. The UNIT number is referenced throughout the file in other type connections to specify input/output connections. There are sometimes comments like "* Model <descriptive name>" with a unit, and those should be used to help name the component instance in a descriptive way. For example, if there is a line like "UNIT 1 TYPE 4001 * Model Pipe 1", the component instance could be named "pipe1 = PipeFlow()". If there is no descriptive comment, use the type and unit number to create a name like "type4001_unit1 = PipeFlow()".
  b. A parameter block starting with "PARAMETERS <number>" indicating that the <number> of lines that follow map to the required parameter list for the TYPE. These should be mapped to the Component.Parameter() members of the corresponding component instance by referencing both the original Fortran TYPE and the converted Python class. 
  c. INPUTS <number> block indicating the number of lines that follow mapping to the required input list for the TYPE. These should be mapped to the Component.Input() members of the corresponding component instance by referencing both the original Fortran TYPE and the converted Python class. Initial input values specified in the deck file should be assigned to the .v member of the corresponding Component.Input() member. The format of each input is <unit number>,<input number>    ! <comment> [sometimes with useful connection information]. Unconnected inputs are denoted 0,0 and should be omitted from conversion. Properly connecting inputs will require that you i) identify the input from the numerical ordering of inputs for a TYPE in the original fortran file, ii) identify the output similarly, and iii) identify the matching inputs/outputs in the converted python files.
  d. Initial values for each input are important. Unconnected values use the initial value specified in the deck file throughout the simulation. Connected values use the initial value only for the first time step, and then are updated based on the connection. Be sure to properly assign initial values to the .v member of the corresponding Component.Input() member. The order of initial values given in the deck file matches the order of inputs. In esclab, the initial value of inputs should be specified in a block after specifying the parameter values.
  d. Optional LABELS block that sometimes specifies the file name of an input or output data file. If mapping isn't obvious, mark with a comment like "# TODO-NEEDS CONVERSION REVIEW: " and a brief description of the issue.
5. Equations blocks in deck files can generally be converted to python script directly. However, these equations are very simple and do not preserve order of operations, so some human review will be needed. In general, err on the side of preserving the equation syntax as it appears in the deck file.
6. Ignore plotters and data file loggers. Ignore anything after the "END" statement, including *!LINK... stuff.

Conversion steps:
1. Create new Python file
2. Generate a list of the components with their unique names and save the file. This will be used as a reference for the next steps. Note any information that you'd need to resume to resume at step 3 below.
3. For each type, create separate file with a mapping template that identifies the parameters, inputs, and outputs for the type based on the original Fortran file. This will be used to map the parameters, inputs, and outputs in the deck file to the corresponding members of the component instances in the Python model file. Structure this template in a way that helps you expedite the mapping process in step 4.
4. Map the parameters, inputs, and outputs for each component instance in the deck file to the corresponding members of the component instances in the Python model file using the mapping template created in step 3.

## RUN TYPE 4

Common TODO flags across files to review:

Resolve the following todo's:
* TODO-NEEDS UNITS CHECK — kPa↔Pa and kJ/kg↔J/kg at eeslib call sites
* TODO-NEEDS LIBRARY — Not sure which of these remain
* TODO-NEEDS CONVERSION REVIEW — dynamic array storage patterns, variable input/output count (Types 4050, 6027)

Review units noting that esclab will always use base SI (K, Pa, J, kg, s). Where clues exist that original fortran units may have been different, correct the units with an appropriate convert() call:
* non-temperature: *convert('<non-SI>', '<SI>')
* temperature: converttemp('<non-SI>', '<SI>', <temp_value>)
Mark implemented units conversions with a comment like "# AUTO UNITS CONVERSION IMPLEMENTED: <description of conversion>" for clarity.

## RUN TYPE 5

This step refactors a directly-converted component (from RUN TYPE 1–4) into a fully native esclab component. It resolves the two main structural problems that remain after a direct conversion: (1) iterative operating-point solvers that belong in the network equation matrix, and (2) state that should live in instance variables rather than Output ports.

**When to apply:** a component has any of the following:
- A secant, Newton, or bisection loop inside `calculate()` to find a flow/pressure operating point
- Output ports used as "previous iteration" memory (e.g., `P1_point_1x`, `P1_point_2x`)
- Many scalar Parameters that are logically a grouped array (e.g., coefficients per pump)
- State that must persist across timesteps (tank mass, pressure, enthalpy, level)

**Do not apply until** all RUN TYPE 1–4 TODO flags are resolved and `get_errors` is clean.

### Rules

**Rule 1 — Read the framework before writing any new code.**
Open `simulate.py` and read the `Component` base class and the `NetworkEquationContext` (`coupled_eqs`) API. Open `circuit_elements.py` and read `Capacitor` and at least one `Pump` class as concrete examples. You need to understand:
- Lifecycle: `presim_setup(**kwargs)` → `calculate()` (iterated) → `is_converged` block inside `calculate()`
- `self.model.is_first_step`, `self.model.is_first_iteration`, `self.model.is_converged`
- `self.coupled_eqs`: `None` for sequential components; a `NetworkEquationContext` for components inside a coupled subnetwork
- `coupled_eqs.add_equation(terms_dict, rhs)` and `coupled_eqs.source(input_port)`

**Rule 2 — Classify every variable into exactly one category before touching the code.**
- *Instance variables* (`self._name`): state that must survive from one timestep to the next (e.g., tank pressure, enthalpy, mass, level, fluid density, temperature). Initialised in `presim_setup()`, advanced in the `is_converged` block.
- *Output ports* (`self.name.v`): values consumed by downstream components or the model file. Also serve as linearisation points between iterations (their values from the previous iteration are read at the top of `calculate()`).
- *Local variables*: intermediate quantities that do not need to survive the call.
Remove all Output ports that served only as solver-state memory (secant-method bracket values, previous-iteration guesses, etc.).

**Rule 3 — Collapse grouped scalar Parameters into array Parameters.**
If a component has N identical sub-elements (pumps, HX passes, etc.) each described by the same K coefficients, replace the N×K scalar Parameters with a single `ndarray` Parameter of shape `(N, K)`. Document the shape, axis ordering, and units in the class docstring. A 3-pump component with [A, B, C] head coefficients becomes one `pump_head_coeffs` Parameter of shape `(3, 3)`.

**Rule 4 — Write `presim_setup()` first, before `calculate()`.**
- Compute initial state from Parameter initial-condition values using property lookups (eeslib, esol_properties).
- Assign results to all `self._*` instance variables.
- Assign physically meaningful initial values to *every* Output port. Never leave a port at 0.0 for a quantity like pressure or enthalpy — use steady-state estimates.
- Call helper methods (e.g., alarm checks) here too so alarms are correct on the first timestep.

**Rule 5 — Extract pure helper methods for reusable sub-calculations.**
Functions that do not touch `self.model`, `self.coupled_eqs`, or Output ports make clean helpers. Alarm/trip logic is a common candidate: `def _check_alarms(self, level): -> (LL_Alarm, LL_Trip, HL_Alarm, HL_Trip)`. Call these from both `presim_setup()` and `calculate()`.

**Rule 6 — Structure `calculate()` in three explicit phases.**

*Phase A — matrix contribution* (guard: `if self.coupled_eqs is not None:`):
  - Read current Output port values as linearisation points.
  - Linearise each nonlinear constraint and call `self.coupled_eqs.add_equation(terms_dict, rhs)` once per equation.
  - Count equations vs. unknowns before coding to confirm the local system is square.
  - End Phase A with a bare `return`.

*Phase B — sequential diagnostics* (no guard; runs only when `coupled_eqs is None`):
  - The network matrix has already written pump/flow Output port values; read them directly.
  - Compute bleed extraction, inlet mixing, per-element power/efficiency, vent enthalpy, etc.
  - Write all Output ports needed by downstream components.
  - End with `if not self.model.is_converged: return`.

*Phase C — convergence block* (runs only on the converged iteration):
  - Perform ODE integration (RK4 or similar) for state variables.
  - Advance all `self._*` instance variables to end-of-timestep values.
  - Write final diagnostic Output ports (level, tank mass, alarm signals, cavitation trips).

**Rule 7 — Linearise hydraulic head curves correctly.**
For a speed-scaled head curve `ΔP = ρg(AQ² + BsQ + Cs²)` linearised around operating point `Q₀` (in m³/s):
- slope (Pa·s/kg): `g × (2A·Q₀ + B·s)`
- intercept (Pa): `ρg × (C·s² − A·Q₀²)`
- Matrix equation: `{P_out: 1, m_dot: −slope}` with `rhs = P_inlet + intercept`

Always clamp `Q₀ = max(m_dot_pi.v, 1e-6) / rho` to prevent the slope from collapsing to zero at start-up. For a constant-speed pump set `s = 1.0`.

**Rule 8 — Count and verify matrix equations before writing Phase A.**
List each unknown Output port the matrix must solve, then write exactly one equation per unknown. Common pattern for a multi-pump component with a shared discharge header and `N` pumps:
- N × (head-curve equation if pump is on, else zero-flow equation)
- 1 × flow summation: `m_dot_out = Σ m_dot_Pi`
- 1 × loop closure: `m_dot_out = source(m_dot_in)` (or static-pressure pin when all pumps are off)
- 1 × enthalpy: `h_out = h_inlet + specific_work`

Total: `N + 3` equations for `N + 3` unknowns (`m_dot_P1…PN`, `m_dot_out`, `P_out`, `h_out`).

**Rule 9 — Use `source(input_port)` for loop-closure equations.**
`self.coupled_eqs.source(self.m_dot_in)` returns a reference to the upstream node in the network graph so the matrix enforces flow continuity without hard-coding a value. Use this for the loop-closure equation whenever the component receives a flow that must equal its own discharge.

**Rule 10 — Write RK4 integration as a nested closure inside the Phase C block.**
Define `def _rk4_rates(P, h): ...` inside `calculate()` after the `is_converged` guard. The closure captures `m_dot_in`, `m_dot_pump`, `m_dot_vent`, and `m_tank` from the surrounding scope, keeping the four RK4 calls concise. Guard denominators with `math.copysign(max(abs(x), ε), x)` to match the stabilisation from the original Fortran code.

**Rule 11 — Delete old code in small, uniquely-identifiable chunks.**
- Insert the new methods first and leave the old body appended below them. The file will have syntax errors but `get_errors` will tell you exactly where.
- Delete old code one logical block at a time (one secant-solver step, one `if` branch, etc.) using `replace_string_in_file` with 3–5 lines of unchanged context on both sides.
- Run `get_errors` after each deletion.
- Never attempt to delete more than ~100 lines in a single call; match text only appears once in the file.
- Finish by running `pytest tests/` to confirm no regressions.