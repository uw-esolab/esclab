# ESCLAB
**Engineering Simulation and Control Lab**

This software provides an environment for time-series simulation and control of complex multi-component systems. The approach used here relies on iterative (successive substitution) solving to converge on values for a given step, allowing compartmentalization of calculations between individual system components. 

The simulation is configured and managed in Python scripts. 

A graphical plotting interface provides real-time simulation results and a network topology rendering tool.

The software is primarily intended for research and teaching activities at the University of Wisconsin-Madison.

## Installation

1. Download and install a Python package manager. This program was developed using [miniconda](https://www.anaconda.com/download/success?reg=skipped), and it's recommended for environment creation and maintenance. 

2. Download and install [VS Code](https://code.visualstudio.com/download) (or your favorite Python IDE)

3. Collect the ESCLab package source code. If you're using git, clone the repository. Otherwise, you can download and unzip the repository. 
```bash
# Change to a convenient location for the repository. A path with no spaces is recommended.
cd C:\repositories
# Clone
git clone https://github.com/uw-esolab/esclab.git
```

4. Open a command console and execute the following:
```bash
# Change to the library directory
cd C:\repositories\esclab
# Create and activate a new conda environment
conda create -n esclab_dev python=3.13
conda activate esclab_dev

# Install the esclab project to the python package you just created. The -e flag makes the installation editable.
pip install -e .
# If you are also modifying the eeslib dependency (published by uw-esolab), then clone and install the editable version to this conda environment.
pip install -e C:\repositories\eeslib
# **otherwise** install using the standard pip call
# pip install eeslib
```

## Getting started

* Source code for constructing and simulating models is found in the /src/esclab directory
* Components are defined in the */src/esclab/components* folder
* Models are defined in the */src/esclab/models* folder

To create a new model, prepare one or more component files that contain instance(s) of Component classes that are units in the system you want to model. Connect components together in a *model* file.

Templates and examples for *Components* and *Models* are provided in their respective folders. 
* `template_model.py` provides a template for building new system models
* `template_component.py` provides a template for a component library
* `sample_circuit_secsub.py` gives and example of a simple solve loop using successive substitution (no coupled equations)
* `sample_circuit.py` gives an equivalent example using the coupled equation approach. `sample_circuit_tee.py`, and `sample_circuit_hilopass.py` also give examples of solving coupled equations for increasingly complex topologies.

The main calculation and plotting scripts are:
* `simulate.py` | core simulation engine
* `network_topology.py` | tools for detecting, organizing, and rendering coupled network systems
* `online_plotter.py` | Qt-based window for real-time plotting 

The `flownetwork` folder is currently under-development code for a CSP thermal-hydraulic simulation. None of it is functioning in this environment!!

Definitions:
| Term               | Definition                                                                                                                       |
| ------------------ | -------------------------------------------------------------------------------------------------------------------------------- |
| Component          | Building block for a system model. Components take input and compute output values, and they connect with other components       |
| Model              | Collection of one or more components together into a simulated system                                                            |
| Simulation         | Time-dependent evaluation of a model                                                                                             |
| Step               | Calculations done at a specific value in time. A simulation contains many steps |
| Coupled equations  | Equations that are designated to be solved together using matrix inversion, outside of the normal successive-substitution method |
| Absolute tolerance | Absolute difference between connection values on succesive iterations                                                            |
| Relative tolerance | Difference relative to the magnitude of the last connection value                                                                                                                                 |
| Learning rate      | Fraction of the difference between the new and old computed values to apply when iterating                                                                                                                                 |

## Authors

Lead author:
Mike Wagner | Associate Professor, University of Wisconsin-Madison | [GitHub: uw-esolab](https://github.com/uw-esolab/) | [Profile](https://engineering.wisc.edu/directory/profile/mike-wagner/)

Contributors:
