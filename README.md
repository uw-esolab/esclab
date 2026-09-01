# ESCLAB
**Engineering Simulation and Control Lab**

This software provides an environment for time-series simulation and control of complex multi-component systems. The approach used here relies on iterative (successive substitution) solving to converge on values for a given step, allowing compartmentalization of calculations between individual system components. 

The simulation is configured and managed in Python scripts. 

A graphical plotting interface provides real-time simulation results and a network topology rendering tool.

The software is structured to support both research and teaching activities and is used in the upper-graduate course "ME 964: Simulation and Optimal Control of Energy Systems" at the University of Wisconsin-Madison.

## Installation

1. Download and install a Python package manager. This program was developed using [miniconda](https://www.anaconda.com/download/success?reg=skipped), and it's recommended for environment creation and maintenance. 

2. Download and install [VS Code](https://code.visualstudio.com/download) (or your favorite Python IDE, though VS Code is currently recommended)

3. Open a command console and create a new conda environment. You can call this what you'd like, but I assume the name 'esclab_dev' in the documentation:

```bash
# Create and activate a new conda environment
conda create -n esclab_dev python=3.13
conda activate esclab_dev
```

### Installing from source code

Use this option when you plan on extending the component or model libraries with code that will closely integrate. This allows direct modification of the ESCLab source code and incorporation of the latest updates that are pushed to the GitHub repository. 

4. Change to a convenient location for the repository. A path with no spaces is recommended.
```bash
cd C:\repositories
```

5. Collect the ESCLab package source code. If you're using git, clone the repository. Otherwise, you can download and unzip the repository. 
```bash
# Clone
git clone https://github.com/uw-esolab/esclab.git
```

6. Install the ESCLab project to the python package you just created. 

```bash
# The -e flag makes the installation editable.
pip install -e .
# ESCLab currently requires the most recent version of EESLib (published by uw-esolab). 
# Clone and install the editable version to this conda environment. 
# Make sure C:\repositories is replace with the path to your repository directory.
pip install -e C:\repositories\eeslib
# Ensure the editable version is correctly installed in esclab_dev.
python -c "import eeslib; print(eeslib.__file__)"
```

### Installing from the package manager (PyPi)

Do this **instead** of the source code option if you're happy with a relatively stable release of the ESCLab core source code and will only be developing models and components for your own use without intending to share them as part of the ESCLab distribution. 

Following step #3 above: 

4. Install the package using `pip`:
```bash
pip install esclab
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

## Uploading a new version to PyPi (developers only!)

To upload a new version of ESCLab to Pypi, follow the steps outlined in the [Python packaging tutorial](https://packaging.python.org/en/latest/tutorials/packaging-projects/). 

The preferred packaging tool is `setuptools`. 

The most relevant steps are as follows:
1. Don't forget to update the code version number in `pyproject.toml` and in the `src/__init__.py` file.
1. Open a command window and navigate to the `esclab` directory, such as
    ```bash
    cd C:\repositories\esclab
    ```
2. Ensure the build and packaging tools are installed in the Conda environment that you're using. 
    The preferred method will install a developer tools, including build, pytest, and twine, specified in the pyproject.toml file in the esclab directory:
    ```bash
    pip install .[dev]
    ```
    Alternatively, you can manually install packages:
    ```bash
    python -m pip install --upgrade build 
    python -m pip install --upgrade twine
    ```
3. Build the Python distributable
    ```bash
    python -m build
    ```
    This should create a folder `dist/` that contains a wheel (.whl) and tar.gz file. 
4. Upload the file to Pypi. You will need to have first created a username and API token, following the packing tutorial instructions. If uploading to the production server, use the command:
    ```bash
    twine upload dist/*
    ```
    If you receive the warning `WARNING  This environment is not supported for trusted publishing`, you can ignore it. 

    If using the test server, use the command:
    ```bash
    python -m twine upload --repository testpypi dist/*
    ```

    To test installation from the test server in a new, temporary environment:
    ```bash
    conda create -n test_esclab python=3.13
    conda activate test_esclab
    pip install -i https://test.pypi.org/simple/ --extra-index-url https://pypi.org/simple/ esclab==0.0.1

    python -c "import esclab; print(esclab.__file__);"
    >> C:\Users\username\AppData\Local\miniconda3\envs\test_esclab\Lib\site-packages\esclab\__init__.py

    conda activate base
    conda env remove -n test_esclab
    ```

5. To install the package, activate the Conda environment (e.g., `conda activate esclab_dev`), and install. 
    If running from the production environment use:
    ```bash
    pip install esclab
    ```

    If you have previous versions of `esclab` already installed, force use of the most recent version using: 
    ```bash
    pip install --force-reinstall --upgrade esclab
    ```


## Authors

Lead author:
Mike Wagner | Associate Professor, University of Wisconsin-Madison | [GitHub: uw-esolab](https://github.com/uw-esolab/) | [Profile](https://engineering.wisc.edu/directory/profile/mike-wagner/)

Contributors:
