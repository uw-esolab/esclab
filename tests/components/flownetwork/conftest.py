"""
Fixtures shared across flownetwork component tests.

Fluid IDs (eeslib / TRNSYS convention)
---------------------------------------
Working fluid - Water
Heat transfer fluid (HTF) - Dowtherm A
Thermal storage salt - Solar Salt

Typical solar-field operating conditions used throughout these tests:
  T_in = 300 °C
  P_in = 500 000 Pa  (5 bar)
  m_dot = 1.0 kg/s
"""

import pytest

from esclab.components.flownetwork.SimplePipe import SimplePipe


# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

FLUID_NAME = 'Water'
HTF_NAME = 'Dowtherm A'
SALT_NAME = 'Solar Salt'  


# ---------------------------------------------------------------------------
# Standard boundary-condition fixtures
# ---------------------------------------------------------------------------

@pytest.fixture
def typical_flow_conditions():
    """Return a dict of standard HTF flow boundary conditions."""
    return dict(
        fluid=FLUID_NAME,
        T_in=300.0,    # °C
        P_in=5e5,      # Pa
        m_dot_in=1.0,  # kg/s
    )


# ---------------------------------------------------------------------------
# Component fixtures
# ---------------------------------------------------------------------------

@pytest.fixture
def simple_pipe():
    """
    Return a SimplePipe configured with representative solar-field parameters.

    Geometry
    --------
    Pipe_ID     = 0.05 m   (DN 50 nominal)
    Length_Pipe = 10 m
    Roughness   = 4.6e-5 m (commercial steel)
    """
    pipe = SimplePipe()
    pipe.Pipe_ID.v = 0.05
    pipe.Length_Pipe.v = 10.0
    pipe.Roughness.v = 4.6e-5
    pipe.fluid.v = FLUID_NAME
    # Set standard operating inputs before presim_setup
    pipe.m_dot_in.v = 1.0
    pipe.T_in.v = 300.0
    pipe.P_in.v = 5e5
    pipe.presim_setup()
    return pipe
