"""
Top-level pytest fixtures shared across all esclab test modules.
"""

import pytest


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def set_inputs(component, **kwargs):
    """
    Convenience helper: set one or more Input values on a component.

    Usage::

        set_inputs(pipe, m_dot_in=1.0, T_in=150.0, P_in=5e5)
    """
    for name, value in kwargs.items():
        getattr(component, name).v = value


def set_params(component, **kwargs):
    """
    Convenience helper: set one or more Parameter values on a component.

    Usage::

        set_params(pipe, Pipe_ID=0.05, Length_Pipe=10.0)
    """
    for name, value in kwargs.items():
        getattr(component, name).v = value
