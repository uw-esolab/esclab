"""
Unit tests for SimplePipe (Type 4001).

Each test follows the standard component workflow:
  1. Configure parameters (done by the ``simple_pipe`` fixture).
  2. Set input values.
  3. Call ``calculate()``.
  4. Assert output values.

Note on class-level descriptors
---------------------------------
``Component.Input``, ``Output``, and ``Parameter`` are class-level attributes,
so all instances of a subclass share the same descriptor objects.  Tests must
not create multiple live instances of the same component simultaneously – the
``simple_pipe`` fixture already handles this correctly through function scope.
"""

import pytest

from esclab.components.flownetwork.SimplePipe import SimplePipe, FricFactor_IC


# =============================================================================
# presim_setup
# =============================================================================

class TestPresimSetup:
    """Outputs are initialised from current inputs after presim_setup()."""

    def test_mass_flow_output_mirrors_input(self, simple_pipe):
        assert simple_pipe.m_dot_out.v == pytest.approx(simple_pipe.m_dot_in.v)

    def test_temperature_output_mirrors_input(self, simple_pipe):
        assert simple_pipe.T_out.v == pytest.approx(simple_pipe.T_in.v)

    def test_pressure_output_mirrors_input(self, simple_pipe):
        assert simple_pipe.P_out.v == pytest.approx(simple_pipe.P_in.v)


# =============================================================================
# Zero / no-flow behaviour
# =============================================================================

class TestZeroFlow:
    """At zero mass flow the component must set DELTA_P = 0."""

    def test_delta_p_is_zero(self, simple_pipe):
        simple_pipe.m_dot_in.v = 0.0
        simple_pipe.calculate()
        assert simple_pipe.DELTA_P.v == pytest.approx(0.0)


# =============================================================================
# Mass and energy conservation
# =============================================================================

class TestConservation:
    """Basic conservation laws that must hold for incompressible pipe flow."""

    def test_mass_flow_conserved(self, simple_pipe):
        """Outlet mass flow equals inlet mass flow (no storage)."""
        simple_pipe.calculate()
        assert simple_pipe.m_dot_out.v == pytest.approx(simple_pipe.m_dot_in.v)

    def test_temperature_unchanged(self, simple_pipe):
        """Temperature is unchanged through an adiabatic incompressible pipe."""
        simple_pipe.calculate()
        assert simple_pipe.T_out.v == pytest.approx(simple_pipe.T_in.v)


# =============================================================================
# Pressure-drop behaviour
# =============================================================================

class TestPressureDrop:
    """Pressure-drop magnitude and direction for positive flow."""

    def test_pressure_drop_is_positive(self, simple_pipe):
        """DELTA_P must be positive for forward flow."""
        simple_pipe.calculate()
        assert simple_pipe.DELTA_P.v > 0.0

    def test_outlet_pressure_lower_than_inlet(self, simple_pipe):
        """Flow direction fixes that outlet pressure < inlet pressure."""
        simple_pipe.calculate()
        assert simple_pipe.P_out.v < simple_pipe.P_in.v

    def test_pressure_balance(self, simple_pipe):
        """P_out must equal P_in - DELTA_P."""
        simple_pipe.calculate()
        assert simple_pipe.P_out.v == pytest.approx(
            simple_pipe.P_in.v - simple_pipe.DELTA_P.v
        )

    @pytest.mark.parametrize("length, expect_larger", [
        (5.0,  False),   # shorter pipe → smaller drop than baseline (10 m)
        (20.0, True),    # longer  pipe → larger  drop
    ])
    def test_pressure_drop_scales_with_length(self, simple_pipe, length, expect_larger):
        """Pressure drop is proportional to pipe length (Darcy–Weisbach)."""
        # Baseline with the fixture's default length (10 m)
        simple_pipe.calculate()
        baseline_dp = simple_pipe.DELTA_P.v

        # Reconfigure length and recalculate
        simple_pipe.Length_Pipe.v = length
        simple_pipe.calculate()
        modified_dp = simple_pipe.DELTA_P.v

        if expect_larger:
            assert modified_dp > baseline_dp
        else:
            assert modified_dp < baseline_dp

    @pytest.mark.parametrize("m_dot", [0.5, 1.0, 2.0, 5.0])
    def test_pressure_drop_increases_with_flow_rate(self, simple_pipe, m_dot):
        """
        DELTA_P grows with mass flow rate (quadratic in the Darcy–Weisbach
        equation for turbulent flow).
        """
        simple_pipe.m_dot_in.v = m_dot
        simple_pipe.calculate()
        assert simple_pipe.DELTA_P.v > 0.0


# =============================================================================
# Friction factor helper
# =============================================================================

class TestFricFactor:
    """Standalone tests for the FricFactor_IC helper function."""

    def test_returns_float_turbulent(self):
        """Returns a finite positive float in the turbulent regime."""
        Re = 10_000
        ff = FricFactor_IC(Rough=4.6e-5 / 0.05, Reynold=Re, guess=0.03)
        assert ff is not None
        assert ff > 0.0

    def test_moody_chart_magnitude(self):
        """
        For smooth commercial steel at Re ≈ 10 000, friction factor should be
        in the range [0.015, 0.045] (standard Moody chart values).
        """
        Re = 10_000
        ff = FricFactor_IC(Rough=4.6e-5 / 0.05, Reynold=Re, guess=0.03)
        assert 0.015 < ff < 0.045
