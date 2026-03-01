"""Simple tee return mixing component model (Type 4016)."""

import numpy as np
from eeslib import fluid_properties as fp
from esclab.components.esol_properties import Incompressible as Inc

from esclab.simulate import Component


class TeeReturn(Component):
    """
    TRNSYS Type 4016: ESOL4016-TeeReturn-Simple.

    Combines two return branches into one stream with mass-weighted energy
    balance and selects outlet pressure based on the pressure-through
    parameter.

    Parameters
    ----------
    Pressure_Through : float
        Pressure selection flag. `0` uses `P_1`; otherwise uses `P_2`.
    Fluid_ID : str
        Fluid identifier used by `eeslib.fluid_properties`.

    Inputs
    ------
    m_dot_1 : float
        Branch 1 mass flow rate [kg/s].
    T_1 : float
        Branch 1 temperature [K].
    P_1 : float
        Branch 1 pressure [Pa].
    m_dot_2 : float
        Branch 2 mass flow rate [kg/s].
    T_2 : float
        Branch 2 temperature [K].
    P_2 : float
        Branch 2 pressure [Pa].
    Mass_Counter_1 : float
        Branch 1 mass counter [kg].
    Mass_Counter_2 : float
        Branch 2 mass counter [kg].

    Outputs
    -------
    m_dot : float
        Outlet mass flow rate [kg/s].
    Pressure : float
        Outlet pressure [Pa].
    Temperature : float
        Outlet temperature [K].
    P_1_in : float
        Echo of branch 1 inlet pressure [Pa].
    P_2_in : float
        Echo of branch 2 inlet pressure [Pa].
    Mass_Counter : float
        Combined mass counter [kg].
    """

    # *** Model Parameters ***
    # Pressure_Through: selects which branch pressure passes through to outlet
    Pressure_Through = Component.Parameter()
    # Fluid_ID: fluid identifier for specific heat property lookups
    Fluid_ID = Component.Parameter()

    # *** Model Inputs ***
    m_dot_1 = Component.Input()
    T_1 = Component.Input()
    P_1 = Component.Input()
    m_dot_2 = Component.Input()
    T_2 = Component.Input()
    P_2 = Component.Input()
    Mass_Counter_1 = Component.Input()
    Mass_Counter_2 = Component.Input()

    # *** Model Outputs ***
    m_dot = Component.Output()
    Pressure = Component.Output()
    Temperature = Component.Output()
    P_1_in = Component.Output()
    P_2_in = Component.Output()
    Mass_Counter = Component.Output()
    _props = Inc()

    def _select_pressure_out(self):
        if self.Pressure_Through.v == 0.0:
            return self.P_1.v
        return self.P_2.v

    def _specific_heat(self, temperature, pressure):
        # Fallback temperature when T was not reasonable.
        t_eval = temperature if temperature > 273.0 else 300.0
        fluid_name = str(self.Fluid_ID.v) if self.Fluid_ID.v == self.Fluid_ID.v else "Nitrate Salt"
        try:
            return float(self._props.specheat(fluid_name, t_eval, pressure))
        except Exception:
            pass
        try:
            return float(fp.specheat(fluid_name, T=t_eval, P=pressure))
        except Exception:
            # Conservative fallback to keep simulation running if Fluid_ID is not recognized.
            return float(fp.specheat("Water", T=t_eval, P=pressure))

    def _solve_mixed_temperature(self):
        # -------------------------------------------------------------------------------------------------------
        # Mass Balance
        # -------------------------------------------------------------------------------------------------------
        m_dot_tot = self.m_dot_1.v + self.m_dot_2.v
        if np.abs(m_dot_tot) < 1.0e-12:
            if np.isfinite(self.Temperature.v):
                return float(self.Temperature.v)
            return float(np.clip(0.5 * (self.T_1.v + self.T_2.v), 273.0, 1000.0))

        # -------------------------------------------------------------------------------------------------------
        # Energy Balance
        # Compute c_p(T1), c_p(T2), then solve for T_out such that c_p(T_out)*T_out
        # matches the mixed enthalpy term from both branches.
        # -------------------------------------------------------------------------------------------------------
        cp_t1 = self._specific_heat(self.T_1.v, self.P_1.v)
        cp_t2 = self._specific_heat(self.T_2.v, self.P_2.v)
        cp_times_t_out = (
            self.m_dot_1.v * cp_t1 * self.T_1.v
            + self.m_dot_2.v * cp_t2 * self.T_2.v
        ) / m_dot_tot

        if np.isfinite(self.Temperature.v):
            t_out_guess = float(np.clip(self.Temperature.v, 273.0, 1000.0))
        else:
            t_out_guess = float(np.clip(0.5 * (self.T_1.v + self.T_2.v), 273.0, 1000.0))

        tol = 100.0
        learning_rate = 0.5
        upperbound = 1000.0
        lowerbound = 273.0
        error = 1000.0
        error_prev = float("nan")
        t_out_guess_prev = float("nan")
        n_iter = 0

        # Iterative root solve with secant-like update and bounded guesses.
        while np.abs(error) > tol and n_iter < 1000:
            n_iter += 1
            cp_guess = self._specific_heat(t_out_guess, self.P_1.v)
            cp_times_t_out_guess = t_out_guess * cp_guess
            error = cp_times_t_out - cp_times_t_out_guess

            if np.abs(error) < tol:
                break

            if n_iter >= 2 and t_out_guess != t_out_guess_prev:
                slope = (error_prev - error) / (t_out_guess_prev - t_out_guess)
                intercept = error - slope * t_out_guess
                if slope != 0.0:
                    t_out_new = np.clip(-intercept / slope, lowerbound, upperbound)
                    t_out_guess_prev = t_out_guess
                    error_prev = error
                    # Relaxed update to improve stability.
                    t_out_guess = t_out_guess + (t_out_new - t_out_guess) * learning_rate
                else:
                    t_out_guess_prev = t_out_guess
                    error_prev = error
                    if error > 0.0:
                        t_out_guess = min(t_out_guess + 1.0, upperbound)
                    else:
                        t_out_guess = max(t_out_guess - 1.0, lowerbound)
            else:
                t_out_guess_prev = t_out_guess
                error_prev = error
                if error > 0.0:
                    t_out_guess = min(t_out_guess + 1.0, upperbound)
                else:
                    t_out_guess = max(t_out_guess - 1.0, lowerbound)

        return float(t_out_guess)

    def calculate(self):
        # -------------------------------------------------------------------------------------------------------
        # Common calculations for every call
        # -------------------------------------------------------------------------------------------------------
        m_dot_out = self.m_dot_1.v + self.m_dot_2.v
        mass_counter_out = self.Mass_Counter_1.v + self.Mass_Counter_2.v
        p_out = self._select_pressure_out()

        if self.model.is_first_step:
            # Equivalent to the Fortran getIsStartTime() block.
            # Keep the original initial outlet temperature value.
            self.m_dot.v = m_dot_out
            self.Pressure.v = p_out
            self.Temperature.v = 500.0
            self.P_1_in.v = self.P_1.v
            self.P_2_in.v = self.P_2.v
            self.Mass_Counter.v = mass_counter_out
            return

        # Main iteration calculations after first timestep.
        t_out = self._solve_mixed_temperature()

        # Set outputs from this model.
        self.m_dot.v = m_dot_out
        self.Pressure.v = p_out
        self.Temperature.v = t_out
        self.P_1_in.v = self.P_1.v
        self.P_2_in.v = self.P_2.v
        self.Mass_Counter.v = mass_counter_out


# Test code to run when executing this file directly
if __name__ == "__main__":
    class _DemoModel:
        is_first_step = True

    comp = TeeReturn()
    comp.model = _DemoModel()

    comp.Pressure_Through.v = 0.0
    comp.Fluid_ID.v = "Water"

    comp.m_dot_1.v = 3.0
    comp.T_1.v = 560.0
    comp.P_1.v = 220000.0
    comp.m_dot_2.v = 2.0
    comp.T_2.v = 520.0
    comp.P_2.v = 210000.0
    comp.Mass_Counter_1.v = 100.0
    comp.Mass_Counter_2.v = 80.0

    comp.calculate()
    print("First-step outputs:")
    print(f"  m_dot={comp.m_dot.v:.3f} kg/s")
    print(f"  Pressure={comp.Pressure.v:.1f} Pa")
    print(f"  Temperature={comp.Temperature.v:.3f} K")
    print(f"  Mass_Counter={comp.Mass_Counter.v:.3f} kg")

    comp.model.is_first_step = False
    comp.calculate()
    print("\nIteration outputs:")
    print(f"  m_dot={comp.m_dot.v:.3f} kg/s")
    print(f"  Pressure={comp.Pressure.v:.1f} Pa")
    print(f"  Temperature={comp.Temperature.v:.3f} K")
    print(f"  Mass_Counter={comp.Mass_Counter.v:.3f} kg")
