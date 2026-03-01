"""Parallel flow solver component model (Type 4050)."""

import numpy as np

from esclab.simulate import Component
from esclab.components.esol_properties import Incompressible as Inc


class ParallelFlowSolver(Component):
    """
    Object: ESOL4050-Parallel-Flow-Solver
    Simulation Studio Model: ESOL4050-Parallel-Flow-Solver

    Author: Matt Tuman
    Editor:
    Date:    February 05, 2024
    last modified: February 05, 2024
    Ported by: GitHub Copilot, March 01, 2026

    Solves parallel branch flow distribution for either a two-level nested
    parallel network (Solver=1, 8-branch system) or a simple single-level
    parallel network (Solver=2, 2-branch system).  Returns the fractional
    flow-split ratios f_23, f_45, and f_67.
    """

    # *** Model Parameters ***
    Solver = Component.Parameter()      # [-Inf;+Inf]
    Fluid_ID = Component.Parameter()   # [-Inf;+Inf]
    Coef_a = Component.Parameter()     # [-Inf;+Inf]
    Coef_b = Component.Parameter()     # [-Inf;+Inf]
    Coef_c = Component.Parameter()     # [-Inf;+Inf]
    N_pumps = Component.Parameter()    # [-Inf;+Inf]

    # *** Model Inputs ***
    m_dot_1 = Component.Input()    # [-Inf;+Inf]
    m_dot_2 = Component.Input()    # [-Inf;+Inf]
    m_dot_3 = Component.Input()    # [-Inf;+Inf]
    m_dot_4 = Component.Input()    # [-Inf;+Inf]
    m_dot_5 = Component.Input()    # [-Inf;+Inf]
    m_dot_6 = Component.Input()    # [-Inf;+Inf]
    m_dot_7 = Component.Input()    # [-Inf;+Inf]
    m_dot_8 = Component.Input()    # [-Inf;+Inf]
    m_dot_9 = Component.Input()    # [-Inf;+Inf]
    P_1 = Component.Input()        # [-Inf;+Inf]
    P_2 = Component.Input()        # [-Inf;+Inf]
    P_3 = Component.Input()        # [-Inf;+Inf]
    P_4 = Component.Input()        # [-Inf;+Inf]
    P_5 = Component.Input()        # [-Inf;+Inf]
    P_6 = Component.Input()        # [-Inf;+Inf]
    P_7 = Component.Input()        # [-Inf;+Inf]
    P_8 = Component.Input()        # [-Inf;+Inf]
    P_9 = Component.Input()        # [-Inf;+Inf]
    P_10 = Component.Input()       # [-Inf;+Inf]
    P_pump_in = Component.Input()  # [-Inf;+Inf]
    P_exp = Component.Input()      # [-Inf;+Inf]
    P_pump = Component.Input()     # [-Inf;+Inf]
    Pump_speed = Component.Input() # [-Inf;+Inf]
    T_1 = Component.Input()        # [-Inf;+Inf]

    # *** Model Outputs ***
    f_23 = Component.Output()  # [-Inf;+Inf]
    f_45 = Component.Output()  # [-Inf;+Inf]
    f_67 = Component.Output()  # [-Inf;+Inf]

    def presim_setup(self, **kwargs):
        super().presim_setup(**kwargs)

        # Allocate memory
        # If using solver for system 1
        if self.Solver.v == 1.0:
            self._A = np.zeros((8, 8))
            self._b = np.zeros(8)
            self._m_dots = np.zeros(8)

            # Enforce continuity in A
            self._A[0, 0] = 1.0
            self._A[0, 1] = 1.0

            self._A[1, 0] = -1.0
            self._A[1, 2] = 1.0
            self._A[1, 3] = 1.0

            self._A[2, 1] = -1.0
            self._A[2, 4] = 1.0
            self._A[2, 5] = 1.0

            self._A[3, 2] = 1.0
            self._A[3, 3] = 1.0
            self._A[3, 6] = -1.0

            self._A[4, 4] = 1.0
            self._A[4, 5] = 1.0
            self._A[4, 7] = -1.0

            self._b[:] = 0.0

        # If using solver for system 2
        elif self.Solver.v == 2.0:
            self._A = np.zeros((2, 2))
            self._b = np.zeros(2)
            self._m_dots = np.zeros(2)

            # Enforce continuity in A
            self._A[0, 0] = 1.0
            self._A[0, 1] = 1.0

            self._b[:] = 0.0

        # Set the Initial Values of the Outputs (#,Value)
        self.f_23.v = 0.5  # f_23
        self.f_45.v = 0.5  # f_45
        self.f_67.v = 0.5  # f_67

        return

    def calculate(self):
        super().calculate()

        # Solves parallel system below
        # If using solver for system 1
        if self.Solver.v == 1.0:
            #             | (1)
            #             |
            #             ^ (f)
            #            / \
            #      (2)  /   \
            #          /     \ (3)
            #         /       \
            #        |         |
            #        |         |
            #        ^         ^
            #       / \       / \
            #      /   \     /   \
            # (4) |     |   |     |
            #     |     |   |     |
            #     |    _|   |  (7)|
            #     |    |    |     |
            #     |    |    |     |
            #     | (5)|    |(6)  |
            #     |   _|    |     |
            #     |   |      |   |
            #      \ /        \ /
            #       U          U
            #       |          |
            #        \        /
            #     (8) \      / (9)
            #          \    /
            #           \  /
            #             U
            #             |
            #             |
            if self.model.timestep_iteration > 0 or self.model.time > self.model.timestep:
                g = 9.81  # gravity
                s = self.Pump_speed.v

                # Compute K values based on previous iteration
                K2 = abs((self.P_1.v - self.P_2.v) / (self.m_dot_2.v) ** 2)
                K3 = abs((self.P_1.v - self.P_3.v) / (self.m_dot_3.v) ** 2)
                K4 = abs((self.P_2.v - self.P_4.v) / (self.m_dot_4.v) ** 2)
                K5 = abs((self.P_2.v - self.P_5.v) / (self.m_dot_5.v) ** 2)
                K6 = abs((self.P_3.v - self.P_6.v) / (self.m_dot_6.v) ** 2)
                K7 = abs((self.P_3.v - self.P_7.v) / (self.m_dot_7.v) ** 2)
                K8 = abs((self.P_4.v - self.P_8.v) / (self.m_dot_8.v) ** 2)
                K9 = abs((self.P_7.v - self.P_9.v) / (self.m_dot_9.v) ** 2)

                # Update A matrix
                # Enforce continuity in A
                self._A[:, :] = 0.0

                self._A[0, 0] = 1.0
                self._A[0, 1] = 1.0

                self._A[1, 0] = -1.0
                self._A[1, 2] = 1.0
                self._A[1, 3] = 1.0

                self._A[2, 1] = -1.0
                self._A[2, 4] = 1.0
                self._A[2, 5] = 1.0

                self._A[3, 2] = 1.0
                self._A[3, 3] = 1.0
                self._A[3, 6] = -1.0

                self._A[4, 4] = 1.0
                self._A[4, 5] = 1.0
                self._A[4, 7] = -1.0

                # Enforce Pressure Drop Relations
                k2_hat = K2 * abs(self.m_dot_2.v)
                self._A[7, 0] = -k2_hat

                k3_hat = K3 * abs(self.m_dot_3.v)
                self._A[7, 1] = k3_hat

                k4_hat = K4 * abs(self.m_dot_4.v)
                self._A[5, 2] = -k4_hat
                self._A[7, 2] = -k4_hat

                self._A[5, 3] = K5 * abs(self.m_dot_5.v)

                k6_hat = K6 * abs(self.m_dot_6.v)
                self._A[6, 4] = -k6_hat
                self._A[7, 4] = k6_hat

                k7_hat = K7 * abs(self.m_dot_7.v)
                self._A[6, 5] = k7_hat

                self._A[7, 6] = -K8 * abs(self.m_dot_8.v)

                self._A[7, 7] = K9 * abs(self.m_dot_9.v)

                # Update b
                self._b[:] = 0.0
                self._b[0] = self.m_dot_1.v

                # Invert matrix and compute mass flow rates
                # TODO-NEEDS CONVERSION REVIEW: SergioScripts matrixinv replaced with numpy.linalg.solve
                # call matrixinv(A, A_inv, 8); m_dots = matmul(A_inv, b)
                self._m_dots = np.linalg.solve(self._A, self._b)

                # Update mass flow rates
                # LR = exp(-Real(getTimestepIteration())+1)
                LR = 0.25

                if self._m_dots[0] < 0.001:
                    m_hold = np.log10(abs(self._m_dots[0]) + 0.0001) * LR + np.log10(abs(self.m_dot_2.v) + 0.0001) * (1 - LR)
                    self.m_dot_2.v = 10 ** m_hold
                else:
                    self.m_dot_2.v = self._m_dots[0] * LR + self.m_dot_2.v * (1 - LR)

                if self._m_dots[1] < 0.001:
                    m_hold = np.log10(abs(self._m_dots[1]) + 0.0001) * LR + np.log10(abs(self.m_dot_3.v) + 0.0001) * (1 - LR)
                    self.m_dot_3.v = 10 ** m_hold
                else:
                    self.m_dot_3.v = self._m_dots[1] * LR + self.m_dot_3.v * (1 - LR)

                if self._m_dots[2] < 0.001:
                    m_hold = np.log10(abs(self._m_dots[2]) + 0.0001) * LR + np.log10(abs(self.m_dot_4.v) + 0.0001) * (1 - LR)
                    self.m_dot_4.v = 10 ** m_hold
                else:
                    self.m_dot_4.v = self._m_dots[2] * LR + self.m_dot_4.v * (1 - LR)

                if self._m_dots[3] < 0.001:
                    m_hold = np.log10(abs(self._m_dots[3]) + 0.0001) * LR + np.log10(abs(self.m_dot_5.v) + 0.0001) * (1 - LR)
                    self.m_dot_5.v = 10 ** m_hold
                else:
                    self.m_dot_5.v = self._m_dots[3] * LR + self.m_dot_5.v * (1 - LR)

                if self._m_dots[4] < 0.001:
                    m_hold = np.log10(abs(self._m_dots[4]) + 0.0001) * LR + np.log10(abs(self.m_dot_6.v) + 0.0001) * (1 - LR)
                    self.m_dot_6.v = 10 ** m_hold
                else:
                    self.m_dot_6.v = self._m_dots[4] * LR + self.m_dot_6.v * (1 - LR)

                if self._m_dots[5] < 0.001:
                    m_hold = np.log10(abs(self._m_dots[5]) + 0.0001) * LR + np.log10(abs(self.m_dot_7.v) + 0.0001) * (1 - LR)
                    self.m_dot_7.v = 10 ** m_hold
                else:
                    self.m_dot_7.v = self._m_dots[5] * LR + self.m_dot_7.v * (1 - LR)

                self.f_23.v = self.m_dot_2.v / (self.m_dot_2.v + self.m_dot_3.v)
                self.f_45.v = self.m_dot_4.v / (self.m_dot_4.v + self.m_dot_5.v)
                self.f_67.v = self.m_dot_6.v / (self.m_dot_6.v + self.m_dot_7.v)

            # Equally distribute flow
            else:
                self.f_23.v = 0.5
                self.f_45.v = 0.5
                self.f_67.v = 0.5
                return

        # Solves parallel system below
        elif self.Solver.v == 2.0:
            #             | (1)
            #             |
            #             ^ (f)
            #            / \
            #      (2)  /   \
            #          /     \ (3)
            #         /       \
            #        /         \
            #       |          |
            #       |          |
            #        \        /
            #         \      /
            #          \    /
            #           \  /
            #             U
            #             |
            #             |

            if self.model.timestep_iteration > 0 or self.model.time > self.model.timestep:

                # Compute K values based on previous iteration
                K2 = abs((self.P_1.v - self.P_2.v) / (self.m_dot_2.v) ** 2)
                K3 = abs((self.P_1.v - self.P_3.v) / (self.m_dot_3.v) ** 2)

                # Update A matrix
                # Enforce continuity in A
                self._A[:, :] = 0.0

                self._A[0, 0] = 1.0
                self._A[0, 1] = 1.0

                # Enforce Pressure Drop Relations
                k2_hat = K2 * abs(self.m_dot_2.v)
                self._A[1, 0] = -k2_hat

                k3_hat = K3 * abs(self.m_dot_3.v)
                self._A[1, 1] = k3_hat

                # Update b
                self._b[:] = 0.0
                self._b[0] = self.m_dot_1.v

                # Invert matrix and compute mass flow rates
                # TODO-NEEDS CONVERSION REVIEW: SergioScripts matrixinv replaced with numpy.linalg.solve
                # A2 = A(1:2,1:2); call matrixinv(A2, A2inv, 2); m_dots(1:2) = matmul(A2inv, b(1:2))
                self._m_dots = np.linalg.solve(self._A, self._b)

                # Update mass flow rates
                LR = 0.5

                # If flow through near field is shut off
                if self._m_dots[0] < 0.001:
                    m_hold = np.log10(self._m_dots[0]) * LR + np.log10(self.m_dot_2.v) * (1 - LR)
                    self.m_dot_2.v = 10 ** m_hold
                else:
                    self.m_dot_2.v = self._m_dots[0] * LR + self.m_dot_2.v * (1 - LR)
                # If flow through far field is shut off
                if self._m_dots[1] < 0.001:
                    m_hold = np.log10(self._m_dots[1]) * LR + np.log10(self.m_dot_3.v) * (1 - LR)
                    self.m_dot_3.v = 10 ** m_hold
                else:
                    self.m_dot_3.v = self._m_dots[1] * LR + self.m_dot_3.v * (1 - LR)

                self.f_23.v = self.m_dot_2.v / (self.m_dot_2.v + self.m_dot_3.v)
                self.f_45.v = 0.0
                self.f_67.v = 0.0

            # Equally distribute flow
            else:
                self.f_23.v = 0.5
                self.f_45.v = 0.0
                self.f_67.v = 0.0
                return

        # Set the Outputs from this Model (#,Value)
        # f_23, f_45, f_67 already written to self.f_23.v, self.f_45.v, self.f_67.v above

    def converged(self):
        # Perform Any "After Convergence" Manipulations That May Be Required at the End of Each Timestep
        # m_hold = gettimestepiteration()  (captured at end of timestep for diagnostics)
        m_hold = self.model.timestep_iteration
