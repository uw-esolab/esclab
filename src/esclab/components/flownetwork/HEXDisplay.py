"""Type 4102 HEX display mapper converted from Fortran."""

from esclab.simulate import Component


class HEXDisplay(Component):
    """
    TRNSYS Type 4102: ESOL4102-HEX-Display.

    Parameters
    ----------
    None.

    Inputs
    ------
    input_1..input_25 : float
        TES mode plus charging/discharging shell/tube in/out streams.

    Outputs
    -------
    output_1..output_16 : float
        Top/bottom display streams and finite-difference temperature derivatives.
    """

    for _idx in range(1, 26):
        locals()[f"input_{_idx}"] = Component.Input()
    for _idx in range(1, 17):
        locals()[f"output_{_idx}"] = Component.Output()

    _prev = [0.0, 0.0, 0.0, 0.0]

    def calculate(self):
        mode = 1.0 if self.input_1.v < 2.0 else 2.0
        if mode == 1.0:
            vals = [
                self.input_2.v, self.input_3.v, self.input_4.v,
                self.input_20.v, self.input_21.v, self.input_22.v,
                self.input_14.v, self.input_15.v, self.input_16.v,
                self.input_8.v, self.input_9.v, self.input_10.v,
            ]
        else:
            vals = [
                self.input_17.v, self.input_18.v, self.input_19.v,
                self.input_11.v, self.input_12.v, self.input_13.v,
                self.input_5.v, self.input_6.v, self.input_7.v,
                self.input_23.v, self.input_24.v, self.input_25.v,
            ]

        for idx, val in enumerate(vals, start=1):
            getattr(self, f"output_{idx}").v = val

        ts = max(self.model.settings.timestep, 1.0)
        curr = [self.output_2.v, self.output_5.v, self.output_8.v, self.output_11.v]
        for i in range(4):
            getattr(self, f"output_{13+i}").v = (curr[i] - self._prev[i]) / ts
        self._prev = curr
