"""Type 6034 TES mode selector converted from Fortran."""

from esclab.simulate import Component


class TESModes(Component):
    """
    TRNSYS Type 6034: TESModes.

    mode 0 = inactive, mode 1 = charging, mode 2 = discharging
    """

    cp_power = Component.Input()
    dp_power = Component.Input()
    htf_cvp_i = Component.Input()
    htf_cvp_a = Component.Input()
    htf_dvp_i = Component.Input()
    htf_dvp_a = Component.Input()

    mode = Component.Output()
    charging_pump_power = Component.Output()
    discharging_pump_power = Component.Output()
    htf_charging_valve_pos = Component.Output()
    htf_discharging_valve_pos = Component.Output()

    def calculate(self):
        mode = self.mode.v if self.mode.v == self.mode.v else 0.0
        cp = 1.0 if self.cp_power.v > 0.0 else 0.0
        dp = 1.0 if self.dp_power.v > 0.0 else 0.0
        cvp_i = max(self.htf_cvp_i.v, 0.0)
        dvp_i = max(self.htf_dvp_i.v, 0.0)
        cvp_a = max(self.htf_cvp_a.v, 0.0)
        dvp_a = max(self.htf_dvp_a.v, 0.0)

        if self.model.is_first_step:
            if cvp_i > 0.0 or cp > 0.0:
                mode = 1.0
                dp = 0.0
                dvp_i = 0.0
            elif dvp_i > 0.0 or dp > 0.0:
                mode = 2.0
                cp = 0.0
                cvp_i = 0.0
            else:
                mode = 0.0
        elif self.model.is_first_iteration:
            if mode == 1.0:
                if cp == 0.0 and cvp_a <= 0.0:
                    mode = 2.0 if (dp > 0.0 or dvp_i > 0.0) else 0.0
                else:
                    dp = 0.0
                    dvp_i = 0.0
            elif mode == 2.0:
                if dp == 0.0 and dvp_a <= 0.0:
                    mode = 1.0 if (cp > 0.0 or cvp_i > 0.0) else 0.0
                else:
                    cp = 0.0
                    cvp_i = 0.0
            else:
                if cp > 0.0 or cvp_i > 0.0:
                    mode = 1.0
                    dp = 0.0
                    dvp_i = 0.0
                elif dp > 0.0 or dvp_i > 0.0:
                    mode = 2.0
                    cp = 0.0
                    cvp_i = 0.0

        self.mode.v = mode
        self.charging_pump_power.v = cp
        self.discharging_pump_power.v = dp
        self.htf_charging_valve_pos.v = cvp_i
        self.htf_discharging_valve_pos.v = dvp_i
