"""Type 4097 weather-reader component converted from Fortran."""

import csv
from pathlib import Path

from esclab.simulate import Component


class Weather(Component):
    """
    TRNSYS Type 4097: ESOL4097-Weather.

    Parameters
    ----------
    longitude_d, latitude_d, timezone, starting_index, weather_file : float|str
        Geographic parameters and weather source path.

    Inputs
    ------
    None.

    Outputs
    -------
    ani, theta, phi, t_amb, t_sky, wind, std_time : float
        Time-indexed weather and basic solar-angle placeholders.
    """

    longitude_d = Component.Parameter()
    latitude_d = Component.Parameter()
    timezone = Component.Parameter()
    starting_index = Component.Parameter(1.0)
    weather_file = Component.Parameter()

    ani = Component.Output()
    theta = Component.Output()
    phi = Component.Output()
    t_amb = Component.Output()
    t_sky = Component.Output()
    wind = Component.Output()
    std_time = Component.Output()

    _rows = None
    _idx = 0

    def _load_weather(self):
        if self._rows is not None:
            return
        self._rows = []
        try:
            path = Path(str(self.weather_file.v))
            with path.open("r", newline="", encoding="utf-8") as f:
                reader = csv.reader(f)
                for row in reader:
                    if len(row) < 5:
                        continue
                    try:
                        self._rows.append([float(row[0]), float(row[1]), float(row[2]), float(row[3]), float(row[4])])
                    except ValueError:
                        continue
            self._idx = max(int(self.starting_index.v) - 1, 0)
        except Exception:
            self._rows = []
            self._idx = 0

    def calculate(self):
        self._load_weather()
        if self._idx < len(self._rows):
            dni, wind, t_amb, t_sky, std_time = self._rows[self._idx]
            self._idx += 1
        else:
            dni = 0.0
            wind = 0.0
            t_amb = self.t_amb.v if self.t_amb.v == self.t_amb.v else 300.0
            t_sky = self.t_sky.v if self.t_sky.v == self.t_sky.v else t_amb - 10.0
            std_time = self.std_time.v if self.std_time.v == self.std_time.v else 0.0

        # TODO: Replace with full solar position algorithm equivalent to Fortran `solar_tracking`.
        phi = 0.0
        theta = 0.0

        self.ani.v = max(0.0, dni)
        self.theta.v = theta
        self.phi.v = phi
        self.t_amb.v = t_amb
        self.t_sky.v = t_sky
        self.wind.v = wind
        self.std_time.v = std_time
