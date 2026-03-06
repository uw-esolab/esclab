"""Weather reader component model (Type 4097)."""

import numpy as np

from esclab.simulate import Component
from esclab.components.flownetwork.solar_position import solar_tracking


class WeatherReader(Component):
    """
    Object: 4097-Weather
    Simulation Studio Model: ESOL4097-Weather

    Author: Matt Tuman
    Editor:
    Date:     January 02, 2024
    last modified: January 02, 2024

    longitudeD  - [-Inf;+Inf]
    latitudeD   - [-Inf;+Inf]
    timezone    - [-Inf;+Inf]

    Model Inputs: (none)

    Model Outputs:
        DNI     - [-Inf;+Inf]
        Theta   - [-Inf;+Inf]
        Phi     - [-Inf;+Inf]
        T_amb   - [-Inf;+Inf]
        T_sky   - [-Inf;+Inf]
        Wind    - [-Inf;+Inf]
        stdTime - [-Inf;+Inf]
    """
    trnsys_type = "4097"

    # PARAMETERS
    longitudeD = Component.Parameter()      # Longitude [deg], Duffie sign convention (west is +, compared to Google Maps where west is -)
    latitudeD = Component.Parameter()       # Latitude [deg]
    timezone = Component.Parameter()        # Time zone
    starting_index = Component.Parameter()  # Starting index into the weather file

    # OUTPUTS
    ANI = Component.Output()        # 1 - Absorbed normal irradiance: max(0, DNI*cos(theta))
    theta = Component.Output()      # 2 - Solar zenith/incidence angle [rad]
    phi = Component.Output()        # 3 - Solar azimuth angle [rad]
    T_amb = Component.Output()      # 4 - Ambient temperature
    T_sky = Component.Output()      # 5 - Sky temperature
    Wind = Component.Output()       # 6 - Wind speed
    stdTime = Component.Output()    # 7 - Standard time

    def calculate(self):
        super().calculate()

        # Get the Global Trnsys Simulation Variables
        Timestep = self.model.timestep * 3600.0  # [s]

        # Do All of the First Timestep Manipulations Here - There Are No Iterations at the Initial Time
        if self.model.is_first_step:

            # Read in parameters
            # longitudeD = getParameterValue(1)  -- accessed via self.longitudeD.v
            # latitudeD  = getParameterValue(2)  -- accessed via self.latitudeD.v
            # timezone   = getParameterValue(3)  -- accessed via self.timezone.v
            # starting_index = getParameterValue(4) -- accessed via self.starting_index.v

            # TODO-NEEDS CONVERSION REVIEW: how getLabel(CurrentUnit, 1) maps in esclab
            Filename = self.model.label

            self._xyzabc = int(self.starting_index.v)

            # Open weather file
            self._file_handle = open(Filename, 'r')
            next(self._file_handle)                         # READ(22, *)  -- skip first header line
            line = next(self._file_handle)                  # Read(22, *), day, month, n_time
            tokens = line.split()
            day = int(tokens[0])
            month = int(tokens[1])
            n_time = float(tokens[2])
            next(self._file_handle)                         # READ(22, *)  -- skip third header line

            # Skip to initial time
            for n in range(1, self._xyzabc):
                next(self._file_handle)                     # unread weather columns

            # Compute Julian Day
            if month == 1:
                self._julian_day = float(day)
            elif month == 2:
                self._julian_day = float(day + 31)
            elif month == 3:
                self._julian_day = float(day + 59)
            elif month == 4:
                self._julian_day = float(day + 90)
            elif month == 5:
                self._julian_day = float(day + 120)
            elif month == 6:
                self._julian_day = float(day + 151)
            elif month == 7:
                self._julian_day = float(day + 181)
            elif month == 8:
                self._julian_day = float(day + 212)
            elif month == 9:
                self._julian_day = float(day + 243)
            elif month == 10:
                self._julian_day = float(day + 273)
            elif month == 11:
                self._julian_day = float(day + 304)
            else:
                self._julian_day = float(day + 334)

            # Check that the weather file is at least longer than the simulation
            # TODO-NEEDS CONVERSION REVIEW: getSimulationStopTime() mapped to self.model.settings.stop_time - verify attribute name
            len_sim = self.model.settings.stop_time * 3600.0  # [s]
            n_time_sim = len_sim / Timestep
            if n_time_sim > n_time:
                # Call FoundBadParameter(1, 'Fatal', 'Number of simulation timesteps is greater than the weather file')
                raise ValueError(
                    "Number of simulation timesteps is greater than the weather file"
                )

            # Initialize storage arrays
            self._ANI_store = np.zeros(30000)
            self._stdTime_store = np.zeros(30000)

            # Initialize local state variables
            self._DNI = 0.0
            self._Wind = 0.0
            self._T_amb = 0.0
            self._T_sky = 0.0
            self._stdTime = 0.0
            self._phi = 0.0
            self._theta = 0.0

            # Set the Initial Values of the Outputs (#,Value)
            self.ANI.v = 0.0        # 1 DNI
            self.theta.v = 0.0      # 2 Theta
            self.phi.v = 0.0        # 3 Phi
            self.T_amb.v = 0.0      # 4 T_amb
            self.T_sky.v = 0.0      # 5 T_sky
            self.Wind.v = 0.0       # 6 Wind
            self.stdTime.v = 0.0    # 7 stdTime

            return

        # Only evaluate during the first iteration
        if self.model.timestep_iteration == 0:
            # Read in weather/time data
            line = next(self._file_handle)
            tokens = line.split()
            self._DNI = float(tokens[0])
            self._Wind = float(tokens[1])
            self._T_amb = float(tokens[2])
            self._T_sky = float(tokens[3])
            self._stdTime = float(tokens[4])

            # Compute solar angles
            # longitudeD uses Duffie sign convention (west is +, compared to Google Maps, where west is -)
            angle_results = solar_tracking(
                self.timezone.v,
                self.longitudeD.v,
                self.latitudeD.v,
                self._julian_day,
                self._stdTime,
            )
            self._phi = angle_results[0]
            self._theta = angle_results[1]

            self._ANI_store[self._xyzabc - 1] = max(0.0, self._DNI * np.cos(self._theta))
            self._stdTime_store[self._xyzabc - 1] = self._stdTime

            self._xyzabc = self._xyzabc + 1

        # Set the Outputs from this Model (#,Value)
        self.ANI.v = max(0.0, self._DNI * np.cos(self._theta))     # 1 ANI
        self.theta.v = self._theta                                   # 2 Theta
        self.phi.v = self._phi                                       # 3 Phi
        self.T_amb.v = self._T_amb                                   # 4 T_amb
        self.T_sky.v = self._T_sky                                   # 5 T_sky
        self.Wind.v = self._Wind                                     # 6 Wind
        self.stdTime.v = self._stdTime                               # 7 stdTime

    def converged(self):
        """Perform Any "After Convergence" Manipulations That May Be Required at the End of Each Timestep."""
        # open(1, file = 'find_66.txt')
        # write(1, '(F9.6, A, I5.1)') getSimulationTime(), " end of timestep ", CurrentUnit
        return

    def finalize(self):
        """Do Any Last Call Manipulations Here."""
        with open('ANI.txt', 'w') as f:
            for n in range(self._xyzabc - 1):
                f.write(f"{self._ANI_store[n]}\n")

        with open('stdTime.txt', 'w') as f:
            for n in range(self._xyzabc - 1):
                f.write(f"{self._stdTime_store[n]}\n")

        self._file_handle.close()
