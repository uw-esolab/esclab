"""
Solar position and tracking helper functions.

Converted from solar_field_modules.f90 :: Solar_Position module.
"""

import math
from typing import Tuple


def solar_tracking(
    TimeZone: float,
    LongD: float,
    LatD: float,
    julian_day: float,
    time: float,
) -> Tuple[float, float]:
    """Compute solar tracking angles for a parabolic trough with N/S axis.

    Computes the tracking angle (phi) and incidence angle (theta) for a
    horizontal collector rotating about the North/South axis based on
    Duffie & Beckman solar angle correlations.

    Parameters
    ----------
    TimeZone : float
        Time zone offset from UTC [degrees, e.g. -7 for MST]
    LongD : float
        Longitude of the site [degrees, positive East]
    LatD : float
        Latitude of the site [degrees, positive North]
    julian_day : float
        Julian day of the year [1-365]
    time : float
        Local standard time [hours, 0-24]

    Returns
    -------
    tuple of float
        (phi, theta) in radians:
        phi   - collector tracking angle [rad]
        theta - solar incidence angle on collector aperture [rad]
    """
    pi = 3.141592653
    n = julian_day
    L_st = 105.0  # Standard meridian [degrees]

    # 3.3 - New B per Duffie & Beckman 1.4.2
    B = (n - 1) * 360.0 / 365.0
    B_rad = B * pi / 180.0

    # 3.4 - Equation of Time in minutes
    EOT = 229.2 * (0.000075 + 0.001868 * math.cos(B_rad) - 0.032077 * math.sin(B_rad)
                   - 0.014615 * math.cos(B_rad * 2) - 0.04089 * math.sin(B_rad * 2.0))

    # 3.5 - Declination (per Duffie & Beckman 1.6.1a)
    dec = 23.45 * math.sin(360.0 * (284.0 + n) / 365.0 * pi / 180.0)
    dec_rad = dec * pi / 180.0

    # Solar Time in Hours
    SolarTime = time + (4.0 * (L_st - LongD) + EOT) / 60.0

    # 3.14 - Calculation of Hour Angle in radians
    HourAngle = (SolarTime - 12.0) * 15.0
    HourAngle_rad = HourAngle * pi / 180.0

    # 3.14 - Solar Altitude (radians)
    Lat_rad = LatD * pi / 180.0
    SolarAlt_rad = math.asin(
        math.sin(dec_rad) * math.sin(Lat_rad)
        + math.cos(Lat_rad) * math.cos(dec_rad) * math.cos(HourAngle_rad)
    )

    # 3.15 - Solar azimuth
    cos_arg = min(
        1.0,
        (math.cos(pi / 2.0 - SolarAlt_rad) * math.sin(Lat_rad) - math.sin(dec_rad))
        / (math.sin(pi / 2.0 - SolarAlt_rad) * math.cos(Lat_rad)),
    )
    SolarAz_rad = math.copysign(1.0, HourAngle_rad) * abs(math.acos(cos_arg))

    # 3.16 - Solar Zenith
    SolarZenith_rad = pi / 2.0 - SolarAlt_rad

    # 3.17 - Collector azimuth (assumes N/S horizontal axis, collector tracks E/W)
    if HourAngle_rad > 0:
        ColAz_rad = 90.0 * pi / 180.0
    else:
        ColAz_rad = -90.0 * pi / 180.0

    # Compute tracking angle (phi)
    phi = math.atan(math.tan(SolarZenith_rad) * math.sin(SolarAz_rad))

    # Compute cos(theta)
    if phi < 0:
        CosTh = 1.0
    else:
        CosTh = (
            -math.sin(SolarZenith_rad) * math.sin(SolarAz_rad) * (-math.sin(phi))
            + math.cos(SolarZenith_rad) * math.cos(phi)
        )

    theta = math.acos(CosTh)

    return phi, theta
