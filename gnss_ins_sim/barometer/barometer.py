# -*- coding: utf-8 -*-
# Filename: barometer.py



# import
# import math
# import numpy as np


def altitude_to_pressure_vec(altitudes_m):
    """
    Converts an array of altitudes (in meters) to atmospheric pressures (in Pascals).

    Parameters:
    - altitudes_m: NumPy array of altitudes in meters

    Returns:
    - NumPy array of pressures in Pascals
    """
    P0 = 101325  # Pressure at sea level (Pa)
    T0 = 288.15  # Temperature at sea level (K)
    g = 9.80665  # Gravitational acceleration (m/s^2)
    L = 0.0065  # Temperature lapse rate (K/m)
    M = 0.0289644  # Molar mass of Earth's air (kg/mol)
    R = 8.3144598  # Universal gas constant (J/(mol·K))

    return P0 * (1 - (L * altitudes_m) / T0) ** ((g * M) / (R * L))


def pressure_to_altitude_vec(pressures_Pa):
    """
    Converts an array of atmospheric pressures (in Pascals) to altitudes (in meters).

    Parameters:
    - pressures_Pa: NumPy array of pressures in Pascals

    Returns:
    - NumPy array of altitudes in meters
    """
    P0 = 101325  # Pressure at sea level (Pa)
    T0 = 288.15  # Temperature at sea level (K)
    g = 9.80665  # Gravitational acceleration (m/s^2)
    L = 0.0065  # Temperature lapse rate (K/m)
    M = 0.0289644  # Molar mass of Earth's air (kg/mol)
    R = 8.3144598  # Universal gas constant (J/(mol·K))

    return (T0 / L) * (1 - (pressures_Pa / P0) ** ((R * L) / (g * M)))

