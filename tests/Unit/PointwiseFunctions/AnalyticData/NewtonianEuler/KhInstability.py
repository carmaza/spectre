# Distributed under the MIT License.
# See LICENSE.txt for details.

import numpy as np


def mass_density(x, adiabatic_index, strip_bimedian, strip_thickness,
                 strip_density, strip_velocity, background_density,
                 background_velocity, pressure, perturbation_amplitude,
                 damping_factor):
    return strip_density \
        if np.absolute(x[x.size - 1] - strip_bimedian) < 0.5 * strip_thickness \
           else background_density


def velocity(x, adiabatic_index, strip_bimedian, strip_thickness,
             strip_density, strip_velocity, background_density,
             background_velocity, pressure, perturbation_amplitude,
             damping_factor):
    dim = x.size
    result = np.zeros(dim)
    result[0] = strip_velocity \
                if (np.absolute(x[x.size - 1] - strip_bimedian) \
                    < 0.5 * strip_thickness) \
                    else background_velocity
    strip_lower_bound = strip_bimedian - 0.5 * strip_thickness
    strip_upper_bound = strip_bimedian + 0.5 * strip_thickness
    result[dim - 1] = (
        np.exp(-0.5 * ((x[dim - 1] - strip_lower_bound) / damping_factor)**2) +
        np.exp(-0.5 * ((x[dim - 1] - strip_upper_bound) / damping_factor)**2))
    result[dim - 1] *= perturbation_amplitude * np.sin(4 * np.pi * x[0])
    return result


def specific_internal_energy(x, adiabatic_index, strip_bimedian,
                             strip_thickness, strip_density, strip_velocity,
                             background_density, background_velocity, pressure,
                             perturbation_amplitude, damping_factor):
    return pressure / (adiabatic_index - 1.0) / strip_density \
        if np.absolute(x[x.size - 1] - strip_bimedian) < 0.5 * strip_thickness \
           else pressure / (adiabatic_index - 1.0) / background_density


def pressure(x, adiabatic_index, strip_bimedian, strip_thickness,
             strip_density, strip_velocity, background_density,
             background_velocity, pressure, perturbation_amplitude,
             damping_factor):
    return pressure
