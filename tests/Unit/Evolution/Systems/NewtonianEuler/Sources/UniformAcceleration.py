# Distributed under the MIT License.
# See LICENSE.txt for details.

import numpy as np


def momentum_density_source(mass_density, momentum_density, acceleration_field):
    return mass_density * np.array(acceleration_field)


def energy_density_source(mass_density, momentum_density, acceleration_field):
    return np.dot(momentum_density, np.array(acceleration_field))


