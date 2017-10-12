import logging
from math import (pi,
                  exp,
                  log)
from random import random

import numpy as np
import sys

from alcor.models.star import GalacticStructureType

logger = logging.getLogger(__name__)


# TODO: take these consts to settings file
def generate_stars(max_stars_count: int = 6_000_000,
                   time_bins_count: int = 5_000,
                   thin_disk_scale_height_kpc: float = 0.25,
                   thick_disk_scale_height_kpc: float = 0.9,
                   thick_disk_age: float = 12.,
                   thick_disk_sfr_param: float = 2.,
                   burst_formation_factor: float = 5.,
                   star_formation_rate_param: float = 25.,
                   sector_radius_kpc: float = 0.05,
                   thin_disk_age: float = 9.2,
                   burst_age: float = 0.6,
                   mass_reduction_factor: float = 0.03,
                   thick_disk_stars_fraction: float = 0.8,
                   initial_mass_function_param: float = -2.35,
                   halo_stars_formation_time: float = 1.,
                   halo_stars_fraction: float = 0.05,
                   halo_age: float = 14.) -> None:
    time_increment = thin_disk_age / time_bins_count
    sector_area = pi * sector_radius_kpc ** 2
    birth_rate = (time_increment * sector_area * 1E6  # TODO: what is 1E6?
                  * mass_reduction_factor
                  * normalization_const(
                        star_formation_rate_param=star_formation_rate_param,
                        thin_disk_age_gyr=thin_disk_age))
    burst_birth_rate = birth_rate * burst_formation_factor
    max_age = max(thin_disk_age, thick_disk_age, halo_age)
    burst_init_time = max_age - burst_age
    thin_disk_birth_init_time = max_age - thin_disk_age
    thick_disk_birth_init_time = max_age - thick_disk_age
    halo_birth_init_time = max_age - halo_age
    # This can be proved by taking derivative from y = t * exp(-t / tau)
    thick_disk_max_sfr_relative_time = thick_disk_birth_init_time
    thick_disk_max_sfr = (thick_disk_max_sfr_relative_time
                          * exp(-thick_disk_max_sfr_relative_time
                                / thick_disk_sfr_param))

    stars_count = 0

    progenitors_masses = []
    galactic_structure_types = []
    birth_times = []
    z_coordinates = []

    for time_bin in range(time_bins_count):
        total_bin_mass = 0.

        current_bin_init_time = (thin_disk_birth_init_time
                                 + time_bin * time_increment)

        # TODO: do smth about this
        if current_bin_init_time >= burst_init_time:
            birth_rate = burst_birth_rate

        while True:
            stars_count += 1

            if stars_count > max_stars_count:
                raise Exception('Number of stars is too high - '
                                'decrease mass reduction factor')

            star_mass = get_mass_from_salpeter_initial_mass_function(
                    initial_mass_function_param=initial_mass_function_param)
            progenitors_masses.append(star_mass)

            if random() <= thick_disk_stars_fraction:
                galactic_structure_types.append(GalacticStructureType.thick)

                # TODO: implement MonteCarlo method function/Smirnov transform
                while True:
                    time_try = thick_disk_age * random()
                    time_try_sfr = (time_try
                                    * exp(-time_try
                                          / thick_disk_sfr_param))
                    sfr_try = thick_disk_max_sfr * random()
                    if sfr_try <= time_try_sfr:
                        birth_times.append(time_try
                                           + thick_disk_birth_init_time)
                        break
            elif (thick_disk_stars_fraction < random()
                  <= thick_disk_stars_fraction + halo_stars_fraction):
                galactic_structure_types.append(GalacticStructureType.halo)
                birth_times.append(halo_birth_init_time
                                   + halo_stars_formation_time * random())
            else:
                galactic_structure_types.append(GalacticStructureType.thin)
                birth_times.append(thin_disk_birth_init_time
                                   + time_bin * time_increment
                                   + time_increment * random())
                total_bin_mass += star_mass

            # TODO: assigning coords should be done in another module
            if galactic_structure_types[-1] == GalacticStructureType.thin:
                scale_height = thin_disk_scale_height_kpc
            elif galactic_structure_types[-1] == GalacticStructureType.thick:
                scale_height = thick_disk_scale_height_kpc

            # Halo stars coords will be generated in polar coords module
            if galactic_structure_types[-1] != GalacticStructureType.halo:
                # Inverse transform sampling for y = exp(-z / H)
                z_coordinate = (-scale_height * log(
                        1. - random() * (1.0 - exp(-sector_radius_kpc
                                                   / scale_height))))
                # TODO: change this
                random_sign = float(1. - 2. * int(2.0 * random()))
                z_coordinates.append(z_coordinate * random_sign)

            # TODO: find out what is going on here
            if total_bin_mass > birth_rate:
                progenitors_masses[-1] -= total_bin_mass - birth_rate
                break


def normalization_const(star_formation_rate_param: float,
                        thin_disk_age_gyr: float,
                        sigma: float = 51.  # TODO: what is sigma?
                        ) -> float:
    return sigma / (star_formation_rate_param
                    * (1 - np.exp(-thin_disk_age_gyr
                                  / star_formation_rate_param)))


# TODO: implement inverse transform sampling
def get_mass_from_salpeter_initial_mass_function(
        initial_mass_function_param: float,
        min_mass: float = 0.4,
        max_mass: float = 50.
        ) -> float:
    y_max = min_mass ** initial_mass_function_param

    mass_amplitude = max_mass - min_mass
    while True:
        # TODO: implement seeds tracking
        y = y_max * np.random.rand()
        mass = min_mass + mass_amplitude * np.random.rand()[0]
        y_imf = mass ** initial_mass_function_param
        if y <= y_imf:
            return mass