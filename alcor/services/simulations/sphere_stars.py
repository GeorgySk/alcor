import logging
import math
from typing import (Union,
                    Dict,
                    List)
from random import random

from alcor.models.star import GalacticStructureType

logger = logging.getLogger(__name__)


# TODO: take these consts to settings file
def generate_stars(*,
                   max_stars_count: int = 6000000,
                   time_bins_count: int = 5000,
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
                   halo_age: float = 14.
                   ) -> Dict[str, List[Union[float, GalacticStructureType]]]:
    time_increment = thin_disk_age / time_bins_count
    sector_area = math.pi * sector_radius_kpc ** 2
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
                          * math.exp(-thick_disk_max_sfr_relative_time
                                     / thick_disk_sfr_param))

    stars_count = 0

    progenitors_masses = []
    galactic_structure_types = []
    birth_times = []

    for time_bin in range(time_bins_count):
        total_bin_mass = 0.

        current_bin_init_time = (thin_disk_birth_init_time
                                 + time_bin * time_increment)

        # TODO: implement birth rate function
        if current_bin_init_time >= burst_init_time:
            birth_rate = burst_birth_rate

        while True:
            stars_count += 1

            if stars_count > max_stars_count:
                raise OverflowError('Number of stars is too high - '
                                    'decrease mass reduction factor.')

            star_mass = salpeter_initial_mass_function(
                    parameter=initial_mass_function_param)
            progenitors_masses.append(star_mass)

            galactic_structure_type = get_galactic_structure_type(
                    thick_disk_stars_fraction=thick_disk_stars_fraction,
                    halo_stars_fraction=halo_stars_fraction)
            galactic_structure_types.append(galactic_structure_type)

            if galactic_structure_type == GalacticStructureType.thick:
                # TODO: implement MonteCarlo method function/Smirnov transform
                birth_times.append(thick_disk_star_birth_time(
                        thick_disk_age=thick_disk_age,
                        thick_disk_sfr_param=thick_disk_sfr_param,
                        thick_disk_max_sfr=thick_disk_max_sfr,
                        thick_disk_birth_init_time=thick_disk_birth_init_time))
            elif galactic_structure_type == GalacticStructureType.halo:
                birth_times.append(halo_star_birth_time(
                        halo_birth_init_time=halo_birth_init_time,
                        halo_stars_formation_time=halo_stars_formation_time))
            else:
                birth_times.append(thin_disk_star_birth_time(
                        thin_disk_birth_init_time=thin_disk_birth_init_time,
                        time_bin=time_bin,
                        time_increment=time_increment))

            if galactic_structure_type != GalacticStructureType.thin:
                continue

            total_bin_mass += star_mass

            if total_bin_mass > birth_rate:
                break

    return dict(progenitors_masses=progenitors_masses,
                galactic_structure_types=galactic_structure_types,
                birth_times=birth_times)


def normalization_const(*,
                        star_formation_rate_param: float,
                        thin_disk_age_gyr: float,
                        sigma: float = 51.  # TODO: what is sigma?
                        ) -> float:
    return sigma / (star_formation_rate_param
                    * (1 - math.exp(-thin_disk_age_gyr
                                    / star_formation_rate_param)))


# TODO: implement inverse transform sampling
def salpeter_initial_mass_function(*,
                                   parameter: float,
                                   min_mass: float = 0.4,
                                   max_mass: float = 50.) -> float:
    y_max = min_mass ** parameter

    mass_amplitude = max_mass - min_mass
    while True:
        # TODO: implement seeds tracking
        y = y_max * random()
        mass = min_mass + mass_amplitude * random()
        y_imf = mass ** parameter
        if y <= y_imf:
            return mass


def thick_disk_star_birth_time(*,
                               thick_disk_age: float,
                               thick_disk_sfr_param: float,
                               thick_disk_max_sfr: float,
                               thick_disk_birth_init_time: float) -> float:
    """
    Return birth time of a thick disk star by using Monte Carlo method.
    SFR - star formation rate. More info at:
    https://www.google.es/search?q=star+formation+rate
    """
    while True:
        time_try = thick_disk_age * random()
        time_try_sfr = time_try * math.exp(-time_try / thick_disk_sfr_param)
        sfr_try = thick_disk_max_sfr * random()
        if sfr_try <= time_try_sfr:
            return time_try + thick_disk_birth_init_time


def halo_star_birth_time(*,
                         halo_birth_init_time: float,
                         halo_stars_formation_time: float) -> float:
    return halo_birth_init_time + halo_stars_formation_time * random()


def thin_disk_star_birth_time(*,
                              thin_disk_birth_init_time: float,
                              time_bin: int,
                              time_increment: float) -> float:
    return (thin_disk_birth_init_time
            + time_bin * time_increment
            + time_increment * random())


def get_galactic_structure_type(*,
                                thick_disk_stars_fraction: float,
                                halo_stars_fraction: float
                                ) -> GalacticStructureType:
    random_number = random()

    if random_number <= thick_disk_stars_fraction:
        return GalacticStructureType.thick

    if (thick_disk_stars_fraction < random_number
            <= thick_disk_stars_fraction + halo_stars_fraction):
        return GalacticStructureType.halo

    return GalacticStructureType.thin


# TODO: move this to 'polar' module
# thin_disk_scale_height_kpc: float = 0.25,
# thick_disk_scale_height_kpc: float = 0.9,
def z_coordinate(*,
                 scale_height: float,
                 sector_radius_kpc: float) -> float:
    # TODO: implement function for inverse transform sampling
    # Inverse transform sampling for y = exp(-z / H)
    coordinate = (-scale_height * math.log(
            1. - random() * (1.0 - math.exp(-sector_radius_kpc
                                            / scale_height))))
    # TODO: find a better way to assign a random sign
    random_sign = float(1. - 2. * int(2.0 * random()))

    return coordinate * random_sign
