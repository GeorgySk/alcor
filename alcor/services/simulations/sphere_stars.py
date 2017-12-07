import math
import random

import numpy as np
import pandas as pd

from alcor.models.star import GalacticStructureType


# TODO: take these consts to settings file
def generate_stars(*,
                   max_stars_count: int = 6000000,
                   time_bins_count: int = 5000,
                   thin_disk_age: float = 9.2,
                   thick_disk_age: float = 12.,
                   halo_age: float = 14.,
                   halo_stars_formation_time: float = 1.,
                   burst_age: float = 0.6,
                   burst_formation_factor: float = 5.,
                   star_formation_rate_param: float = 25.,
                   thick_disk_sfr_param: float = 2.,
                   thick_disk_stars_fraction: float = 0.8,
                   halo_stars_fraction: float = 0.05,
                   sector_radius_kpc: float = 0.05,
                   mass_reduction_factor: float = 0.03,
                   initial_mass_function_param: float = -2.35) -> pd.DataFrame:
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

    thin_disk_stars = generate_thin_disk_stars(
            thin_disk_birth_init_time=thin_disk_birth_init_time,
            max_age=max_age,
            time_bins_count=time_bins_count,
            burst_init_time=burst_init_time,
            birth_rate=birth_rate,
            burst_birth_rate=burst_birth_rate,
            initial_mass_function_parameter=initial_mass_function_param,
            time_increment=time_increment,
            max_stars_count=max_stars_count)

    thick_disk_stars = generate_thick_disk_stars(
            thin_disk_stars_fraction=(1. - thick_disk_stars_fraction
                                      - halo_stars_fraction),
            thick_disk_stars_fraction=thick_disk_stars_fraction,
            thin_disk_stars_count=thin_disk_stars.shape[0],
            initial_mass_function_parameter=initial_mass_function_param,
            thick_disk_age=thick_disk_age,
            thick_disk_birth_init_time=thick_disk_birth_init_time,
            thick_disk_max_sfr=thick_disk_max_sfr,
            thick_disk_sfr_param=thick_disk_sfr_param)

    halo_stars = generate_halo_stars(
            thin_disk_stars_count=thin_disk_stars.shape[0],
            halo_stars_fraction=halo_stars_fraction,
            thin_disk_stars_fraction=(1. - thick_disk_stars_fraction
                                      - halo_stars_fraction),
            initial_mass_function_parameter=initial_mass_function_param,
            halo_birth_init_time=halo_birth_init_time,
            halo_stars_formation_time=halo_stars_formation_time)

    return pd.concat([thin_disk_stars, thick_disk_stars, halo_stars])


def generate_halo_stars(*,
                        thin_disk_stars_count: float,
                        halo_stars_fraction: float,
                        thin_disk_stars_fraction: float,
                        initial_mass_function_parameter: float,
                        halo_birth_init_time: float,
                        halo_stars_formation_time: float) -> pd.DataFrame:
    halo_stars_count = int(thin_disk_stars_count
                           * halo_stars_fraction
                           / thin_disk_stars_fraction)
    progenitors_masses = []
    birth_times = []

    for _ in range(halo_stars_count):
        progenitors_masses.append(
                initial_star_mass_by_salpeter(initial_mass_function_parameter))
        birth_times.append(halo_star_birth_time(
                    birth_initial_time=halo_birth_init_time,
                    formation_time=halo_stars_formation_time))

    halo_stars = pd.DataFrame(dict(progenitor_mass=progenitors_masses,
                                   birth_time=birth_times))
    halo_stars['galactic_structure_type'] = GalacticStructureType.halo

    return halo_stars


def generate_thick_disk_stars(*,
                              thin_disk_stars_fraction: float,
                              thick_disk_stars_fraction: float,
                              thin_disk_stars_count: int,
                              initial_mass_function_parameter: float,
                              thick_disk_age: float,
                              thick_disk_birth_init_time: float,
                              thick_disk_max_sfr: float,
                              thick_disk_sfr_param: float) -> pd.DataFrame:
    thick_disk_stars_count = int(thin_disk_stars_count
                                 * thick_disk_stars_fraction
                                 / thin_disk_stars_fraction)
    progenitors_masses = []
    birth_times = []

    for _ in range(thick_disk_stars_count):
        progenitors_masses.append(
                initial_star_mass_by_salpeter(initial_mass_function_parameter))
        birth_times.append(thick_disk_star_birth_time(
                age=thick_disk_age,
                birth_initial_time=thick_disk_birth_init_time,
                max_formation_rate=thick_disk_max_sfr,
                formation_rate_parameter=thick_disk_sfr_param))

    thick_disk_stars = pd.DataFrame(dict(progenitor_mass=progenitors_masses,
                                         birth_time=birth_times))
    thick_disk_stars['galactic_structure_type'] = GalacticStructureType.thick

    return thick_disk_stars


def generate_thin_disk_stars(*,
                             thin_disk_birth_init_time: float,
                             max_age: float,
                             time_bins_count: int,
                             burst_init_time: float,
                             birth_rate: float,
                             burst_birth_rate: float,
                             initial_mass_function_parameter: float,
                             time_increment: float,
                             max_stars_count: int) -> pd.DataFrame:
    time_bins_initial_times = np.linspace(start=thin_disk_birth_init_time,
                                          stop=max_age,
                                          num=time_bins_count,
                                          endpoint=False)
    birth_rates = get_birth_rates(time_bins_initial_times,
                                  burst_init_time=burst_init_time,
                                  birth_rate=birth_rate,
                                  burst_birth_rate=burst_birth_rate)

    progenitors_masses = []
    birth_times = []

    for bin_initial_time, birth_rate in np.column_stack(
            (time_bins_initial_times, birth_rates)):
        total_bin_mass = 0.

        while total_bin_mass < birth_rate:
            birth_times.append(thin_disk_star_birth_time(
                    bin_initial_time=bin_initial_time,
                    time_increment=time_increment))

            star_mass = initial_star_mass_by_salpeter(
                    initial_mass_function_parameter)
            progenitors_masses.append(star_mass)

            total_bin_mass += star_mass

        if len(progenitors_masses) > max_stars_count:
            raise OverflowError('Number of stars is too high - '
                                'decrease mass reduction factor.')

    thin_disk_stars = pd.DataFrame(dict(progenitor_mass=progenitors_masses,
                                        birth_time=birth_times))
    thin_disk_stars['galactic_structure_type'] = GalacticStructureType.thin

    return thin_disk_stars


def get_birth_rates(times: np.ndarray,
                    *,
                    burst_init_time: float,
                    birth_rate: float,
                    burst_birth_rate: float) -> np.ndarray:
    burst_times_mask = times >= burst_init_time
    return np.piecewise(times,
                        [~burst_times_mask, burst_times_mask],
                        [birth_rate, burst_birth_rate])


def normalization_const(*,
                        star_formation_rate_param: float,
                        thin_disk_age_gyr: float,
                        sigma: float = 51.  # TODO: what is sigma?
                        ) -> float:
    return sigma / (star_formation_rate_param
                    * (1 - math.exp(-thin_disk_age_gyr
                                    / star_formation_rate_param)))


# TODO: implement inverse transform sampling
def initial_star_mass_by_salpeter(exponent: float,
                                  *,
                                  min_mass: float = 0.4,
                                  max_mass: float = 50.) -> float:
    y_max = min_mass ** exponent

    while True:
        mass = random.uniform(min_mass, max_mass)
        y_mass = mass ** exponent
        if random.uniform(0, y_max) <= y_mass:
            return mass


def thick_disk_star_birth_time(*,
                               age: float,
                               formation_rate_parameter: float,
                               max_formation_rate: float,
                               birth_initial_time: float) -> float:
    """
    Return birth time of a thick disk star by using Monte Carlo method.
    SFR - star formation rate. More info at:
    https://www.google.es/search?q=star+formation+rate
    """
    while True:
        time_try = age * random.random()
        time_try_sfr = time_try * math.exp(-time_try
                                           / formation_rate_parameter)
        sfr_try = max_formation_rate * random.random()
        if sfr_try <= time_try_sfr:
            return time_try + birth_initial_time


def halo_star_birth_time(*,
                         birth_initial_time: float,
                         formation_time: float) -> float:
    return birth_initial_time + formation_time * random.random()


def thin_disk_star_birth_time(*,
                              bin_initial_time: float,
                              time_increment: float) -> float:
    return bin_initial_time + time_increment * random.random()


def get_galactic_structure_type(*,
                                thick_disk_stars_fraction: float,
                                halo_stars_fraction: float
                                ) -> GalacticStructureType:
    random_number = random.random()

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
            1. - random.random() * (1.0 - math.exp(-sector_radius_kpc
                                                   / scale_height))))
    # TODO: find a better way to assign a random sign
    random_sign = float(1. - 2. * int(2.0 * random.random()))

    return coordinate * random_sign
