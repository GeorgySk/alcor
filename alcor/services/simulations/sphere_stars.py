import math
import random
from functools import partial
from itertools import chain
from typing import (Any,
                    Callable,
                    Iterable,
                    Iterator,
                    Tuple)

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
                   thick_disk_formation_rate_param: float = 2.,
                   thick_disk_stars_fraction: float = 0.8,
                   halo_stars_fraction: float = 0.05,
                   initial_mass_function_param: float = -2.35,
                   sector_radius_kpc: float = 0.05,
                   burst_formation_factor: float = 5.,
                   star_formation_rate_param: float = 25.,
                   mass_reduction_factor: float = 0.03,
                   generator: Callable[[float, float], float] = random.uniform
                   ) -> pd.DataFrame:
    max_age = max(thin_disk_age, thick_disk_age, halo_age)
    thin_disk_stars_fraction = (1. - thick_disk_stars_fraction
                                - halo_stars_fraction)

    thin_disk_stars = generate_thin_disk_stars(
            max_age=max_age,
            time_bins_count=time_bins_count,
            burst_age=burst_age,
            initial_mass_function_parameter=initial_mass_function_param,
            disk_age=thin_disk_age,
            max_stars_count=max_stars_count,
            sector_radius_kpc=sector_radius_kpc,
            burst_formation_factor=burst_formation_factor,
            star_formation_rate_param=star_formation_rate_param,
            mass_reduction_factor=mass_reduction_factor,
            generator=generator)

    thick_disk_stars_count = int(thin_disk_stars.shape[0]
                                 * thick_disk_stars_fraction
                                 / thin_disk_stars_fraction)
    thick_disk_stars = generate_thick_disk_stars(
            stars_count=thick_disk_stars_count,
            initial_mass_function_parameter=initial_mass_function_param,
            disk_age=thick_disk_age,
            birth_initial_time=max_age - thick_disk_age,
            formation_rate_parameter=thick_disk_formation_rate_param,
            generator=generator)

    halo_stars_count = int(thin_disk_stars.shape[0] * halo_stars_fraction
                           / thin_disk_stars_fraction)
    halo_stars = generate_halo_stars(
            stars_count=halo_stars_count,
            initial_mass_function_parameter=initial_mass_function_param,
            birth_initial_time=max_age - halo_age,
            formation_time=halo_stars_formation_time,
            generator=generator)

    return pd.concat([thin_disk_stars, thick_disk_stars, halo_stars])


def generate_halo_stars(*,
                        stars_count: int,
                        initial_mass_function_parameter: float,
                        birth_initial_time: float,
                        formation_time: float,
                        generator: Callable[[float, float], float]
                        ) -> pd.DataFrame:
    progenitors_masses = [
        initial_star_mass_by_salpeter(initial_mass_function_parameter,
                                      generator=generator)
        for _ in range(stars_count)]
    birth_times = [halo_star_birth_time(
            birth_initial_time=birth_initial_time,
            formation_time=formation_time,
            generator=generator)
        for _ in range(stars_count)]

    return pd.DataFrame(dict(
            progenitor_mass=progenitors_masses,
            birth_time=birth_times,
            galactic_structure_type=GalacticStructureType.halo))


def generate_thick_disk_stars(*,
                              stars_count: int,
                              initial_mass_function_parameter: float,
                              disk_age: float,
                              birth_initial_time: float,
                              formation_rate_parameter: float,
                              generator: Callable[[float, float], float]
                              ) -> pd.DataFrame:
    # This can be proved by taking derivative from y = t * exp(-t / tau)
    max_formation_rate = formation_rate_parameter / math.e

    progenitors_masses = [
        initial_star_mass_by_salpeter(initial_mass_function_parameter,
                                      generator=generator)
        for _ in range(stars_count)]
    birth_times = [thick_disk_star_birth_time(
            disk_age=disk_age,
            birth_initial_time=birth_initial_time,
            max_formation_rate=max_formation_rate,
            formation_rate_parameter=formation_rate_parameter,
            generator=generator)
        for _ in range(stars_count)]

    return pd.DataFrame(dict(
            progenitor_mass=progenitors_masses,
            birth_time=birth_times,
            galactic_structure_type=GalacticStructureType.thick))


def generate_thin_disk_stars(*,
                             max_age: float,
                             time_bins_count: int,
                             burst_age: float,
                             initial_mass_function_parameter: float,
                             disk_age: float,
                             max_stars_count: int,
                             sector_radius_kpc: float,
                             burst_formation_factor: float,
                             star_formation_rate_param: float,
                             mass_reduction_factor: float,
                             generator: Callable[[float, float], float]
                             ) -> pd.DataFrame:
    time_increment = disk_age / time_bins_count
    sector_area = math.pi * sector_radius_kpc ** 2
    birth_rate = (
        time_increment * sector_area * 1E6  # TODO: what is 1E6?
        * mass_reduction_factor
        * normalization_const(
                star_formation_rate_param=star_formation_rate_param,
                thin_disk_age_gyr=disk_age))
    burst_birth_rate = birth_rate * burst_formation_factor
    burst_initial_time = max_age - burst_age
    birth_initial_time = max_age - disk_age

    time_bins_initial_times = np.linspace(start=birth_initial_time,
                                          stop=max_age,
                                          num=time_bins_count,
                                          endpoint=False)
    birth_rates = get_birth_rates(time_bins_initial_times,
                                  burst_initial_time=burst_initial_time,
                                  birth_rate=birth_rate,
                                  burst_birth_rate=burst_birth_rate)

    star_mass_generator = partial(initial_star_mass_by_salpeter,
                                  initial_mass_function_parameter,
                                  generator=generator)
    birth_time_generators = [partial(thin_disk_star_birth_time,
                                     bin_initial_time=bin_initial_time,
                                     time_increment=time_increment,
                                     generator=generator)
                             for bin_initial_time in time_bins_initial_times]

    columns = birth_rates, birth_time_generators
    progenitors_masses, birth_times = zip(*chain.from_iterable(
            progenitors_masses_births_times(
                    birth_rate=birth_rate,
                    birth_time_generator=birth_time_generator,
                    star_mass_generator=star_mass_generator,
                    max_stars_count=max_stars_count)
            for birth_rate, birth_time_generator in np.column_stack(columns)))

    return pd.DataFrame(dict(
            progenitor_mass=progenitors_masses,
            birth_time=birth_times,
            galactic_structure_type=GalacticStructureType.thin))


def progenitors_masses_births_times(*,
                                    birth_rate: float,
                                    birth_time_generator: Callable[[], float],
                                    star_mass_generator: Callable[[], float],
                                    max_stars_count: int
                                    ) -> Iterator[Tuple[float, float]]:
    total_bin_mass = 0.
    for _ in range(max_stars_count):
        birth_time = birth_time_generator()
        star_mass = star_mass_generator()
        yield birth_time, star_mass

        total_bin_mass += star_mass
        if total_bin_mass >= birth_rate:
            return
    else:
        raise OverflowError('Number of stars is too high - '
                            'decrease mass reduction factor.')


# TODO: find out the meaning
def get_birth_rates(times: np.ndarray,
                    *,
                    burst_initial_time: float,
                    birth_rate: float,
                    burst_birth_rate: float) -> np.ndarray:
    burst_times_mask = times >= burst_initial_time
    return np.piecewise(times,
                        [~burst_times_mask, burst_times_mask],
                        [birth_rate, burst_birth_rate])


# TODO: find out the meaning
def normalization_const(*,
                        star_formation_rate_param: float,
                        thin_disk_age_gyr: float,
                        sigma: float = 51.  # TODO: what is sigma?
                        ) -> float:
    return sigma / (star_formation_rate_param
                    * (1 - math.exp(-thin_disk_age_gyr
                                    / star_formation_rate_param)))


# TODO: implement inverse transform sampling
def initial_star_mass_by_salpeter(
        exponent: float,
        *,
        min_mass: float = 0.4,
        max_mass: float = 50.,
        generator: Callable[[float, float], float]) -> float:
    y_max = min_mass ** exponent

    while True:
        mass = generator(min_mass, max_mass)
        y_mass = mass ** exponent
        if generator(0, y_max) <= y_mass:
            return mass


def thick_disk_star_birth_time(
        *,
        disk_age: float,
        formation_rate_parameter: float,
        max_formation_rate: float,
        birth_initial_time: float,
        generator: Callable[[float, float], float]
        ) -> float:
    """
    Returns birth time of a thick disk star
    calculated by using Monte Carlo method.
    SFR - star formation rate. More info at:
    https://www.google.es/search?q=star+formation+rate
    """
    while True:
        time_try = generator(0, disk_age)
        time_try_formation_rate = time_try * math.exp(
                -time_try / formation_rate_parameter)
        if generator(0, max_formation_rate) <= time_try_formation_rate:
            return time_try + birth_initial_time


def halo_star_birth_time(
        *,
        birth_initial_time: float,
        formation_time: float,
        generator: Callable[[float, float], float]
        ) -> float:
    return generator(birth_initial_time, birth_initial_time + formation_time)


def thin_disk_star_birth_time(
        *,
        bin_initial_time: float,
        time_increment: float,
        generator: Callable[[float, float], float]
        ) -> float:
    return generator(bin_initial_time, bin_initial_time + time_increment)


# TODO: move this to 'polar' module
# thin_disk_scale_height_kpc: float = 0.25,
# thick_disk_scale_height_kpc: float = 0.9,
def z_coordinate(*,
                 scale_height: float,
                 sector_radius_kpc: float,
                 rng: Callable[[], float] = random.random,
                 rng_choice: Callable[[Iterable[Any]], Any] = random.choice
                 ) -> float:
    # TODO: implement function for inverse transform sampling
    # Inverse transform sampling for y = exp(-z / H)
    coordinate = (-scale_height * math.log(
            1. - rng() * (1.0 - math.exp(-sector_radius_kpc / scale_height))))
    random_sign = rng_choice([-1, 1])

    return coordinate * random_sign
