import math
from functools import partial
from typing import Callable

import numpy as np
import pandas as pd

from alcor.types import (GeneratorType,
                         UnitRangeGeneratorType)

random_signs = partial(np.random.choice, [-1, 1])


def assign_polar_coordinates(
        stars: pd.DataFrame,
        *,
        sector_radius: float,
        scale_length: float,
        solar_galactocentric_distance: float,
        thin_disk_scale_height: float,
        thick_disk_scale_height: float,
        halo_core_radius: float = 5.,
        generator: GeneratorType = np.random.uniform,
        unit_range_generator: UnitRangeGeneratorType = np.random.rand,
        signs_generator: Callable[[int], np.ndarray] = random_signs
        ) -> pd.DataFrame:
    min_sector_radius = solar_galactocentric_distance - sector_radius
    max_sector_radius = solar_galactocentric_distance + sector_radius
    squared_max_sector_radius = max_sector_radius ** 2
    squared_min_sector_radius = min_sector_radius ** 2
    squared_radii_difference = (squared_max_sector_radius
                                - squared_min_sector_radius)
    angle_covering_sector = 2. * math.asin(sector_radius
                                           / solar_galactocentric_distance)
    radial_distrib_max = math.exp(-min_sector_radius / scale_length)

    stars['theta_cylindrical'] = thetas_cylindrical(
            size=stars.shape[0],
            angle_covering_sector=angle_covering_sector,
            generator=generator)

    halo_stars_mask = stars['galactic_disk_type'] == 'halo'

    halo_stars = stars[halo_stars_mask]
    disks_stars = stars[~halo_stars_mask]

    halo_stars['r_cylindrical'] = halo_r_cylindrical_coordinates(
            size=halo_stars.shape[0],
            min_sector_radius=min_sector_radius,
            max_sector_radius=max_sector_radius,
            halo_core_radius=halo_core_radius,
            squared_min_sector_radius=squared_min_sector_radius,
            squared_radii_difference=squared_radii_difference,
            generator=unit_range_generator)

    disks_stars['r_cylindrical'] = disks_r_cylindrical(
            size=disks_stars.shape[0],
            min_sector_radius=min_sector_radius,
            max_sector_radius=max_sector_radius,
            scale_length=scale_length,
            radial_distrib_max=radial_distrib_max,
            squared_min_sector_radius=squared_min_sector_radius,
            squared_radii_difference=squared_radii_difference,
            generator=generator)

    thin_disk_stars_mask = disks_stars['galactic_disk_type'] == 'thin'

    thin_disk_stars = disks_stars[thin_disk_stars_mask]
    thick_disk_stars = disks_stars[~thin_disk_stars_mask]

    halo_stars['z_coordinate'] = halo_z_coordinates(
            angle_covering_sector=angle_covering_sector,
            r_cylindrical=disks_stars['r_cylindrical'].values,
            generator=generator)

    thin_disk_stars['z_coordinate'] = disk_z_coordinates(
            size=thin_disk_stars.shape[0],
            scale_height=thin_disk_scale_height,
            sector_radius=sector_radius,
            generator=unit_range_generator,
            signs_generator=signs_generator)

    thick_disk_stars['z_coordinate'] = disk_z_coordinates(
            size=thick_disk_stars.shape[0],
            scale_height=thick_disk_scale_height,
            sector_radius=sector_radius,
            generator=unit_range_generator,
            signs_generator=signs_generator)

    # TODO: here we had filtering out too close stars to the Sun.
    # Do we really need to do that? (alpha_centauri_distance = 1.5e-6)

    return pd.concat(halo_stars, thin_disk_stars, thick_disk_stars)


def disk_z_coordinates(*,
                       size: int,
                       scale_height: float,
                       sector_radius: float,
                       generator: UnitRangeGeneratorType,
                       signs_generator: Callable[[int], np.ndarray]
                       ) -> np.ndarray:
    abs_z_coordinates = (-scale_height
                         * np.log(1. - generator(size)
                                  * (1. - math.exp(-sector_radius
                                                   / scale_height))))
    return np.multiply(abs_z_coordinates,
                       signs_generator(size=size))


def halo_z_coordinates(*,
                       angle_covering_sector: float,
                       r_cylindrical: np.ndarray,
                       generator: GeneratorType) -> np.ndarray:
    random_angles = generator(low=-angle_covering_sector / 2.,
                              high=angle_covering_sector / 2.,
                              size=r_cylindrical.size)
    return r_cylindrical * np.sin(random_angles)


def disks_r_cylindrical(*,
                        size: int,
                        min_sector_radius: float,
                        max_sector_radius: float,
                        scale_length: float,
                        radial_distrib_max: float,
                        squared_min_sector_radius: float,
                        squared_radii_difference: float,
                        generator: GeneratorType) -> np.ndarray:
    radii_tries = disks_stars_radii_tries(
            size=size,
            min_sector_radius=min_sector_radius,
            max_sector_radius=max_sector_radius,
            scale_length=scale_length,
            radial_distrib_max=radial_distrib_max,
            generator=generator)

    # TODO: put this in a function
    # Inverse transform sampling method for generating stars
    # uniformly in a circle sector in polar coordinates
    random_values = (radii_tries - min_sector_radius) / (max_sector_radius
                                                         - min_sector_radius)
    return np.sqrt(squared_radii_difference * random_values
                   + squared_min_sector_radius)


def halo_r_cylindrical_coordinates(*,
                                   size: int,
                                   min_sector_radius: float,
                                   max_sector_radius: float,
                                   halo_core_radius: float,
                                   squared_min_sector_radius: float,
                                   squared_radii_difference: float,
                                   generator: UnitRangeGeneratorType
                                   ) -> np.ndarray:
    radii_tries = halo_stars_radii_tries(size=size,
                                         min_sector_radius=min_sector_radius,
                                         max_sector_radius=max_sector_radius,
                                         halo_core_radius=halo_core_radius,
                                         generator=generator)

    # Inverse transform sampling method for generating stars
    # uniformly in a circle sector in polar coordinates
    random_values = (radii_tries - min_sector_radius) / (max_sector_radius
                                                         - min_sector_radius)
    return np.sqrt(squared_radii_difference * random_values
                   + squared_min_sector_radius)


def halo_stars_radii_tries(*,
                           size: int,
                           min_sector_radius: float,
                           max_sector_radius: float,
                           halo_core_radius: float,
                           generator: UnitRangeGeneratorType) -> np.ndarray:
    """
    Inverse transform sampling for halo distribution.
    See (4) at "Simulating Gaia performances on white
    dwarfs" by Torres et al. 2005
    """
    min_atan = math.atan(min_sector_radius / halo_core_radius)
    max_atan = math.atan(max_sector_radius / halo_core_radius)
    return np.multiply(halo_core_radius,
                       np.tan(generator(size) * (max_atan - min_atan)
                              + min_atan))


def disks_stars_radii_tries(*,
                            size: int,
                            min_sector_radius: float,
                            max_sector_radius: float,
                            scale_length: float,
                            radial_distrib_max: float,
                            generator: GeneratorType) -> np.ndarray:
    disk_stars_radii_function = partial(disk_stars_radius,
                                        scale_length=scale_length)
    disks_radii_distribution = partial(monte_carlo_generator,
                                       function=disk_stars_radii_function,
                                       min_x=min_sector_radius,
                                       max_x=max_sector_radius,
                                       max_y=radial_distrib_max,
                                       generator=generator)
    return np.array([disks_radii_distribution()
                     for _ in range(size)])


def disk_stars_radius(value: float,
                      *,
                      scale_length: float) -> float:
    return math.exp(-value / scale_length)


def monte_carlo_generator(*,
                          function: Callable[[float], float],
                          min_x: float,
                          max_x: float,
                          max_y: float,
                          generator: Callable[[float, float], float],
                          min_y: float = 0.,
                          max_iterations_count: int = 10 ** 9
                          ) -> float:
    for _ in range(max_iterations_count):
        x_value = generator(min_x, max_x)
        if generator(min_y, max_y) <= function(x_value):
            return x_value
    else:
        raise OverflowError('Exceeded maximum number of iterations '
                            'in Monte Carlo generator for "{function}"'
                            .format(function=function.__qualname__))


def thetas_cylindrical(*,
                       size: int,
                       angle_covering_sector: float,
                       generator: GeneratorType) -> np.ndarray:
    thetas = generator(low=-angle_covering_sector / 2.,
                       high=angle_covering_sector / 2.,
                       size=size)
    thetas[thetas < 0.] += 2. * np.pi
    return thetas
