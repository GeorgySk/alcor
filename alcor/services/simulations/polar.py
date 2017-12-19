import math
from functools import partial
from random import random
from typing import (Callable,
                    Tuple)

import numpy as np
import pandas as pd


def assign_polar_coordinates(
        stars: pd.DataFrame,
        sector_radius: float,
        scale_length: float,
        solar_galactocentric_distance: float,
        thin_disk_scale_height: float,
        thick_disk_scale_height: float,
        halo_core_radius: float = 5.,
        alpha_centauri_distance: float = 1.5e-6,
        generator: Callable[[float, float, float], np.ndarray] = (
                np.random.uniform),
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

    halo_stars['r_cylindrical'] = halo_r_cylindrical(
            size=halo_stars.shape[0],
            min_sector_radius=min_sector_radius,
            max_sector_radius=max_sector_radius,
            halo_core_radius=halo_core_radius,
            squared_min_sector_radius=squared_min_sector_radius,
            squared_radii_difference=squared_radii_difference)

    disks_stars['r_cylindrical'] = disks_r_cylindrical(
            size=disks_stars.shape[0],
            min_sector_radius=min_sector_radius,
            max_sector_radius=max_sector_radius,
            scale_length=scale_length,
            radial_distrib_max=radial_distrib_max,
            squared_min_sector_radius=squared_min_sector_radius,
            squared_radii_difference=squared_radii_difference)

    set_z_coordinate(stars,
                     angle_covering_sector=angle_covering_sector,
                     thin_disk_scale_height=thin_disk_scale_height,
                     thick_disk_scale_height=thick_disk_scale_height,
                     sector_radius=sector_radius)
    stars = replace_too_close_to_the_sun_stars(
            stars,
            min_distance=alpha_centauri_distance,
            solar_galactocentric_distance=solar_galactocentric_distance)

    return stars


def disks_r_cylindrical(*,
                        size: int,
                        min_sector_radius: float,
                        max_sector_radius: float,
                        scale_length: float,
                        radial_distrib_max: float,
                        squared_min_sector_radius: float,
                        squared_radii_difference: float) -> np.ndarray:
    radii_tries = disks_stars_radii_tries(
            size=size,
            min_sector_radius=min_sector_radius,
            max_sector_radius=max_sector_radius,
            scale_length=scale_length,
            radial_distrib_max=radial_distrib_max)

    # TODO: put this in a function
    # Inverse transform sampling method for generating stars
    # uniformly in a circle sector in polar coordinates
    random_values = (radii_tries - min_sector_radius) / (max_sector_radius
                                                         - min_sector_radius)
    return np.sqrt(squared_radii_difference * random_values
                   + squared_min_sector_radius)


def halo_r_cylindrical(*,
                       size: int,
                       min_sector_radius: float,
                       max_sector_radius: float,
                       halo_core_radius: float,
                       squared_min_sector_radius: float,
                       squared_radii_difference: float) -> np.ndarray:
    radii_tries = halo_stars_radii_tries(size=size,
                                         min_sector_radius=min_sector_radius,
                                         max_sector_radius=max_sector_radius,
                                         halo_core_radius=halo_core_radius,
                                         generator=np.random.rand)

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
                           generator: Callable[[Tuple[int, ...]], np.ndarray]
                           ) -> np.ndarray:
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
                            radial_distrib_max: float) -> np.ndarray:
    disk_stars_radii_function = partial(disk_stars_radius,
                                        scale_length=scale_length)

    disks_radii_distribution = partial(monte_carlo_generator,
                                       function=disk_stars_radii_function,
                                       min_x=min_sector_radius,
                                       max_x=max_sector_radius,
                                       max_y=radial_distrib_max,
                                       generator=np.random.uniform)

    return np.array([disks_radii_distribution()
                     for _ in range(size)])


# TODO: figure out workflow for these stars
def replace_too_close_to_the_sun_stars(stars: pd.DataFrame,
                                       min_distance: float,
                                       solar_galactocentric_distance: float
                                       ) -> pd.DataFrame:
    return stars


def set_z_coordinate(stars: pd.DataFrame,
                     *,
                     angle_covering_sector: float,
                     thin_disk_scale_height: float,
                     thick_disk_scale_height: float,
                     sector_radius: float) -> None:
    stars['z_coordinate'] = np.nan

    halo_stars_mask = stars['galactic_disk_type'] == 'halo'
    halo_stars_count = stars[halo_stars_mask].shape[0]

    random_angles = (angle_covering_sector * np.random.rand(halo_stars_count)
                     - angle_covering_sector / 2)
    stars.loc[halo_stars_mask, 'z_coordinate'] = (stars['r_cylindrical']
                                                  * np.sin(random_angles))

    thin_disk_stars_mask = stars['galactic_disk_type'] == 'thin'
    thin_disk_stars_count = stars[thin_disk_stars_mask].shape[0]

    abs_z_coordinates = (-thin_disk_scale_height
                         * np.log(1.
                                  - np.random.rand(thin_disk_stars_count)
                                    * (1.
                                       - math.exp(-sector_radius
                                                  / thin_disk_scale_height))))
    stars.loc[thin_disk_stars_mask, 'z_coordinate'] = (
        abs_z_coordinates * random.choice([1., -1.]))

    thick_disk_stars_mask = stars['galactic_disk_type'] == 'thick'
    thick_disk_stars_count = stars[thick_disk_stars_mask].shape[0]

    abs_z_coordinates = (-thick_disk_scale_height
                         * np.log(1.
                                  - np.random.rand(thick_disk_stars_count)
                                  * (1.
                                     - math.exp(-sector_radius
                                                / thick_disk_scale_height))))
    stars.loc[thick_disk_stars_mask, 'z_coordinate'] = (
        abs_z_coordinates * random.choice([1., -1.]))


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


def thetas_cylindrical(
        *,
        size: int,
        angle_covering_sector: float,
        generator: Callable[[float, float, float], np.ndarray]) -> np.ndarray:
    thetas = generator(low=-angle_covering_sector / 2.,
                       high=angle_covering_sector / 2.,
                       size=size)
    thetas[thetas < 0.] += 2. * np.pi

    return thetas
