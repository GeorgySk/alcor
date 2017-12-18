import math
from random import random
from typing import Callable

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
    sector_diameter = max_sector_radius - min_sector_radius
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
    set_r_cylindrical(stars,
                      min_sector_radius=min_sector_radius,
                      max_sector_radius=max_sector_radius,
                      halo_core_radius=halo_core_radius,
                      sector_diameter=sector_diameter,
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


def set_r_cylindrical(stars: pd.DataFrame,
                      min_sector_radius: float,
                      max_sector_radius: float,
                      halo_core_radius: float,
                      sector_diameter: float,
                      scale_length: float,
                      radial_distrib_max: float,
                      squared_min_sector_radius: float,
                      squared_radii_difference: float) -> None:
    radii_tries = get_radii_tries(
            stars,
            min_sector_radius=min_sector_radius,
            max_sector_radius=max_sector_radius,
            halo_core_radius=halo_core_radius,
            sector_diameter=sector_diameter,
            scale_length=scale_length,
            radial_distrib_max=radial_distrib_max)
    # Inverse transform sampling method for generating stars
    # uniformly in a circle sector in polar coordinates
    random_values = (radii_tries - min_sector_radius) / sector_diameter
    stars['r_cylindrical'] = np.sqrt(squared_radii_difference * random_values
                                     + squared_min_sector_radius)


def get_radii_tries(stars: pd.DataFrame,
                    min_sector_radius: float,
                    max_sector_radius: float,
                    halo_core_radius: float,
                    sector_diameter: float,
                    scale_length: float,
                    radial_distrib_max: float) -> np.ndarray:
    radii_tries = np.empty(stars.shape[0])
    halo_stars_mask = stars['galactic_disk_type'] == 'halo'
    disks_stars_mask = ~halo_stars_mask
    halo_stars_count = stars[halo_stars_mask].shape[0]
    disk_stars_count = stars[disks_stars_mask].shape[0]

    # Inverse transform sampling for halo distribution.
    # See (4) at "Simulating Gaia performances on white
    # dwarfs" by Torres et al. 2005
    min_atan = math.atan(min_sector_radius / halo_core_radius)
    max_atan = math.atan(max_sector_radius / halo_core_radius)
    radii_tries[halo_stars_mask] = (halo_core_radius * np.tan(
            np.random.rand(halo_stars_count) * (max_atan - min_atan)
            + min_atan))

    # Accepting-rejecting method
    disk_stars_radii_tries = []
    while len(disk_stars_radii_tries) != disk_stars_count:
        radius_try = min_sector_radius + sector_diameter * random()
        radius_try_distrib = math.exp(-radius_try / scale_length)
        radial_distrib_random = radial_distrib_max * random()
        if radial_distrib_random <= radius_try_distrib:
            disk_stars_radii_tries.append(radius_try)
    radii_tries[disks_stars_mask] = disk_stars_radii_tries

    return radii_tries


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
