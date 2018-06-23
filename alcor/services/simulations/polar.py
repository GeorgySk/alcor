import math
from functools import partial
from typing import Union

import numpy as np
import pandas as pd
import sympy as sym

from alcor.models.star import GalacticDiskType
from alcor.services import inverse_transform
from alcor.types import GeneratorType

random_signs = partial(np.random.choice, [-1, 1])


def assign_polar_coordinates(
        stars: pd.DataFrame,
        *,
        sector_radius: float,
        scale_length: float,
        solar_galactocentric_distance: float,
        thin_disk_scale_height: float,
        thick_disk_scale_height: float,
        halo_core_radius: float,
        generator: GeneratorType = np.random.uniform) -> pd.DataFrame:
    min_sector_radius = solar_galactocentric_distance - sector_radius
    max_sector_radius = solar_galactocentric_distance + sector_radius

    angle_covering_sector = 2. * math.asin(sector_radius
                                           / solar_galactocentric_distance)
    stars['theta_cylindrical'] = thetas_cylindrical(
            size=stars.shape[0],
            angle_covering_sector=angle_covering_sector,
            generator=generator)

    halo_stars_mask = stars['galactic_disk_type'] == GalacticDiskType.halo

    halo_stars = stars[halo_stars_mask].copy()
    disks_stars = stars[~halo_stars_mask].copy()

    halo_star_distances = halo_star_distance(min_distance=min_sector_radius,
                                             max_distance=max_sector_radius,
                                             core_radius=halo_core_radius)
    halo_stars['r_cylindrical'] = inverse_transform.sample(
            halo_star_distances,
            size=halo_stars.shape[0])

    disk_star_distances = disk_star_distance(min_distance=min_sector_radius,
                                             max_distance=max_sector_radius,
                                             scale_length=scale_length)
    disks_stars['r_cylindrical'] = inverse_transform.sample(
            disk_star_distances,
            size=disks_stars.shape[0])

    thin_disk_stars_mask = (disks_stars['galactic_disk_type']
                            == GalacticDiskType.thin)

    thin_disk_stars = disks_stars[thin_disk_stars_mask].copy()
    thick_disk_stars = disks_stars[~thin_disk_stars_mask].copy()

    halo_stars['z_coordinate'] = halo_z_coordinates(
            angle_covering_sector=angle_covering_sector,
            r_cylindrical=halo_stars['r_cylindrical'].values,
            generator=generator)

    thin_disk_z_coordinates = disk_z_coordinate(
            radius=sector_radius,
            scale_height=thin_disk_scale_height)
    thin_disk_stars['z_coordinate'] = inverse_transform.sample(
            thin_disk_z_coordinates,
            size=thin_disk_stars.shape[0])

    thick_disk_z_coordinates = disk_z_coordinate(
            radius=sector_radius,
            scale_height=thick_disk_scale_height)
    thick_disk_stars['z_coordinate'] = inverse_transform.sample(
            thick_disk_z_coordinates,
            size=thick_disk_stars.shape[0])

    # TODO: here we had filtering out too close stars to the Sun.
    # Do we really need to do that? (alpha_centauri_distance = 1.5e-6)

    return pd.concat([halo_stars, thin_disk_stars, thick_disk_stars])


def thetas_cylindrical(*,
                       size: int,
                       angle_covering_sector: float,
                       generator: GeneratorType) -> Union[float, np.ndarray]:
    thetas = generator(low=-angle_covering_sector / 2.,
                       high=angle_covering_sector / 2.,
                       size=size)
    thetas[thetas < 0.] += 2. * np.pi
    return thetas


def halo_star_distance(*,
                       min_distance: float,
                       max_distance: float,
                       core_radius: float) -> sym.Piecewise:
    """
    Equation from:
    Simulating Gaia performances on white dwarfs - by Torres et al. 2005
    """
    x = sym.Symbol('x')
    return sym.Piecewise((0, x < min_distance),
                         (1 / (core_radius ** 2 + x ** 2), x < max_distance),
                         (0, True))


def disk_star_distance(*,
                       min_distance: float,
                       max_distance: float,
                       scale_length: float) -> sym.Piecewise:
    x = sym.Symbol('x')
    return sym.Piecewise((0, x < min_distance),
                         (sym.exp(-x / scale_length), x < max_distance),
                         (0, True))


def halo_z_coordinates(*,
                       angle_covering_sector: float,
                       r_cylindrical: np.ndarray,
                       generator: GeneratorType) -> Union[float, np.ndarray]:
    random_angles = generator(low=-angle_covering_sector / 2.,
                              high=angle_covering_sector / 2.,
                              size=r_cylindrical.size)
    return r_cylindrical * np.sin(random_angles)


def disk_z_coordinate(*,
                      radius: float,
                      scale_height: float) -> sym.Piecewise:
    x = sym.Symbol('x')
    return sym.Piecewise((0, x < -radius),
                         (sym.exp(x / scale_height), x < 0),
                         (sym.exp(-x / scale_height), x < radius),
                         (0, True))
