from typing import (Callable,
                    Tuple)

import numpy as np

from alcor.services.simulations.polar import (thetas_cylindrical,
                                              halo_stars_radii_tries,
                                              disks_stars_radii_tries,
                                              halo_r_cylindrical,
                                              disks_r_cylindrical)


def test_theta_cylindrical(
        size: int,
        angle_covering_sector: float,
        generator: Callable[[float, float, float], np.ndarray]) -> None:
    thetas = thetas_cylindrical(
            size=size,
            angle_covering_sector=angle_covering_sector,
            generator=generator)

    assert isinstance(thetas, np.ndarray)
    assert thetas.size == size


def test_halo_stars_radii_tries(
        size: int,
        min_sector_radius: float,
        max_sector_radius: float,
        halo_core_radius: float,
        generator: Callable[[Tuple[int, ...]], np.ndarray]) -> None:
    radii_tries = halo_stars_radii_tries(size=size,
                                         min_sector_radius=min_sector_radius,
                                         max_sector_radius=max_sector_radius,
                                         halo_core_radius=halo_core_radius,
                                         generator=generator)

    assert (isinstance(radii_tries, np.ndarray) or
            isinstance(radii_tries, float))
    assert radii_tries.size == size


def test_disks_stars_radii_tries(
        size: int,
        min_sector_radius: float,
        max_sector_radius: float,
        scale_length: float,
        radial_distrib_max: float) -> None:
    radii_tries = disks_stars_radii_tries(
            size=size,
            min_sector_radius=min_sector_radius,
            max_sector_radius=max_sector_radius,
            scale_length=scale_length,
            radial_distrib_max=radial_distrib_max)

    assert (isinstance(radii_tries, float) or
            isinstance(radii_tries, np.ndarray))
    assert radii_tries.size == size


def test_halo_r_cylindrical(size: int,
                            min_sector_radius: float,
                            max_sector_radius: float,
                            halo_core_radius: float,
                            squared_min_sector_radius: float,
                            squared_radii_difference: float) -> None:
    coordinates = halo_r_cylindrical(
            size=size,
            min_sector_radius=min_sector_radius,
            max_sector_radius=max_sector_radius,
            halo_core_radius=halo_core_radius,
            squared_min_sector_radius=squared_min_sector_radius,
            squared_radii_difference=squared_radii_difference)

    assert (isinstance(coordinates, float) or
            isinstance(coordinates, np.ndarray))
    assert coordinates.size == size


def test_disks_r_cylindrical(size: int,
                             min_sector_radius: float,
                             max_sector_radius: float,
                             scale_length: float,
                             radial_distrib_max: float,
                             squared_min_sector_radius: float,
                             squared_radii_difference: float) -> None:
    coordinates = disks_r_cylindrical(
            size=size,
            min_sector_radius=min_sector_radius,
            max_sector_radius=max_sector_radius,
            scale_length=scale_length,
            radial_distrib_max=radial_distrib_max,
            squared_min_sector_radius=squared_min_sector_radius,
            squared_radii_difference=squared_radii_difference)

    assert (isinstance(coordinates, float) or
            isinstance(coordinates, np.ndarray))
    assert coordinates.size == size
