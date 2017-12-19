from typing import (Callable,
                    Tuple)

import numpy as np

from alcor.services.simulations.polar import (thetas_cylindrical,
                                              halo_stars_radii_tries,
                                              disks_stars_radii_tries,
                                              halo_r_cylindrical,
                                              disks_r_cylindrical,
                                              halo_z_coordinates,
                                              disk_z_coordinates)


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

    assert (isinstance(radii_tries, float) or
            isinstance(radii_tries, np.ndarray))
    assert radii_tries.size == size


def test_disks_stars_radii_tries(
        size: int,
        min_sector_radius: float,
        max_sector_radius: float,
        scale_length: float,
        radial_distrib_max: float,
        generator: Callable[[float, float, float], np.ndarray]) -> None:
    radii_tries = disks_stars_radii_tries(
            size=size,
            min_sector_radius=min_sector_radius,
            max_sector_radius=max_sector_radius,
            scale_length=scale_length,
            radial_distrib_max=radial_distrib_max,
            generator=generator)

    assert (isinstance(radii_tries, float) or
            isinstance(radii_tries, np.ndarray))
    assert radii_tries.size == size


def test_halo_r_cylindrical(
        size: int,
        min_sector_radius: float,
        max_sector_radius: float,
        halo_core_radius: float,
        squared_min_sector_radius: float,
        squared_radii_difference: float,
        unit_range_generator: Callable[[Tuple[int, ...]], np.ndarray]) -> None:
    coordinates = halo_r_cylindrical(
            size=size,
            min_sector_radius=min_sector_radius,
            max_sector_radius=max_sector_radius,
            halo_core_radius=halo_core_radius,
            squared_min_sector_radius=squared_min_sector_radius,
            squared_radii_difference=squared_radii_difference,
            generator=unit_range_generator)

    assert (isinstance(coordinates, float) or
            isinstance(coordinates, np.ndarray))
    assert coordinates.size == size


def test_disks_r_cylindrical(
        size: int,
        min_sector_radius: float,
        max_sector_radius: float,
        scale_length: float,
        radial_distrib_max: float,
        squared_min_sector_radius: float,
        squared_radii_difference: float,
        generator: Callable[[float, float, float], np.ndarray]) -> None:
    coordinates = disks_r_cylindrical(
            size=size,
            min_sector_radius=min_sector_radius,
            max_sector_radius=max_sector_radius,
            scale_length=scale_length,
            radial_distrib_max=radial_distrib_max,
            squared_min_sector_radius=squared_min_sector_radius,
            squared_radii_difference=squared_radii_difference,
            generator=generator)

    assert (isinstance(coordinates, float) or
            isinstance(coordinates, np.ndarray))
    assert coordinates.size == size


def test_halo_z_coordinates(
        angle_covering_sector: float,
        r_cylindrical: np.ndarray,
        generator: Callable[[float, float, float], np.ndarray]) -> None:
    coordinates = halo_z_coordinates(
            angle_covering_sector=angle_covering_sector,
            r_cylindrical=r_cylindrical,
            generator=generator)

    assert (isinstance(coordinates, float) or
            isinstance(coordinates, np.ndarray))
    assert coordinates.size == r_cylindrical.size


def test_disk_z_coordinates(
        size: int,
        scale_height: float,
        sector_radius: float,
        unit_range_generator: Callable[[Tuple[int, ...]], np.ndarray],
        sign_generator: Callable[[int], np.ndarray]) -> None:
    coordinates = disk_z_coordinates(
            size=size,
            scale_height=scale_height,
            sector_radius=sector_radius,
            generator=unit_range_generator,
            sign_generator=sign_generator)

    assert (isinstance(coordinates, float) or
            isinstance(coordinates, np.ndarray))
    assert coordinates.size == size
