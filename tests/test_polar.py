from typing import Callable

import numpy as np

from alcor.services.simulations.polar import (thetas_cylindrical,
                                              halo_z_coordinates)


def test_theta_cylindrical(
        size: int,
        angle_covering_sector: float,
        array_generator: Callable[[float, float, float], np.ndarray]) -> None:
    thetas = thetas_cylindrical(
            size=size,
            angle_covering_sector=angle_covering_sector,
            generator=array_generator)

    assert isinstance(thetas, np.ndarray)
    assert thetas.size == size


def test_halo_z_coordinates(
        angle_covering_sector: float,
        r_cylindrical: np.ndarray,
        array_generator: Callable[[float, float, float], np.ndarray]) -> None:
    coordinates = halo_z_coordinates(
            angle_covering_sector=angle_covering_sector,
            r_cylindrical=r_cylindrical,
            generator=array_generator)

    assert isinstance(coordinates, (float, np.ndarray))
    assert coordinates.size == r_cylindrical.size
