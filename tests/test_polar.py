from typing import Callable

import numpy as np

from alcor.services.simulations.polar import (thetas_cylindrical,
                                              get_radii_tries)


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


def test_get_radii_tries(galactic_structures: np.ndarray,
                         min_sector_radius: float,
                         max_sector_radius: float,
                         halo_core_radius: float,
                         sector_diameter: float,
                         scale_length: float,
                         radial_distrib_max: float) -> None:
    radii_tries = get_radii_tries(galactic_structures=galactic_structures,
                                  min_sector_radius=min_sector_radius,
                                  max_sector_radius=max_sector_radius,
                                  halo_core_radius=halo_core_radius,
                                  sector_diameter=sector_diameter,
                                  scale_length=scale_length,
                                  radial_distrib_max=radial_distrib_max)

    assert isinstance(radii_tries, np.ndarray)
    assert radii_tries.size == galactic_structures.size
