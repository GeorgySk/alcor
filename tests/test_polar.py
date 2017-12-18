from typing import Callable

import numpy as np

from alcor.services.simulations.polar import thetas_cylindrical


def test_theta_cylindrical(
        size: int,
        angle_covering_sector: float,
        generator: Callable[[float, float, float], np.ndarray]) -> None:
    thetas = thetas_cylindrical(
            size=size,
            angle_covering_sector=angle_covering_sector,
            generator=generator)

    assert thetas.size == size
