from typing import Callable

import numpy as np
import pandas as pd

from alcor.services.simulations.polar import set_thetas_cylindrical


def test_set_theta_cylindrical(
        stars_without_theta: pd.DataFrame,
        angle_covering_sector: float,
        generator: Callable[[float, float, float], np.ndarray]) -> None:
    thetas = set_thetas_cylindrical(
            stars_without_theta,
            angle_covering_sector=angle_covering_sector,
            generator=generator)

    assert isinstance(stars_without_theta, pd.DataFrame)
    assert thetas.size == stars_without_theta.shape[0]
