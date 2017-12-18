from typing import Callable

import numpy as np
import pandas as pd

from alcor.services.simulations.polar import set_thetas_cylindrical


def test_set_theta_cylindrical(
        stars_without_theta: pd.DataFrame,
        angle_covering_sector: float,
        generator: Callable[[float, float, float], np.ndarray]) -> None:
    columns_before = stars_without_theta.columns
    set_thetas_cylindrical(stars_without_theta,
                           angle_covering_sector=angle_covering_sector,
                           generator=generator)
    columns_after = stars_without_theta.columns

    assert isinstance(stars_without_theta, pd.DataFrame)
    assert 'theta_cylindrical' not in columns_before
    assert 'theta_cylindrical' in columns_after
