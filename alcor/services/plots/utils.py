from functools import partial
from typing import Callable

# TODO: sort out the mess with numpy and/or pandas
import numpy as np
import pandas as pd

SOLAR_ABSOLUTE_BOLOMETRIC_MAGNITUDE = 4.75

nan_array = partial(np.full,
                    fill_value=np.nan)


def bolometric_indexer(*,
                       min_magnitude: float,
                       stars_bin_size: float) -> Callable[[np.ndarray],
                                                          np.ndarray]:
    def bolometric_index(magnitude: np.ndarray) -> np.ndarray:
        magnitude_amplitude = magnitude - min_magnitude
        return np.floor(magnitude_amplitude / stars_bin_size).astype(np.int32)

    return bolometric_index


def bolometric_magnitude(luminosity: pd.Series) -> pd.Series:
    # More info at
    # https://en.wikipedia.org/wiki/Absolute_magnitude#Bolometric_magnitude
    return 2.5 * luminosity + SOLAR_ABSOLUTE_BOLOMETRIC_MAGNITUDE