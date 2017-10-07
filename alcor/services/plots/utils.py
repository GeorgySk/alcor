from functools import partial
from math import radians
from typing import Callable, Tuple

# TODO: sort out the mess with numpy and/or pandas
import numpy as np
import pandas as pd

SOLAR_ABSOLUTE_BOLOMETRIC_MAGNITUDE = 4.75
DEC_GPOLE = radians(27.128336)
RA_GPOLE = radians(192.859508)
AUX_ANGLE = radians(122.932)

nan_array = partial(np.full,
                    fill_value=np.nan)


def bolometric_indexer(*,
                       min_magnitude: float,
                       stars_bin_size: float) -> Callable[[np.ndarray],
                                                          np.ndarray]:
    def bolometric_index(magnitudes: np.ndarray) -> np.ndarray:
        magnitude_amplitudes = magnitudes - min_magnitude
        return np.floor(magnitude_amplitudes / stars_bin_size).astype(np.int32)

    return bolometric_index


def bolometric_magnitude(luminosities: pd.Series) -> pd.Series:
    # More info at
    # https://en.wikipedia.org/wiki/Absolute_magnitude#Bolometric_magnitude
    return 2.5 * luminosities + SOLAR_ABSOLUTE_BOLOMETRIC_MAGNITUDE


def to_cartesian_from_equatorial(stars: pd.DataFrame) -> Tuple[pd.Series,
                                                               pd.Series,
                                                               pd.Series]:
    latitudes = (np.arcsin(np.cos(stars['declination']) * np.cos(DEC_GPOLE)
                           * np.cos(stars['right_ascension'] - RA_GPOLE)
                           + np.sin(stars['declination']) * np.sin(DEC_GPOLE)))
    x = np.sin(stars['declination']) - np.sin(latitudes) * np.sin(DEC_GPOLE)
    y = (np.cos(stars['declination'])
         * np.sin(stars['right_ascension'] - RA_GPOLE) * np.cos(DEC_GPOLE))
    longitudes = np.arctan(x / y) + AUX_ANGLE - np.pi / 2.
    longitudes[((x > 0.) & (y < 0.)) | ((x <= 0.) & (y <= 0.))] += np.pi

    x_coordinates = stars['distance'] * np.cos(latitudes) * np.cos(longitudes)
    y_coordinates = stars['distance'] * np.cos(latitudes) * np.sin(longitudes)
    z_coordinates = stars['distance'] * np.sin(latitudes)
    return x_coordinates, y_coordinates, z_coordinates
