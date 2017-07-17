from math import (ceil,
                  log10,
                  sqrt)
from typing import (List,
                    Iterable)

from alcor.cassandra_models import (CGroup,
                                    CStar)
from alcor.cassandra_models.luminosity_function import CPoint
from alcor.models import Group
from alcor.models.luminosity_function import Point
from alcor.types import StarsBinsType

MIN_BOLOMETRIC_MAGNITUDE = 6.
MAX_BOLOMETRIC_MAGNITUDE = 21.
BIN_SIZE = 0.5
BOLOMETRIC_MAGNITUDE_AMPLITUDE = (MAX_BOLOMETRIC_MAGNITUDE
                                  - MIN_BOLOMETRIC_MAGNITUDE)
BINS_COUNT = int(BOLOMETRIC_MAGNITUDE_AMPLITUDE / BIN_SIZE)

OBSERVATIONAL_DATA_TRUSTED_BINS_OBJECT_COUNT = 220
FORTY_PARSEC_NORTHERN_HEMISPHERE_VOLUME = 134041.29


def generate_stars_bins(stars: List[CStar]) -> StarsBinsType:
    stars_bins = [[] for _ in range(BINS_COUNT)]
    for star in stars:
        index = get_stars_bin_index(star)
        stars_bins[index].append(star)
    return stars_bins


def get_stars_bin_index(star: CStar) -> int:
    return int(ceil((float(star.bolometric_magnitude)
                     - MIN_BOLOMETRIC_MAGNITUDE)
                    / BIN_SIZE))


def points(*,
           stars_bins: StarsBinsType,
           group: CGroup,
           normalization_factor: float,
           cls) -> Iterable[CPoint]:
    non_empty_bins = filter(None, stars_bins)
    for stars_bin_index, stars_bin in enumerate(non_empty_bins):
        stars_count = len(stars_bin)
        avg_bin_magnitude = (MIN_BOLOMETRIC_MAGNITUDE
                             + BIN_SIZE * (stars_bin_index - 0.5))
        stars_count_logarithm = get_stars_count_logarithm(
            stars_count=stars_count,
            normalization_factor=normalization_factor)
        upper_error_bar = get_upper_error_bar(stars_count)
        lower_error_bar = get_lower_error_bar(stars_count)
        yield cls(group_id=group.id,
                  avg_bin_magnitude=avg_bin_magnitude,
                  stars_count_logarithm=stars_count_logarithm,
                  upper_error_bar=upper_error_bar,
                  lower_error_bar=lower_error_bar)


def get_stars_count_logarithm(stars_count: int,
                              normalization_factor: float) -> float:
    try:
        return log10(stars_count / normalization_factor)
    except ValueError:
        return 0.


def get_upper_error_bar(stars_count: int) -> float:
    if stars_count == 0:
        return 0.
    else:
        return (log10((stars_count + sqrt(stars_count))
                      / FORTY_PARSEC_NORTHERN_HEMISPHERE_VOLUME)
                - log10(stars_count / FORTY_PARSEC_NORTHERN_HEMISPHERE_VOLUME))


def get_lower_error_bar(stars_count: int) -> float:
    if stars_count == 0:
        return 0.
    elif stars_count == 1:
        return 5.  # Random number so that error bar would go below the plot
    else:
        return (log10(stars_count / FORTY_PARSEC_NORTHERN_HEMISPHERE_VOLUME)
                - log10((stars_count - sqrt(stars_count))
                        / FORTY_PARSEC_NORTHERN_HEMISPHERE_VOLUME))
