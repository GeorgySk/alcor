from math import (log10,
                  sqrt)
from typing import (Iterable,
                    Iterator)

from alcor.models import Group
from alcor.models.luminosity_function import Point
from alcor.services.common import (FORTY_PARSEC_NORTHERN_HEMISPHERE_VOLUME,
                                   bolometric_indexer,
                                   star_bolometric_indexer,
                                   stars_packer)
from alcor.types import StarsBinsType

MIN_BOLOMETRIC_MAGNITUDE = 6.
MAX_BOLOMETRIC_MAGNITUDE = 21.
STARS_BIN_SIZE = 0.5
bolometric_index = bolometric_indexer(min_magnitude=MIN_BOLOMETRIC_MAGNITUDE,
                                      stars_bin_size=STARS_BIN_SIZE)

STARS_BINS_COUNT = bolometric_index(MAX_BOLOMETRIC_MAGNITUDE)

# More info at (The white dwarf population within 40 pc of the Sun,
# Torres et al., 2016)
SMALLEST_ERROR_BARS_BINS_INDEXES = {15, 16, 17}

OBSERVATIONAL_DATA_TRUSTED_BINS_OBJECT_COUNT = 220

star_bolometric_index = star_bolometric_indexer(bolometric_index)
pack_stars = stars_packer(stars_bins_count=STARS_BINS_COUNT,
                          star_bolometric_index=star_bolometric_index)


def stars_bins_normalization_factor(
        stars_bins: StarsBinsType,
        stars_bins_indexes: Iterable[int] = SMALLEST_ERROR_BARS_BINS_INDEXES
) -> float:
    smallest_error_bars_bins_lengths_sum = (
        sum(len(stars_bins[index])
            for index in stars_bins_indexes))
    return (FORTY_PARSEC_NORTHERN_HEMISPHERE_VOLUME
            * smallest_error_bars_bins_lengths_sum
            / OBSERVATIONAL_DATA_TRUSTED_BINS_OBJECT_COUNT)


def graph_points(*,
                 stars_counts: Iterable[int],
                 group: Group,
                 normalization_factor: float) -> Iterator[Point]:
    group_id = group.id
    for index, stars_count in enumerate(stars_counts):
        if not stars_count:
            continue

        avg_bin_magnitude = (MIN_BOLOMETRIC_MAGNITUDE
                             + STARS_BIN_SIZE * (index + 0.5))
        logarithm = to_logarithm(stars_count=stars_count,
                                 normalization_factor=normalization_factor)
        upper_error_bar = to_upper_error_bar(stars_count)
        lower_error_bar = to_lower_error_bar(stars_count)
        yield Point(group_id=group_id,
                    avg_bin_magnitude=avg_bin_magnitude,
                    stars_count_logarithm=logarithm,
                    upper_error_bar=upper_error_bar,
                    lower_error_bar=lower_error_bar)


def to_logarithm(stars_count: int,
                 normalization_factor: float) -> float:
    try:
        return log10(stars_count / normalization_factor)
    except ValueError:
        return 0.


def to_upper_error_bar(stars_count: int) -> float:
    if stars_count == 0:
        return 0.
    else:
        return (log10((stars_count + sqrt(stars_count))
                      / FORTY_PARSEC_NORTHERN_HEMISPHERE_VOLUME)
                - log10(stars_count / FORTY_PARSEC_NORTHERN_HEMISPHERE_VOLUME))


def to_lower_error_bar(stars_count: int) -> float:
    if stars_count == 0:
        return 0.
    elif stars_count == 1:
        return 5.  # Some number so that error bar would go below the plot
    else:
        return (log10(stars_count / FORTY_PARSEC_NORTHERN_HEMISPHERE_VOLUME)
                - log10((stars_count - sqrt(stars_count))
                        / FORTY_PARSEC_NORTHERN_HEMISPHERE_VOLUME))
