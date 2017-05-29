import csv
import logging
from math import (ceil,
                  log10,
                  sqrt)
from typing import (List,
                    Iterable)

from alcor.models.star import Star
from alcor.types import StarsBinsType, RowType

logging.basicConfig(format='%(filename)s %(funcName)s '
                           '%(levelname)s: %(message)s',
                    level=logging.DEBUG)
logger = logging.getLogger(__name__)

MIN_BOLOMETRIC_MAGNITUDE = 6.
MAX_BOLOMETRIC_MAGNITUDE = 21.
BIN_SIZE = 0.5
BOLOMETRIC_MAGNITUDE_AMPLITUDE = (MAX_BOLOMETRIC_MAGNITUDE
                                  - MIN_BOLOMETRIC_MAGNITUDE)
BINS_COUNT = int(BOLOMETRIC_MAGNITUDE_AMPLITUDE / BIN_SIZE)

OBSERVATIONAL_DATA_TRUSTED_BINS_OBJECT_COUNT = 220
FORTY_PARSEC_NORTHERN_HEMISPHERE_VOLUME = 134041.29


def write_luminosity_function_data(stars: List[Star]) -> None:
    bins = generate_bins(stars)
    # TODO: comment on the choice of these particular bins
    normalization_factor = (FORTY_PARSEC_NORTHERN_HEMISPHERE_VOLUME
                            * (len(bins[16]) + len(bins[17]) + len(bins[18]))
                            / OBSERVATIONAL_DATA_TRUSTED_BINS_OBJECT_COUNT)

    with open('luminosity_function.csv', 'w') as output_file:
        output_writer = csv.writer(output_file, delimiter=' ')
        output_writer.writerow(['normalization_factor:',
                                normalization_factor])
        output_writer.writerow(['average_bin_magnitude',
                                'star_count_logarithm',
                                'upper_errorbar',
                                'lower_errorbar'])
        for row in rows(bins=bins,
                        normalization_factor=normalization_factor):
            output_writer.writerow(row)


def generate_bins(stars: List[Star]) -> StarsBinsType:
    bins = [[] for _ in range(BINS_COUNT)]
    for star in stars:
        stars_bin = get_magnitude_bin_index(star)
        bins[stars_bin].append(star)
    return bins


def get_magnitude_bin_index(star: Star) -> int:
    return int(ceil((float(star.bolometric_magnitude)
                     - MIN_BOLOMETRIC_MAGNITUDE)
                    / BIN_SIZE))


def rows(bins: StarsBinsType,
         normalization_factor: float) -> Iterable[RowType]:
    non_empty_bins = filter(None, bins)
    for stars_bin_index, stars_bin in enumerate(non_empty_bins):
        stars_count = len(stars_bin)
        average_magnitude = (MIN_BOLOMETRIC_MAGNITUDE
                             + BIN_SIZE * (stars_bin_index - 0.5))
        stars_count_logarithm = get_stars_count_logarithm(
            stars_count=stars_count,
            normalization_factor=normalization_factor)
        upper_errorbar = get_upper_errorbar(stars_count)
        lower_errorbar = get_lower_errorbar(stars_count)
        yield (average_magnitude,
               stars_count_logarithm,
               upper_errorbar,
               lower_errorbar)


def get_stars_count_logarithm(stars_count: int,
                              normalization_factor: float) -> float:
    try:
        return log10(stars_count / normalization_factor)
    except ValueError:
        return 0.


def get_upper_errorbar(stars_count: int) -> float:
    if stars_count == 0:
        return 0.
    else:
        return (log10((stars_count + sqrt(stars_count))
                      / FORTY_PARSEC_NORTHERN_HEMISPHERE_VOLUME)
                - log10(stars_count / FORTY_PARSEC_NORTHERN_HEMISPHERE_VOLUME))


def get_lower_errorbar(stars_count: int) -> float:
    if stars_count == 0:
        return 0.
    elif stars_count == 1:
        return 5.  # Random number so that errorbar would go below the plot
    else:
        return (log10(stars_count / FORTY_PARSEC_NORTHERN_HEMISPHERE_VOLUME)
                - log10((stars_count - sqrt(stars_count))
                        / FORTY_PARSEC_NORTHERN_HEMISPHERE_VOLUME))
