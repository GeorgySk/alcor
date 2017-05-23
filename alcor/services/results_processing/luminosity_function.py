import csv
import logging
from math import (ceil,
                  log10,
                  sqrt)
from typing import (List,
                    Iterable)

from decimal import Decimal

from alcor.models.star import Star
from alcor.types import StarsBinsType, RowType

logging.basicConfig(format='%(filename)s %(funcName)s '
                           '%(levelname)s: %(message)s',
                    level=logging.DEBUG)
logger = logging.getLogger(__name__)

MIN_BOLOMETRIC_MAGNITUDE = Decimal(6.0)
MAX_BOLOMETRIC_MAGNITUDE = Decimal(21.0)
BIN_SIZE = Decimal(0.5)
BOLOMETRIC_MAGNITUDE_AMPLITUDE = (MAX_BOLOMETRIC_MAGNITUDE
                                  - MIN_BOLOMETRIC_MAGNITUDE)
BINS_COUNT = int(BOLOMETRIC_MAGNITUDE_AMPLITUDE / BIN_SIZE)

OBSERVATIONAL_DATA_TRUSTED_BINS_OBJECT_COUNT = Decimal(220)
FORTY_PARSEC_NORTHERN_HEMISPHERE_VOLUME = Decimal(134041.29)


def write_luminosity_function_data(stars: List[Star]) -> None:
    bins = apply_binning_distribution(stars)
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


def apply_binning_distribution(stars: List[Star]) -> StarsBinsType:
    bins = [[] for _ in range(BINS_COUNT)]
    for star in stars:
        stars_bin = get_magnitude_bin_index(star)
        bins[stars_bin].append(star)
    return bins


def get_magnitude_bin_index(star: Star) -> int:
    return int(ceil((star.bolometric_magnitude
                     - MIN_BOLOMETRIC_MAGNITUDE)
                    / BIN_SIZE))


def rows(bins: StarsBinsType,
         normalization_factor: float) -> Iterable[RowType]:
    for stars_bin_index, stars_bin in enumerate(bins):
        average_magnitude = (MIN_BOLOMETRIC_MAGNITUDE
                             + BIN_SIZE * Decimal(stars_bin_index - 0.5))
        stars_count_logarithm = get_stars_count_logarithm(
            stars_count=len(stars_bin),
            normalization_factor=normalization_factor)
        upper_errorbar = (log10(Decimal(len(stars_bin) + sqrt(len(stars_bin)))
                                / FORTY_PARSEC_NORTHERN_HEMISPHERE_VOLUME)
                          - len(stars_bin))
        lower_errorbar = (len(stars_bin)
                          - log10(Decimal(len(stars_bin)
                                          - sqrt(len(stars_bin)))
                                  / FORTY_PARSEC_NORTHERN_HEMISPHERE_VOLUME))
        yield(average_magnitude,
              stars_count_logarithm,
              upper_errorbar,
              lower_errorbar)


def get_stars_count_logarithm(stars_count: int,
                              normalization_factor: float) -> float:
    try:
        return log10(stars_count / normalization_factor)
    except ValueError:
        return 0.
