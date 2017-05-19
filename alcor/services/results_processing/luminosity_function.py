import csv
from math import (ceil,
                  log10,
                  sqrt)
from typing import List, Iterable

from alcor.models.star import Star
from alcor.types import StarsBinsType, RowType

MIN_BOLOMETRIC_MAGNITUDE = 6.0
MAX_BOLOMETRIC_MAGNITUDE = 21.0
BIN_SIZE = 0.5
BOLOMETRIC_MAGNITUDE_AMPLITUDE = (MAX_BOLOMETRIC_MAGNITUDE
                                  - MIN_BOLOMETRIC_MAGNITUDE)
BINS_COUNT = BOLOMETRIC_MAGNITUDE_AMPLITUDE // BIN_SIZE

OBSERVATIONAL_DATA_TRUSTED_BINS_OBJECT_COUNT = 220
FORTY_PARSEC_NORTHERN_HEMISPHERE_VOLUME = 134041.29


def distribute_into_bins_for_luminosity_function(star: Star,
                                                 bins: StarsBinsType
                                                 ) -> None:
    stars_bin = get_magnitude_bin_index(star)
    bins[stars_bin].append(star)


def get_magnitude_bin_index(star: Star) -> int:
    return int(ceil((star.bolometric_magnitude
                     - MIN_BOLOMETRIC_MAGNITUDE)
                    / BIN_SIZE))


def write_bins_luminosity_function_info(bins: List[List[Star]]) -> None:
    # TODO: comment on the choice of these particular bins
    normalization_factor = (FORTY_PARSEC_NORTHERN_HEMISPHERE_VOLUME
                            * (len(bins[16]) + len(bins[17]) + len(bins[18]))
                            / OBSERVATIONAL_DATA_TRUSTED_BINS_OBJECT_COUNT)

    with open('luminosity_function.csv', 'w') as output_file:
        output_writer = csv.writer(output_file, delimiter='  ')
        output_writer.writerow('normalization_factor:',
                               normalization_factor)
        output_writer.writerow('average_bin_magnitude',
                               'star_count_logarithm',
                               'upper_errorbar',
                               'lower_errorbar')
        for row in rows(bins):
            output_writer.writerow(row)


def rows(bins: StarsBinsType,
         normalization_factor: float) -> Iterable[RowType]:
    for stars_bin_index, stars_bin in enumerate(bins):
        average_bin_magnitude = (MIN_BOLOMETRIC_MAGNITUDE
                                 + BIN_SIZE * (stars_bin_index - 0.5))
        stars_count_logarithm = get_stars_count_logarithm(
            stars_count=len(stars_bin),
            normalization_factor=normalization_factor)
        upper_errorbar = (log10((len(stars_bin) + sqrt(len(stars_bin)))
                                / FORTY_PARSEC_NORTHERN_HEMISPHERE_VOLUME)
                          - len(stars_bin))
        lower_errorbar = (len(stars_bin)
                          - log10((len(stars_bin) -
                                   sqrt(len(stars_bin)))
                                  / FORTY_PARSEC_NORTHERN_HEMISPHERE_VOLUME))
        yield(average_bin_magnitude,
              stars_count_logarithm,
              upper_errorbar,
              lower_errorbar)


def get_stars_count_logarithm(stars_count: int,
                              normalization_factor: float) -> float:
    try:
        return log10(stars_count / normalization_factor)
    except ValueError:
        return 0.
