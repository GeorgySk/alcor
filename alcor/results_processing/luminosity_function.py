import csv
from math import (ceil,
                  log10,
                  sqrt)
from typing import List

from alcor.models.star import Star

MIN_BOLOMETRIC_MAGNITUDE = 6.0
MAX_BOLOMETRIC_MAGNITUDE = 21.0
BIN_SIZE = 0.5
BOLOMETRIC_MAGNITUDE_AMPLITUDE = (MAX_BOLOMETRIC_MAGNITUDE
                                  - MIN_BOLOMETRIC_MAGNITUDE)
BINS_COUNT = int(BOLOMETRIC_MAGNITUDE_AMPLITUDE / BIN_SIZE)


def distribute_into_bins(star: Star,
                         bins: List[List[Star]]) -> None:
    bolometric_magnitude = 2.5 * star.luminosity + 4.75
    data_bin = int(ceil((bolometric_magnitude - MIN_BOLOMETRIC_MAGNITUDE)
                   / BIN_SIZE))
    bins[data_bin].append(star)


def write_bins_luminosity_function_info(bins: List[List[Star]]) -> None:
    observational_data_trusted_bins_object_count = 220
    forty_parsec_northern_hemisphere_volume = 134041.29
    normalization_factor = (forty_parsec_northern_hemisphere_volume
                            * (len(bins[16]) + len(bins[17]) + len(bins[18]))
                            / observational_data_trusted_bins_object_count)

    with open('luminosity_function.csv', 'w') as output_file:
        output_writer = csv.writer(output_file, delimiter='  ')
        output_writer.writerow(normalization_factor)

        for data_bin_index, data_bin in enumerate(bins):
            normalized_stars_number_in_bin = (len(data_bin)
                                              / normalization_factor)
            average_bin_magnitude = (MIN_BOLOMETRIC_MAGNITUDE
                                     + BIN_SIZE * (data_bin_index - 0.5))
            if normalized_stars_number_in_bin == 0:
                star_count_logarithm = 0.0
            else:
                star_count_logarithm = log10(normalized_stars_number_in_bin)
            upper_errorbar = (log10((len(data_bin) + sqrt(len(data_bin)))
                                    / forty_parsec_northern_hemisphere_volume)
                              - len(data_bin))
            lower_errorbar = (len(data_bin) - log10((len(data_bin)
                                                     - sqrt(len(data_bin)))
                                               / forty_parsec_northern_hemisphere_volume))
            output_writer.writerow(average_bin_magnitude,
                                   star_count_logarithm,
                                   upper_errorbar,
                                   lower_errorbar)
        output_writer.writerow('\n')
