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


def distribute_into_bins_for_luminosity_function(star: Star,
                                                 bins: List[List[Star]]) -> (
        None):
    stars_bin = int(ceil((star.bolometric_magnitude - MIN_BOLOMETRIC_MAGNITUDE)
                    / BIN_SIZE))
    bins[stars_bin].append(star)


def write_bins_luminosity_function_info(bins: List[List[Star]]) -> None:
    observational_data_trusted_bins_object_count = 220
    forty_parsec_northern_hemisphere_volume = 134041.29
    # TODO: comment on the choice of these particular bins
    normalization_factor = (forty_parsec_northern_hemisphere_volume
                            * (len(bins[16]) + len(bins[17]) + len(bins[18]))
                            / observational_data_trusted_bins_object_count)

    with open('luminosity_function.csv', 'w') as output_file:
        output_writer = csv.writer(output_file, delimiter='  ')
        output_writer.writerow('normalization_factor:',
                               normalization_factor)
        output_writer.writerow('average_bin_magnitude',
                               'star_count_logarithm',
                               'upper_errorbar',
                               'lower_errorbar')

        for stars_bin_index, stars_bin in enumerate(bins):
            normalized_stars_number_in_bin = (len(stars_bin)
                                              / normalization_factor)
            average_bin_magnitude = (MIN_BOLOMETRIC_MAGNITUDE
                                     + BIN_SIZE * (stars_bin_index - 0.5))
            if normalized_stars_number_in_bin == 0:
                star_count_logarithm = 0.0
            else:
                star_count_logarithm = log10(normalized_stars_number_in_bin)
            upper_errorbar = (log10((len(stars_bin) + sqrt(len(stars_bin)))
                                    / forty_parsec_northern_hemisphere_volume)
                              - len(stars_bin))
            lower_errorbar = (
                len(stars_bin)
                - log10((len(stars_bin) - sqrt(len(stars_bin)))
                        / forty_parsec_northern_hemisphere_volume))
            output_writer.writerow(average_bin_magnitude,
                                   star_count_logarithm,
                                   upper_errorbar,
                                   lower_errorbar)
        output_writer.writerow('\n')
