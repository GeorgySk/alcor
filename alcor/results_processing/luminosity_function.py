from math import (ceil,
                  log10,
                  sqrt)

from numpy.ma.core import std

from alcor.models.star import Star

MIN_BOLOMETRIC_MAGNITUDE = 6.0
MAX_BOLOMETRIC_MAGNITUDE = 21.0
BIN_SIZE = 0.5
BOLOMETRIC_MAGNITUDE_AMPLITUDE = (MAX_BOLOMETRIC_MAGNITUDE
                                  - MIN_BOLOMETRIC_MAGNITUDE)
BINS_COUNT = int(BOLOMETRIC_MAGNITUDE_AMPLITUDE / BIN_SIZE)


def distribute_into_bins(star: Star,
                         bins: list[list[Star]]) -> None:
    bolometric_magnitude = 2.5 * star.luminosity + 4.75
    bin = int(ceil((bolometric_magnitude - MIN_BOLOMETRIC_MAGNITUDE)
                   / BIN_SIZE))
    bins[bin].append(star)


def write_bins_luminosity_function_info(bins: list[list[Star]]) -> None:
    observational_data_trusted_bins_object_count = 220
    forty_parsec_northern_hemisphere_volume = 134041.29
    normalization_factor = (forty_parsec_northern_hemisphere_volume
                            * (len(bins[16]) + len(bins[17]) + len(bins[18]))
                            / observational_data_trusted_bins_object_count)

    with open('luminosity_function.res', 'w') as output_file:
        output_file.write(normalization_factor)

        for bin_index, bin in enumerate(bins):
            normalized_stars_number_in_bin = (len(bin) / normalization_factor)
            average_bin_magnitude = (MIN_BOLOMETRIC_MAGNITUDE
                                     + BIN_SIZE * (bin_index - 0.5))
            if normalized_stars_number_in_bin == 0:
                star_count_logarithm = 0.0
            else:
                star_count_logarithm = log10(normalized_stars_number_in_bin)
            upper_errorbar = (log10((len(bin) + sqrt(len(bin)))
                                    / forty_parsec_northern_hemisphere_volume)
                              - len(bin))
            lower_errorbar = (len(bin) - log10((len(bin) - sqrt(len(bin)))
                                               / forty_parsec_northern_hemisphere_volume))
            output_file.write(average_bin_magnitude,
                              star_count_logarithm,
                              upper_errorbar,
                              lower_errorbar)
