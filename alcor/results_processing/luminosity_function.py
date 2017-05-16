from math import ceil

from numpy.ma.core import std

from alcor.models.star import Star

MIN_BOLOMETRIC_MAGNITUDE = 6.0
MAX_BOLOMETRIC_MAGNITUDE = 21.0
BIN_SIZE = 0.5
BOLOMETRIC_MAGNITUDE_AMPLITUDE = (MAX_BOLOMETRIC_MAGNITUDE
                                  - MIN_BOLOMETRIC_MAGNITUDE)
BINS_COUNT = BOLOMETRIC_MAGNITUDE_AMPLITUDE / BIN_SIZE


def distribute_into_bins(star: Star,
                         bins: list[list[Star]]) -> None:
    bolometric_magnitude = 2.5 * star.luminosity + 4.75
    bin = int(ceil((bolometric_magnitude - MIN_BOLOMETRIC_MAGNITUDE)
                   / BIN_SIZE))
    bins[bin].append(star)

def write_bins_info(bins: list[list[Star]]) -> None:
    with open('magnitude_bins.res', 'w') as output_file:
        for bin_index, bin in enumerate(bins):
            average_bin_magnitude = (MIN_BOLOMETRIC_MAGNITUDE
                                     + BIN_SIZE * (bin_index - 0.5))
            average_bin_velocity_u = sum(bin(:).velocity_u) / len(bin)
            average_bin_velocity_v = sum(bin(:).velocity_v) / len(bin)
            average_bin_velocity_w = sum(bin(:).velocity_w) / len(bin)
            bin_standard_deviation_of_velocity_u = std.(bin(:).velocity_u)
            bin_standard_deviation_of_velocity_v = std.(bin(:).velocity_v)
            bin_standard_deviation_of_velocity_w = std.(bin(:).velocity_w)
            output_file.write(average_bin_magnitude,
                              average_bin_velocity_u,
                              average_bin_velocity_v,
                              average_bin_velocity_w,
                              bin_standard_deviation_of_velocity_u,
                              bin_standard_deviation_of_velocity_v,
                              bin_standard_deviation_of_velocity_w)
            