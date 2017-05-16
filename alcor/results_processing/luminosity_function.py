from math import ceil

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
    bin_count = int(ceil((bolometric_magnitude - MIN_BOLOMETRIC_MAGNITUDE)
                         / BIN_SIZE))
    bins[bin_count].append(star)
