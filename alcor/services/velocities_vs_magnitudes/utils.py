from alcor.models import Star

MIN_BOLOMETRIC_MAGNITUDE = 6.0
MAX_BOLOMETRIC_MAGNITUDE = 30.0
BIN_SIZE = 0.5
BOLOMETRIC_MAGNITUDE_AMPLITUDE = (MAX_BOLOMETRIC_MAGNITUDE
                                  - MIN_BOLOMETRIC_MAGNITUDE)
BINS_COUNT = int(BOLOMETRIC_MAGNITUDE_AMPLITUDE / BIN_SIZE)
DEFAULT_VELOCITY_STD = 100.


def get_stars_bin_index(star: Star) -> int:
    return int((float(star.bolometric_magnitude)
                - MIN_BOLOMETRIC_MAGNITUDE) / BIN_SIZE)
