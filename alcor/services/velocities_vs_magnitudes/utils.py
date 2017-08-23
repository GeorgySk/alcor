from alcor.models import Star

MIN_BOLOMETRIC_MAGNITUDE = 6.0
MAX_BOLOMETRIC_MAGNITUDE = 30.0
STARS_BIN_SIZE = 0.5


def bolometric_magnitude_to_stars_bin_index(
        magnitude: float,
        *,
        min_magnitude: float = MIN_BOLOMETRIC_MAGNITUDE,
        stars_bin_size: float = STARS_BIN_SIZE) -> int:
    magnitude_amplitude = magnitude - min_magnitude
    return int(magnitude_amplitude / stars_bin_size)


STARS_BINS_COUNT = bolometric_magnitude_to_stars_bin_index(
        MAX_BOLOMETRIC_MAGNITUDE)
DEFAULT_VELOCITY_STD = 100.


def stars_bin_index(star: Star) -> int:
    return bolometric_magnitude_to_stars_bin_index(star.bolometric_magnitude)
