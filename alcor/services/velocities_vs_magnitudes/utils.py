from alcor.services.common import (bolometric_indexer,
                                   star_bolometric_indexer,
                                   stars_packer)

MIN_BOLOMETRIC_MAGNITUDE = 6.0
MAX_BOLOMETRIC_MAGNITUDE = 30.0
STARS_BIN_SIZE = 0.5

bolometric_index = bolometric_indexer(min_magnitude=MIN_BOLOMETRIC_MAGNITUDE,
                                      stars_bin_size=STARS_BIN_SIZE)
STARS_BINS_COUNT = bolometric_index(MAX_BOLOMETRIC_MAGNITUDE)
DEFAULT_VELOCITY_STD = 100.

star_bolometric_index = star_bolometric_indexer(bolometric_index)
pack_stars = stars_packer(stars_bins_count=STARS_BINS_COUNT,
                          star_bolometric_index=star_bolometric_index)
