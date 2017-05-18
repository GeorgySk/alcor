from collections import Counter
from functools import partial
from itertools import filterfalse
from typing import Dict

from .kinematics import (write_data_for_velocity_clouds,
                         write_data_for_velocities_vs_magnitude)
from .luminosity_function import (BINS_COUNT,
                                  distribute_into_bins_for_luminosity_function,
                                  write_bins_luminosity_function_info)
from .sampling import (write_elimination_stats,
                       check_elimination)

from alcor.utils import parse_stars


def run_processing(*,
                   settings: Dict) -> None:
    with open('output.res', 'r') as output_file:
        full_stars_sample = list(parse_stars(output_file, 'test'))

    eliminations_counter = Counter(int)
    apply_elimination_criteria = partial(check_elimination,
                                         eliminations_counter
                                         =eliminations_counter)
    restricted_stars_sample = filterfalse(apply_elimination_criteria,
                                          full_stars_sample)

    bins = [[] for _ in range(BINS_COUNT)]

    for star in restricted_stars_sample:
        star.set_radial_velocity_to_zero()
        distribute_into_bins_for_luminosity_function(star, bins)

    write_bins_luminosity_function_info(bins)
    write_elimination_stats(len(full_stars_sample),
                            len(restricted_stars_sample),
                            eliminations_counter)
    write_data_for_velocity_clouds(restricted_stars_sample)
    write_data_for_velocities_vs_magnitude(restricted_stars_sample)
