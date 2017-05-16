from collections import Counter
from functools import partial
from itertools import filterfalse

from .kinematics import (write_bins_kinematic_info,
                         write_data_for_velocity_clouds)
from .luminosity_function import (BINS_COUNT,
                                  distribute_into_bins,
                                  write_bins_luminosity_function_info)
from .sampling import (write_elimination_stats,
                       check_elimination)

from alcor.utils import parse_stars


def main() -> None:
    with open('output.res', 'r') as output_file:
        full_stars_sample = list(parse_stars(output_file, 'test'))

    elimination_counters = Counter(int)
    apply_elimination_criteria = partial(check_elimination,
                                         elimination_counters
                                         =elimination_counters)
    restricted_stars_sample = filterfalse(apply_elimination_criteria,
                                          full_stars_sample)

    bins = [] * BINS_COUNT

    for star in restricted_stars_sample:
        star.set_radial_velocity_to_zero()
        distribute_into_bins(star, bins)

    write_bins_kinematic_info(bins)
    write_bins_luminosity_function_info(bins)
    write_elimination_stats(len(full_stars_sample),
                            len(restricted_stars_sample),
                            elimination_counters)
    write_data_for_velocity_clouds(restricted_stars_sample)


if __name__ == '__main__':
    main()
