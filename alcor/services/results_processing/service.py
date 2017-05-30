import logging
import uuid

from collections import Counter
from functools import partial
from itertools import filterfalse

from .kinematics import (write_velocity_clouds_data,
                         write_velocities_vs_magnitude_data)
from .luminosity_function import write_luminosity_function_data
from .sampling import (write_elimination_stats,
                       check_elimination)
from alcor.utils import parse_stars

logging.basicConfig(format='%(filename)s %(funcName)s '
                           '%(levelname)s: %(message)s',
                    level=logging.DEBUG)
logger = logging.getLogger(__name__)


def run_processing(*,
                   data_path: str,
                   method: str,
                   nullify_radial_velocity: bool,
                   luminosity_function: bool,
                   velocity_clouds: bool,
                   velocities_vs_magnitude: bool,
                   lepine_criterion: bool) -> None:
    group_id = uuid.uuid4()
    with open(data_path, 'r') as input_file:
        stars = list(parse_stars(input_file, group_id))

    eliminations_counter = Counter()
    apply_elimination_criteria = partial(
        check_elimination,
        eliminations_counter=eliminations_counter,
        method=method)
    if method in {'restricted', 'full'}:
        filtered_stars = list(filterfalse(apply_elimination_criteria,
                                          stars))
    else:
        filtered_stars = stars
    raw_sample_stars_count = len(list(stars))
    filtered_stars_count = len(filtered_stars)
    write_elimination_stats(raw_sample_stars_count=raw_sample_stars_count,
                            filtered_stars_count=filtered_stars_count,
                            eliminations_counter=eliminations_counter)

    if nullify_radial_velocity:
        for star in filtered_stars:
            star.set_radial_velocity_to_zero()

    if luminosity_function:
        write_luminosity_function_data(filtered_stars)

    if velocity_clouds:
        write_velocity_clouds_data(filtered_stars,
                                   lepine_criterion)

    if velocities_vs_magnitude:
        write_velocities_vs_magnitude_data(filtered_stars,
                                           lepine_criterion)
