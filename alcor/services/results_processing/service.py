from collections import Counter
from functools import partial
from itertools import filterfalse
from typing import Dict

from alcor.services.results_processing.kinematics import (
    write_velocity_clouds_data,
    write_velocities_vs_magnitude_data)
from alcor.services.results_processing.luminosity_function import (
    write_luminosity_function_data)
from .sampling import (write_elimination_stats,
                       check_elimination)
from alcor.utils import parse_stars


def run_processing(*, settings: Dict[str, str]) -> None:
    input_path = settings['input_path']
    with open(input_path, 'r') as input_file:
        raw_stars_sample = list(parse_stars(input_file, ''))

    methods = settings['methods']
    plots = settings['plots']

    sampling_method = methods['sampling_method']
    eliminations_counter = Counter(int)
    apply_elimination_criteria = partial(check_elimination,
                                         eliminations_counter
                                         =eliminations_counter,
                                         method=sampling_method)
    if sampling_method == 'full':
        stars = filterfalse(apply_elimination_criteria,
                            raw_stars_sample)
    elif sampling_method == 'restricted':
        stars = filterfalse(apply_elimination_criteria,
                            raw_stars_sample)
    write_elimination_stats(len(raw_stars_sample),
                            len(stars),
                            eliminations_counter)

    if methods['nullify_radial_velocity']:
        for star in stars:
            star.set_radial_velocity_to_zero()

    lepine_criterion_applied = methods['lepine_criterion_applied']

    if plots['luminosity_function']:
        write_luminosity_function_data(stars)
    if plots['velocity_clouds']:
        write_velocity_clouds_data(stars,
                                   lepine_criterion_applied)
    if plots['velocities_vs_magnitude']:
        write_velocities_vs_magnitude_data(stars,
                                           lepine_criterion_applied)
