import logging
from collections import Counter
from functools import partial
from itertools import filterfalse
from typing import List

from cassandra.cluster import Session

from alcor.models import (Group,
                          Star)
from alcor.models.eliminations import StarsCounter
from alcor.services.data_access import (insert,
                                        model_insert_statement)
from alcor.services.luminosity_function import (
    process_stars_group_luminosity_function)
from alcor.services.velocities import process_stars_group_velocities_clouds
from alcor.services.velocities_vs_magnitudes import (
    process_stars_group_velocities_vs_magnitudes)
from .sampling import check_elimination

logging.basicConfig(format='%(filename)s %(funcName)s '
                           '%(levelname)s: %(message)s',
                    level=logging.DEBUG)
logger = logging.getLogger(__name__)


def process_stars_group(*,
                        stars: List[Star],
                        group: Group,
                        filtration_method: str,
                        nullify_radial_velocity: bool,
                        luminosity_function: bool,
                        velocities_clouds: bool,
                        velocities_vs_magnitude: bool,
                        lepine_criterion: bool,
                        session: Session) -> None:
    stars_count = len(stars)
    logger.info('Starting processing stars, '
                f'objects number: {stars_count}.')
    eliminations_counter = Counter()
    apply_elimination_criteria = partial(
        check_elimination,
        eliminations_counter=eliminations_counter,
        filtration_method=filtration_method)
    if filtration_method in {'restricted', 'full'}:
        stars = list(filterfalse(apply_elimination_criteria,
                                 stars))

    counter = StarsCounter(
        group_id=group.id,
        raw=stars_count,
        by_parallax=eliminations_counter['parallax'],
        by_declination=eliminations_counter['declination'],
        by_velocity=eliminations_counter['velocity'],
        by_proper_motion=eliminations_counter['proper_motion'],
        by_reduced_proper_motion=eliminations_counter['reduced_proper_motion'],
        by_apparent_magnitude=eliminations_counter['apparent_magnitude'])
    statement = model_insert_statement(StarsCounter)
    insert(instances=[counter],
           statement=statement,
           session=session)

    if nullify_radial_velocity:
        for star in stars:
            star.set_radial_velocity_to_zero()

    if luminosity_function:
        process_stars_group_luminosity_function(
            stars=stars,
            group=group,
            session=session)

    if velocities_clouds:
        process_stars_group_velocities_clouds(
            stars=stars,
            group=group,
            lepine_criterion=lepine_criterion,
            session=session)

    if velocities_vs_magnitude:
        process_stars_group_velocities_vs_magnitudes(
            stars=stars,
            group=group,
            lepine_criterion=lepine_criterion,
            session=session)
