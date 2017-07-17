import logging
from collections import Counter
from functools import partial
from itertools import filterfalse
from typing import List

from sqlalchemy.orm.session import Session

from alcor.models import (Group,
                          Star)
from alcor.models.eliminations import StarsCounter
from alcor.services import (luminosity_function,
                            velocities,
                            velocities_vs_magnitudes)
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
                        w_luminosity_function: bool,
                        w_velocities_clouds: bool,
                        w_velocities_vs_magnitude: bool,
                        w_lepine_criterion: bool,
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

    session.add(counter)

    if nullify_radial_velocity:
        for star in stars:
            star.set_radial_velocity_to_zero()

    if w_luminosity_function:
        luminosity_function.process_stars_group(
            stars=stars,
            group=group,
            session=session)

    if w_velocities_clouds:
        velocities.process_stars_group(
            stars=stars,
            group=group,
            w_lepine_criterion=w_lepine_criterion,
            session=session)

    if w_velocities_vs_magnitude:
        velocities_vs_magnitudes.process_stars_group(
            stars=stars,
            group=group,
            w_lepine_criterion=w_lepine_criterion,
            session=session)
