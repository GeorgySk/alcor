import logging
import uuid
from collections import Counter
from functools import partial
from itertools import filterfalse

from sqlalchemy.orm.session import Session

from alcor.models import (Group,
                          Star)
from alcor.models.eliminations import StarsCounter
from alcor.models.processed_star_association import ProcessedStarsAssociation
from alcor.services import (luminosity_function,
                            velocities,
                            velocities_vs_magnitudes)
from alcor.services.data_access import fetch_group_stars
from . import elimination
from .utils import copy_velocities

logger = logging.getLogger(__name__)


def process_stars_group(*,
                        group: Group,
                        filtration_method: str,
                        nullify_radial_velocity: bool,
                        w_luminosity_function: bool,
                        w_velocities_clouds: bool,
                        w_velocities_vs_magnitude: bool,
                        w_lepine_criterion: bool,
                        session: Session) -> None:
    stars = fetch_group_stars(group_id=group.id,
                              session=session)
    stars_count = len(stars)
    logger.info('Starting processing stars, '
                f'objects number: {stars_count}.')

    eliminations_counter = Counter()
    if filtration_method in {'restricted', 'full'}:
        is_eliminated = partial(elimination.check,
                                eliminations_counter=eliminations_counter,
                                filtration_method=filtration_method)
        stars = list(filterfalse(is_eliminated, stars))

    counter = StarsCounter(
            group_id=group.id,
            raw=stars_count,
            by_parallax=eliminations_counter['parallax'],
            by_declination=eliminations_counter['declination'],
            by_velocity=eliminations_counter['velocity'],
            by_proper_motion=eliminations_counter['proper_motion'],
            by_reduced_proper_motion=eliminations_counter[
                'reduced_proper_motion'],
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

    original_id = group.id
    processed_group_id = uuid.uuid4()
    processed_group = Group(
            id=processed_group_id,
            original_id=original_id)
    session.add(processed_group)

    processed_stars = [Star(group_id=processed_group_id)
                       for _ in range(len(stars))]
    if nullify_radial_velocity:
        copy_velocities(stars, processed_stars)
    session.add_all(processed_stars)

    session.commit()

    for star, processed_star in zip(stars, processed_stars):
        processed_star_association = ProcessedStarsAssociation(
                original_star_id=star.id,
                processed_star_id=processed_star.id)
        session.add(processed_star_association)

    session.commit()
