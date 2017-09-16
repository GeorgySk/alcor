import uuid
from collections import Counter
from functools import partial
from itertools import filterfalse
from typing import (Iterator,
                    List)

from sqlalchemy.orm.session import Session

from alcor.models import (Group,
                          Star,
                          eliminations,
                          associations)
from alcor.services import (luminosity_function,
                            velocities,
                            velocities_vs_magnitudes)
from alcor.services.data_access import fetch_group_stars
from . import elimination
from .utils import copy_velocities


def process(*,
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

    eliminations_counter = Counter()
    if filtration_method in {'restricted', 'full'}:
        is_eliminated = partial(elimination.check,
                                eliminations_counter=eliminations_counter,
                                filtration_method=filtration_method)
        stars = list(filterfalse(is_eliminated, stars))

    counter = eliminations.StarsCounter(group_id=group.id,
                                        raw=stars_count,
                                        **eliminations_counter)

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
    processed_group = Group(id=processed_group_id,
                            original_id=original_id)
    session.add(processed_group)

    processed_stars = [Star(group_id=processed_group_id)
                       for _ in range(len(stars))]
    if nullify_radial_velocity:
        copy_velocities(stars, processed_stars)
    session.add_all(processed_stars)

    # We are flushing to set ids
    # more at
    # https://stackoverflow.com/questions/4201455/sqlalchemy-whats-the-difference-between-flush-and-commit
    session.flush(processed_stars)

    stars_associations = processed_stars_associations(stars, processed_stars)
    session.add_all(stars_associations)

    session.commit()


def processed_stars_associations(original_stars: List[Star],
                                 processed_stars: List[Star]
                                 ) -> Iterator[associations.ProcessedStars]:
    for original_star, processed_star in zip(original_stars,
                                             processed_stars):
        yield associations.ProcessedStars(
            original_id=original_star.id,
            processed_id=processed_star.id)
