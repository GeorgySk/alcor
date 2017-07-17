import logging
from collections import Counter
from functools import partial
from itertools import filterfalse
from typing import List

from cassandra.cluster import Session as CassandraSession
from sqlalchemy.orm.session import Session

from alcor.cassandra_models import (CGroup,
                                    CStar)
from alcor.cassandra_models.eliminations import CStarsCounter
from alcor.models import Group, Star
from alcor.models.eliminations import StarsCounter
from alcor.services import (luminosity_function,
                            velocities,
                            velocities_vs_magnitudes)
from alcor.services.data_access import (insert,
                                        model_insert_statement)
from .sampling import check_elimination

logging.basicConfig(format='%(filename)s %(funcName)s '
                           '%(levelname)s: %(message)s',
                    level=logging.DEBUG)
logger = logging.getLogger(__name__)


def process_stars_group(*,
                        c_stars: List[CStar],
                        stars: List[Star],
                        c_group: CGroup,
                        group: Group,
                        filtration_method: str,
                        nullify_radial_velocity: bool,
                        w_luminosity_function: bool,
                        w_velocities_clouds: bool,
                        w_velocities_vs_magnitude: bool,
                        w_lepine_criterion: bool,
                        cassandra_session: CassandraSession,
                        session: Session) -> None:
    stars_count = len(c_stars)
    logger.info('Starting processing c_stars, '
                f'objects number: {stars_count}.')
    c_eliminations_counter = Counter()
    eliminations_counter = Counter()
    c_apply_elimination_criteria = partial(
        check_elimination,
        eliminations_counter=c_eliminations_counter,
        filtration_method=filtration_method)
    apply_elimination_criteria = partial(
        check_elimination,
        eliminations_counter=eliminations_counter,
        filtration_method=filtration_method)
    if filtration_method in {'restricted', 'full'}:
        c_stars = list(filterfalse(c_apply_elimination_criteria,
                                   c_stars))
        stars = list(filterfalse(apply_elimination_criteria,
                                 stars))

    c_counter = CStarsCounter(
        group_id=c_group.id,
        raw=stars_count,
        by_parallax=c_eliminations_counter['parallax'],
        by_declination=c_eliminations_counter['declination'],
        by_velocity=c_eliminations_counter['velocity'],
        by_proper_motion=c_eliminations_counter['proper_motion'],
        by_reduced_proper_motion=c_eliminations_counter[
            'reduced_proper_motion'],
        by_apparent_magnitude=c_eliminations_counter['apparent_magnitude'])
    counter = StarsCounter(
        group_id=group.id,
        raw=stars_count,
        by_parallax=eliminations_counter['parallax'],
        by_declination=eliminations_counter['declination'],
        by_velocity=eliminations_counter['velocity'],
        by_proper_motion=eliminations_counter['proper_motion'],
        by_reduced_proper_motion=eliminations_counter['reduced_proper_motion'],
        by_apparent_magnitude=eliminations_counter['apparent_magnitude'])

    statement = model_insert_statement(CStarsCounter)
    insert(instances=[c_counter],
           statement=statement,
           session=cassandra_session)

    session.add(counter)

    if nullify_radial_velocity:
        for c_star in c_stars:
            c_star.set_radial_velocity_to_zero()
        for star in stars:
            star.set_radial_velocity_to_zero()

    if w_luminosity_function:
        luminosity_function.process_stars_group(
            c_stars=c_stars,
            stars=stars,
            c_group=c_group,
            group=group,
            cassandra_session=cassandra_session,
            session=session)

    if w_velocities_clouds:
        velocities.process_stars_group(
            stars=c_stars,
            group=c_group,
            w_lepine_criterion=w_lepine_criterion,
            cassandra_session=cassandra_session,
            session=session)

    if w_velocities_vs_magnitude:
        velocities_vs_magnitudes.process_stars_group(
            c_stars=c_stars,
            stars=stars,
            c_group=c_group,
            group=group,
            w_lepine_criterion=w_lepine_criterion,
            cassandra_session=cassandra_session,
            session=session)
