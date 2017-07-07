import logging
import uuid
from collections import Counter
from functools import partial
from itertools import filterfalse
from math import (pi,
                  cos)
from typing import (List,
                    Tuple)

from sqlalchemy.orm.session import Session

from alcor.models import (Group,
                          Star)
from alcor.models.eliminations import StarsCounter
from alcor.models.simulation import Parameter
from alcor.services import (luminosity_function,
                            velocities,
                            velocities_vs_magnitudes)
from alcor.services.data_access.reading import fetch
from .sampling import check_elimination

logging.basicConfig(format='%(filename)s %(funcName)s '
                           '%(levelname)s: %(message)s',
                    level=logging.DEBUG)
logger = logging.getLogger(__name__)


DELTA_LATITUDE = 2.64 * pi / 180.0


def process_stars_group(*,
                        stars: List[Star],
                        group: Group,
                        filtration_method: str,
                        nullify_radial_velocity: bool,
                        w_luminosity_function: bool,
                        w_velocities_clouds: bool,
                        w_velocities_vs_magnitude: bool,
                        w_lepine_criterion: bool,
                        no_overlapping_cones: bool,
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

    # TODO: put this in correct place (where?)
    if no_overlapping_cones:
        stars = eliminate_stars_lying_in_prev_cones(group_id=group.id,
                                                    stars=stars,
                                                    session=session)

    original_unprocessed_group_id = group.id
    processed_group_id = uuid.uuid4()
    processed_group = Group(
        id=processed_group_id,
        original_unprocessed_group_id=original_unprocessed_group_id,
        processed=True)
    insert_groups_statement = model_insert_statement(Group)
    insert(instances=[processed_group],
           statement=insert_groups_statement,
           session=session)

    for star in stars:
        star.group_id = processed_group_id
        star.id = uuid.uuid4()
    insert_stars_statement = model_insert_statement(Star)
    insert(instances=stars,
           statement=insert_stars_statement,
           session=session)


def eliminate_stars_lying_in_prev_cones(group_id: uuid.uuid4,
                                        stars: List[Star],
                                        session: Session) -> List[Star]:
    (min_longitude,
     max_longitude,
     min_latitude,
     max_latitude) = get_cone_angles_ranges(group_id=group_id,
                                            session=session)

    (min_longitudes,
     max_longitudes,
     min_latitudes,
     max_latitudes) = get_overlapping_cone_angles_ranges(
        min_longitude=min_longitude,
        max_longitude=max_longitude,
        min_latitude=min_latitude,
        max_latitude=min_latitude,
        session=session)

    if min_longitudes:
        stars = [star
                 for star in stars
                 for (min_longitude,
                      max_longitude,
                      min_latitude,
                      max_latitude) in zip(min_longitudes,
                                           max_longitudes,
                                           min_latitudes,
                                           max_latitudes)
                 if not (min_longitude < star.galactic_longitude
                         < max_longitude
                         and min_latitude < star.galactic_latitude
                         < max_latitude)]
    logger.info(f'Number of stars left in the group after checking for '
                f'overlapping: {len(stars)}')
    return stars


def get_overlapping_cone_angles_ranges(min_longitude: float,
                                       max_longitude: float,
                                       min_latitude: float,
                                       max_latitude: float,
                                       session: Session
                                       ) -> Tuple[List[float], ...]:
    processed_groups = fetch_processed_groups(session=session)

    min_longitudes = []
    max_longitudes = []
    min_latitudes = []
    max_latitudes = []
    for group in processed_groups:
        (processed_group_min_longitude,
         processed_group_max_longitude,
         processed_group_min_latitude,
         processed_group_max_latitude) = get_cone_angles_ranges(
            group_id=group.original_unprocessed_group_id,
            session=session)
        # Cross-zero overlap cases:
        if (min_longitude < (processed_group_max_longitude - 2 * pi)
                and min_longitude < 0):
            processed_group_min_longitude -= 2 * pi
            processed_group_max_longitude -= 2 * pi
        if (processed_group_min_longitude < (max_longitude - 2 * pi)
                and processed_group_min_longitude < 0):
            min_longitude -= 2 * pi
            max_longitude -= 2 * pi

        # TODO: use this: https://stackoverflow.com/questions/2953967/built-in-function-for-computing-overlap-in-python
        longitude_overlapping = ((processed_group_min_longitude
                                  <= min_longitude
                                  <= processed_group_max_longitude)
                                 or (processed_group_min_longitude
                                     <= max_longitude
                                     <= processed_group_max_longitude)
                                 or (min_longitude
                                     <= processed_group_min_longitude
                                     <= max_longitude)
                                 or (min_longitude
                                     <= processed_group_max_longitude
                                     <= max_longitude))
        latitude_overlapping = ((processed_group_min_latitude
                                 <= min_latitude
                                 <= processed_group_max_latitude)
                                or (processed_group_min_latitude
                                    <= max_latitude
                                    <= processed_group_max_latitude)
                                or (min_latitude
                                    <= processed_group_min_latitude
                                    <= max_latitude)
                                or (min_latitude
                                    <= processed_group_max_latitude
                                    <= max_latitude))
        if longitude_overlapping and latitude_overlapping:
            min_longitudes.append(processed_group_min_longitude)
            max_longitudes.append(processed_group_max_longitude)
            min_latitudes.append(processed_group_min_latitude)
            max_latitudes.append(processed_group_max_latitude)
    return (min_longitudes,
            max_longitudes,
            min_latitudes,
            max_latitudes)


def fetch_processed_groups(*,
                           session: Session
                           ) -> List[Group]:
    query = (Group.objects
             .filter(Group.processed == True)
             .limit(None))
    records = fetch(query=query,
                    session=session)
    return [Group(**record)
            for record in records]


def get_cone_angles_ranges(group_id: uuid.uuid4,
                           session: Session) -> Tuple[float, ...]:
    simulation_parameters = fetch_model_by_group_id(model=Parameter,
                                                    id_list=[group_id],
                                                    session=session)
    # TODO: all angles should be recorded in radians when reading csv!
    for parameter in simulation_parameters:
        if parameter.name == 'longitude':
            height_longitude = float(parameter.value) * pi / 180
        if parameter.name == 'latitude':
            height_latitude = float(parameter.value) * pi / 180

    delta_longitude = DELTA_LATITUDE / cos(height_latitude)

    min_longitude = height_longitude - delta_longitude / 2
    max_longitude = height_longitude + delta_longitude / 2

    min_latitude = height_latitude - DELTA_LATITUDE / 2
    max_latitude = height_latitude + DELTA_LATITUDE / 2

    return (min_longitude,
            max_longitude,
            min_latitude,
            max_latitude)


def fetch_model_by_group_id(*,
                            session: Session,
                            model: Model,
                            id_list: List[uuid.uuid4]
                            ) -> List[Model]:
    # TODO: next 3 statements should be placed in separate function fetch-all
    query = (model.objects.all().limit(None))
    result_set = fetch(query=query,
                       session=session)
    all_records = [model(**result)
                   for result in result_set]
    filtered_records = [record
                        for record in all_records
                        if record.group_id in id_list]
    return filtered_records
