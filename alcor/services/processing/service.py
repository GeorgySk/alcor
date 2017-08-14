import logging
import uuid
from typing import List

from sqlalchemy.orm.session import Session

from alcor.models import Group
from alcor.services.group_processing import process_stars_group

logger = logging.getLogger(__name__)


def run_processing(*,
                   filtration_method: str,
                   nullify_radial_velocity: bool,
                   w_luminosity_function: bool,
                   w_velocities_clouds: bool,
                   w_velocities_vs_magnitude: bool,
                   w_lepine_criterion: bool,
                   last_groups_count: int,
                   unprocessed_groups: bool,
                   id_groups: uuid.uuid4,
                   session: Session) -> None:
    # TODO: add fetching by id
    if unprocessed_groups:
        groups = fetch_unprocessed_groups(session=session)
    elif last_groups_count:
        groups = fetch_last_groups(count=last_groups_count,
                                   session=session)
    elif id_groups:
        groups = fetch_groups_by_id(id_list=id_groups,
                                    session=session)
    for group in groups:
        process_stars_group(
            group=group,
            filtration_method=filtration_method,
            nullify_radial_velocity=nullify_radial_velocity,
            w_luminosity_function=w_luminosity_function,
            w_velocities_clouds=w_velocities_clouds,
            w_velocities_vs_magnitude=w_velocities_vs_magnitude,
            w_lepine_criterion=w_lepine_criterion,
            session=session)


# TODO: move this to reading
def fetch_unprocessed_groups(*,
                             session: Session) -> List[Group]:
    query = (session.query(Group)
             .filter(Group.original_unprocessed_group_id.is_(None)))
    return query.all()


def fetch_last_groups(*,
                      count: int,
                      session: Session) -> List[Group]:
    query = (session.query(Group).
             order_by(Group.updated_timestamp.desc()).limit(count))
    return query.all()


def fetch_groups_by_id(*,
                       id_list: uuid.uuid4,
                       session: Session) -> List[Group]:
    query = (session.query(Group)
             .filter(Group.id.is_(id_list)))
    return query.all()
