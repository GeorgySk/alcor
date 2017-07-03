from datetime import datetime
import logging
from typing import List

from sqlalchemy.orm.session import Session

from alcor.models import (Group,
                          Star)
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
                   session: Session) -> None:
    # TODO: add other options here
    if unprocessed_groups:
        groups = fetch_unprocessed_groups(session=session)
    for group in groups:
        stars = fetch_stars(group=group,
                            session=session)
        update_group(group,
                     session=session)
        process_stars_group(stars=stars,
                            group=group,
                            filtration_method=filtration_method,
                            nullify_radial_velocity=nullify_radial_velocity,
                            w_luminosity_function=w_luminosity_function,
                            w_velocities_clouds=w_velocities_clouds,
                            w_velocities_vs_magnitude=w_velocities_vs_magnitude,
                            w_lepine_criterion=w_lepine_criterion,
                            session=session)


def update_group(group: Group,
                 *,
                 session: Session) -> None:
    current_date_time = datetime.now()
    query = (session.query(Group)
             .filter(Group.id == group.id))
    query.update({Group.processed: True,
                  Group.updated_timestamp: current_date_time})


def fetch_unprocessed_groups(*,
                             session: Session) -> List[Group]:
    query = (session.query(Group)
             .filter(Group.processed.is_(False)))
    return query.all()


def fetch_stars(*,
                group: Group,
                session: Session
                ) -> List[Star]:
    query = (session.query(Star)
             .filter(Star.group_id == group.id))
    return query.all()
