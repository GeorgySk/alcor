import logging
import uuid
from typing import Optional

from sqlalchemy.orm.session import Session

from alcor.services.data_access import (fetch_unprocessed_groups,
                                        fetch_last_groups,
                                        fetch_group_by_id)
from alcor.services import stars_group

logger = logging.getLogger(__name__)


def run_processing(*,
                   filtration_method: str,
                   nullify_radial_velocity: bool,
                   w_luminosity_function: bool,
                   w_velocities_clouds: bool,
                   w_velocities_vs_magnitude: bool,
                   w_lepine_criterion: bool,
                   last_groups_count: Optional[int],
                   unprocessed_groups: bool,
                   group_id: Optional[uuid.UUID],
                   session: Session) -> None:
    if unprocessed_groups:
        groups = fetch_unprocessed_groups(session=session)
    elif last_groups_count:
        groups = fetch_last_groups(limit=last_groups_count,
                                   session=session)
    else:
        groups = [fetch_group_by_id(group_id=group_id,
                                    session=session)]
    for group in groups:
        stars_group.process(
            group=group,
            filtration_method=filtration_method,
            nullify_radial_velocity=nullify_radial_velocity,
            w_luminosity_function=w_luminosity_function,
            w_velocities_clouds=w_velocities_clouds,
            w_velocities_vs_magnitude=w_velocities_vs_magnitude,
            w_lepine_criterion=w_lepine_criterion,
            session=session)
