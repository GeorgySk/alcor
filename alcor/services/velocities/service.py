from typing import List

from sqlalchemy.orm.session import Session

from alcor.models import (Group,
                          Star)
from alcor.models.velocities import Cloud
from .utils import star_cloud


def process_stars_group(*,
                        stars: List[Star],
                        group: Group,
                        w_lepine_criterion: bool,
                        session: Session
                        ) -> None:
    if w_lepine_criterion:
        clouds = map(star_cloud, stars)
    else:
        clouds = (Cloud(group_id=group.id,
                        u_velocity=star.u_velocity,
                        v_velocity=star.v_velocity,
                        w_velocity=star.w_velocity)
                  for star in stars)
    session.add_all(clouds)
    session.commit()
