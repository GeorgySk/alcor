from typing import List

from sqlalchemy.orm.session import Session

from alcor.models import (Group,
                          Star)
from alcor.models.velocities import Cloud
from .utils import generate_clouds


def process_stars_group(*,
                        stars: List[Star],
                        group: Group,
                        w_lepine_criterion: bool,
                        session: Session
                        ) -> None:
    if w_lepine_criterion:
        clouds = generate_clouds(stars)
    else:
        clouds = (Cloud(group_id=group.id,
                        velocity_u=star.velocity_u,
                        velocity_v=star.velocity_v,
                        velocity_w=star.velocity_w)
                  for star in stars)
    session.add_all(clouds)
    session.commit()
