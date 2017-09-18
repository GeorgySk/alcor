from typing import List

from sqlalchemy.orm.session import Session

from alcor.models import (Group,
                          Star)
from alcor.models.velocities_vs_magnitudes import Cloud
from . import (raw,
               lepine)


def process_stars_group(*,
                        stars: List[Star],
                        group: Group,
                        w_lepine_criterion: bool,
                        session: Session) -> None:
    if w_lepine_criterion:
        clouds = list(lepine.clouds(stars=stars,
                                    group=group))
        bins = lepine.bins(stars=stars,
                           group=group)
    else:
        clouds = [Cloud(group_id=group.id,
                        bolometric_magnitude=star.bolometric_magnitude,
                        u_velocity=star.u_velocity,
                        v_velocity=star.v_velocity,
                        w_velocity=star.w_velocity)
                  for star in stars]
        bins = raw.bins(stars=stars,
                        group=group)
    session.add_all(clouds)
    session.add_all(bins)
    session.commit()
