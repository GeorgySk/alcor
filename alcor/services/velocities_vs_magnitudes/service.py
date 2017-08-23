from itertools import chain
from typing import List

from sqlalchemy.orm.session import Session

from alcor.models import (Group,
                          Star)
from alcor.models.velocities_vs_magnitudes import Cloud
from .utils import (stars_clouds,
                    lepine_stars_bins,
                    generate_u_bins,
                    generate_v_bins,
                    generate_w_bins,
                    raw_stars_bins,
                    generate_bins)


def process_stars_group(*,
                        stars: List[Star],
                        group: Group,
                        w_lepine_criterion: bool,
                        session: Session) -> None:
    if w_lepine_criterion:
        clouds = list(stars_clouds(stars=stars,
                                   group=group))
        u_stars_bins, v_stars_bins, w_stars_bins = lepine_stars_bins(stars)
        u_bins = generate_u_bins(stars_bins=u_stars_bins,
                                 group=group)
        v_bins = generate_v_bins(stars_bins=v_stars_bins,
                                 group=group)
        w_bins = generate_w_bins(stars_bins=w_stars_bins,
                                 group=group)
        bins = chain(u_bins, v_bins, w_bins)
    else:
        clouds = [Cloud(group_id=group.id,
                        bolometric_magnitude=star.bolometric_magnitude,
                        u_velocity=star.u_velocity,
                        v_velocity=star.v_velocity,
                        w_velocity=star.w_velocity)
                  for star in stars]
        stars_bins = raw_stars_bins(stars)
        bins = generate_bins(stars_bins=stars_bins,
                             group=group)
    session.add_all(clouds)
    session.add_all(bins)
    session.commit()
