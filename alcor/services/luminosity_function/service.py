from typing import List

from sqlalchemy.orm.session import Session

from alcor.models import (Group,
                          Star)
from .utils import (stars_bins_normalization_factor,
                    pack_stars,
                    graph_points)


def process_stars_group(*,
                        stars: List[Star],
                        group: Group,
                        session: Session) -> None:
    stars_bins = pack_stars(stars)
    normalization_factor = stars_bins_normalization_factor(stars_bins)
    stars_counts = map(len, stars_bins)
    points = graph_points(stars_counts=stars_counts,
                          group=group,
                          normalization_factor=normalization_factor)
    session.add_all(points)
    session.commit()
