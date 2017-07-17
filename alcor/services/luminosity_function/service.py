from typing import List

from sqlalchemy.orm.session import Session

from alcor.models import (Group,
                          Star)
from .utils import (FORTY_PARSEC_NORTHERN_HEMISPHERE_VOLUME,
                    OBSERVATIONAL_DATA_TRUSTED_BINS_OBJECT_COUNT,
                    generate_stars_bins,
                    points)

# TODO: find out the name for set of indexes
MAGIC_STARS_BINS_INDEXES = {16, 17, 18}


def process_stars_group(*,
                        stars: List[Star],
                        group: Group,
                        session: Session) -> None:
    stars_bins = generate_stars_bins(stars)
    # TODO: find out the meaning of sum
    magic_bins_lengths_sum = sum(len(stars_bins[index])
                                 for index in MAGIC_STARS_BINS_INDEXES)
    normalization_factor = (FORTY_PARSEC_NORTHERN_HEMISPHERE_VOLUME
                            * magic_bins_lengths_sum
                            / OBSERVATIONAL_DATA_TRUSTED_BINS_OBJECT_COUNT)
    graph_points = points(stars_bins=stars_bins,
                          group=group,
                          normalization_factor=normalization_factor)
    session.add_all(graph_points)
    session.commit()
