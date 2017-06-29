from typing import List

from cassandra.cluster import Session

from alcor.models import (Group,
                          Star)
from alcor.models.luminosity_function import Point
from alcor.services.data_access import (insert,
                                        model_insert_statement)
from .utils import (FORTY_PARSEC_NORTHERN_HEMISPHERE_VOLUME,
                    OBSERVATIONAL_DATA_TRUSTED_BINS_OBJECT_COUNT,
                    generate_stars_bins,
                    points)

# TODO: find out the name for set of indexes
MAGIC_STARS_BINS_INDEXES = {16 - 1, 17 - 1, 18 - 1}  # Thanks, Python!


def process_stars_group_luminosity_function(*,
                                            stars: List[Star],
                                            group: Group,
                                            session: Session
                                            ) -> None:
    stars_bins = generate_stars_bins(stars)
    # TODO: find out the meaning of sum
    magic_bins_lengths_sum = sum(len(stars_bins[index])
                                 for index in MAGIC_STARS_BINS_INDEXES)
    normalization_factor = (FORTY_PARSEC_NORTHERN_HEMISPHERE_VOLUME
                            * magic_bins_lengths_sum
                            / OBSERVATIONAL_DATA_TRUSTED_BINS_OBJECT_COUNT)
    graph_points = points(bins=stars_bins,
                          group=group,
                          normalization_factor=normalization_factor)
    statement = model_insert_statement(Point)
    insert(instances=graph_points,
           statement=statement,
           session=session)
