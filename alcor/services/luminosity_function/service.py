from typing import List

from cassandra.cluster import Session as CassandraSession
from sqlalchemy.orm.session import Session

from alcor.cassandra_models import (CGroup,
                                    CStar)
from alcor.cassandra_models.luminosity_function import CPoint
from alcor.models import Group, Star
from alcor.models.luminosity_function import Point
from alcor.services.data_access import (insert,
                                        model_insert_statement)
from .utils import (FORTY_PARSEC_NORTHERN_HEMISPHERE_VOLUME,
                    OBSERVATIONAL_DATA_TRUSTED_BINS_OBJECT_COUNT,
                    generate_stars_bins,
                    points)

# TODO: find out the name for set of indexes
MAGIC_STARS_BINS_INDEXES = {16, 17, 18}


def process_stars_group(*,
                        c_stars: List[CStar],
                        stars: List[Star],
                        c_group: CGroup,
                        group: Group,
                        cassandra_session: CassandraSession,
                        session: Session
                        ) -> None:
    c_stars_bins = generate_stars_bins(c_stars)
    stars_bins = generate_stars_bins(stars)
    # TODO: find out the meaning of sum
    с_magic_bins_lengths_sum = sum(len(c_stars_bins[index])
                                   for index in MAGIC_STARS_BINS_INDEXES)
    magic_bins_lengths_sum = sum(len(stars_bins[index])
                                 for index in MAGIC_STARS_BINS_INDEXES)
    с_normalization_factor = (FORTY_PARSEC_NORTHERN_HEMISPHERE_VOLUME
                              * с_magic_bins_lengths_sum
                              / OBSERVATIONAL_DATA_TRUSTED_BINS_OBJECT_COUNT)
    normalization_factor = (FORTY_PARSEC_NORTHERN_HEMISPHERE_VOLUME
                            * magic_bins_lengths_sum
                            / OBSERVATIONAL_DATA_TRUSTED_BINS_OBJECT_COUNT)
    с_graph_points = list(points(stars_bins=c_stars_bins,
                                 group=c_group,
                                 normalization_factor=с_normalization_factor,
                                 cls=CPoint))
    graph_points = list(points(stars_bins=stars_bins,
                               group=group,
                               normalization_factor=normalization_factor,
                               cls=Point))

    statement = model_insert_statement(CPoint)
    insert(instances=с_graph_points,
           statement=statement,
           session=cassandra_session)

    session.add_all(graph_points)
    session.commit()
