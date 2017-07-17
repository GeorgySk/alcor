from itertools import chain
from typing import List

from cassandra.cluster import Session as CassandraSession
from sqlalchemy.orm.session import Session

from alcor.cassandra_models import (CGroup,
                                    CStar)
from alcor.cassandra_models.velocities_vs_magnitudes import (CCloud,
                                                             CLepineCaseUCloud,
                                                             CLepineCaseVCloud,
                                                             CLepineCaseWCloud,
                                                             CBin,
                                                             CLepineCaseUBin,
                                                             CLepineCaseVBin,
                                                             CLepineCaseWBin)
from alcor.models import (Group,
                          Star)
from alcor.models.velocities_vs_magnitudes import (Cloud,
                                                   LepineCaseUCloud,
                                                   LepineCaseVCloud,
                                                   LepineCaseWCloud,
                                                   Bin,
                                                   LepineCaseUBin,
                                                   LepineCaseVBin,
                                                   LepineCaseWBin)
from alcor.services.data_access import (insert,
                                        model_insert_statement)
from .utils import (generate_clouds,
                    generate_stars_bins,
                    generate_u_bins,
                    generate_v_bins,
                    generate_w_bins,
                    raw_stars_bins,
                    generate_bins)


def process_stars_group(*,
                        c_stars: List[CStar],
                        stars: List[Star],
                        c_group: CGroup,
                        group: Group,
                        w_lepine_criterion: bool,
                        cassandra_session: CassandraSession,
                        session: Session
                        ) -> None:
    if w_lepine_criterion:
        c_u_clouds, c_v_clouds, c_w_clouds = generate_clouds(
            stars=c_stars,
            group=c_group,
            u_cls=CLepineCaseUCloud,
            v_cls=CLepineCaseVCloud,
            w_cls=CLepineCaseWCloud)
        u_clouds, v_clouds, w_clouds = generate_clouds(stars=stars,
                                                       group=group,
                                                       u_cls=LepineCaseUCloud,
                                                       v_cls=LepineCaseVCloud,
                                                       w_cls=LepineCaseWCloud)
        c_clouds = {CLepineCaseUCloud: c_u_clouds,
                    CLepineCaseVCloud: c_v_clouds,
                    CLepineCaseWCloud: c_w_clouds}
        clouds = u_clouds + v_clouds + w_clouds
        u_stars_bins, v_stars_bins, w_stars_bins = generate_stars_bins(c_stars)
        c_u_bins = generate_u_bins(stars_bins=u_stars_bins,
                                   group=c_group,
                                   cls=CLepineCaseUBin)
        c_v_bins = generate_v_bins(stars_bins=v_stars_bins,
                                   group=c_group,
                                   cls=CLepineCaseVBin)
        c_w_bins = generate_w_bins(stars_bins=w_stars_bins,
                                   group=c_group,
                                   cls=CLepineCaseWBin)
        c_bins = {CLepineCaseUBin: c_u_bins,
                  CLepineCaseVBin: c_v_bins,
                  CLepineCaseWBin: c_w_bins}
        u_bins = generate_u_bins(stars_bins=u_stars_bins,
                                 group=group,
                                 cls=LepineCaseUBin)
        v_bins = generate_v_bins(stars_bins=v_stars_bins,
                                 group=group,
                                 cls=LepineCaseVBin)
        w_bins = generate_w_bins(stars_bins=w_stars_bins,
                                 group=group,
                                 cls=LepineCaseWBin)
        bins = chain(u_bins, v_bins, w_bins)
    else:
        c_clouds = {CCloud: [CCloud(
            group_id=c_group.id,
            bolometric_magnitude=star.bolometric_magnitude,
            velocity_u=star.velocity_u,
            velocity_v=star.velocity_v,
            velocity_w=star.velocity_w)
                             for star in c_stars]}
        clouds = [Cloud(group_id=group.id,
                        bolometric_magnitude=star.bolometric_magnitude,
                        velocity_u=star.velocity_u,
                        velocity_v=star.velocity_v,
                        velocity_w=star.velocity_w)
                  for star in stars]
        c_stars_bins = raw_stars_bins(c_stars)
        stars_bins = raw_stars_bins(stars)
        c_bins = {CBin: generate_bins(stars_bins=c_stars_bins,
                                      group=c_group,
                                      cls=CBin)}
        bins = generate_bins(stars_bins=stars_bins,
                             group=group,
                             cls=Bin)

    for model, instances in c_clouds.items():
        statement = model_insert_statement(model)
        insert(instances=instances,
               statement=statement,
               session=cassandra_session)

    for model, instances in c_bins.items():
        statement = model_insert_statement(model)
        insert(instances=instances,
               statement=statement,
               session=cassandra_session)

    session.add_all(clouds)
    session.add_all(bins)
    session.commit()
