from typing import List
from uuid import UUID

from cassandra.cluster import Session

from alcor.models import Star
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


def process_stars_group_velocities_vs_magnitudes(*,
                                                 stars: List[Star],
                                                 group_id: UUID,
                                                 lepine_criterion: bool,
                                                 session: Session
                                                 ) -> None:
    if lepine_criterion:
        u_clouds, v_clouds, w_clouds = generate_clouds(stars=stars,
                                                       group_id=group_id)
        clouds = {LepineCaseUCloud: u_clouds,
                  LepineCaseVCloud: v_clouds,
                  LepineCaseWCloud: w_clouds}
        u_stars_bins, v_stars_bins, w_stars_bins = generate_stars_bins(stars)
        u_bins = generate_u_bins(stars_bins=u_stars_bins,
                                 group_id=group_id)
        v_bins = generate_v_bins(stars_bins=v_stars_bins,
                                 group_id=group_id)
        w_bins = generate_w_bins(stars_bins=w_stars_bins,
                                 group_id=group_id)
        bins = {LepineCaseUBin: u_bins,
                LepineCaseVBin: v_bins,
                LepineCaseWBin: w_bins}
    else:
        clouds = {Cloud: [Cloud(group_id=group_id,
                                bolometric_magnitude=star.bolometric_magnitude,
                                velocity_u=star.velocity_u,
                                velocity_v=star.velocity_v,
                                velocity_w=star.velocity_w)
                          for star in stars]}
        stars_bins = raw_stars_bins(stars)
        bins = {Bin: generate_bins(stars_bins=stars_bins,
                                   group_id=group_id)}

    for model, instances in clouds.items():
        statement = model_insert_statement(model)
        insert(instances=instances,
               statement=statement,
               session=session)

    for model, instances in bins.items():
        statement = model_insert_statement(model)
        insert(instances=instances,
               statement=statement,
               session=session)
