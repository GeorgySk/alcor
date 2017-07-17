from typing import List

from cassandra.cluster import Session as CassandraSession
from sqlalchemy.orm.session import Session

from alcor.cassandra_models import (CGroup,
                                    CStar)
from alcor.cassandra_models.velocities import (CCloud,
                                               CLepineCaseUVCloud,
                                               CLepineCaseUWCloud,
                                               CLepineCaseVWCloud)
from alcor.models.velocities import (Cloud,
                                     LepineCaseUVCloud,
                                     LepineCaseUWCloud,
                                     LepineCaseVWCloud)
from alcor.services.data_access import (insert,
                                        model_insert_statement)
from .utils import generate_clouds


def process_stars_group(*,
                        stars: List[CStar],
                        group: CGroup,
                        w_lepine_criterion: bool,
                        cassandra_session: CassandraSession,
                        session: Session
                        ) -> None:
    if w_lepine_criterion:
        c_uv_clouds, c_uw_clouds, c_vw_clouds = generate_clouds(
            stars,
            uv_cls=CLepineCaseUVCloud,
            uw_cls=CLepineCaseUWCloud,
            vw_cls=CLepineCaseVWCloud)
        uv_clouds, uw_clouds, vw_clouds = generate_clouds(
            stars,
            uv_cls=LepineCaseUVCloud,
            uw_cls=LepineCaseUWCloud,
            vw_cls=LepineCaseVWCloud)
        с_clouds = {CLepineCaseUVCloud: c_uv_clouds,
                    CLepineCaseUWCloud: c_uw_clouds,
                    CLepineCaseVWCloud: c_vw_clouds}
        clouds = uv_clouds + uw_clouds + vw_clouds
    else:
        с_clouds = {CCloud: [CCloud(group_id=group.id,
                                    velocity_u=star.velocity_u,
                                    velocity_v=star.velocity_v,
                                    velocity_w=star.velocity_w)
                             for star in stars]}
        clouds = [Cloud(group_id=group.id,
                        velocity_u=star.velocity_u,
                        velocity_v=star.velocity_v,
                        velocity_w=star.velocity_w)
                  for star in stars]
    for model, instances in с_clouds.items():
        statement = model_insert_statement(model)
        insert(instances=instances,
               statement=statement,
               session=cassandra_session)
    session.add_all(clouds)
    session.commit()
