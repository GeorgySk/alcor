from typing import List

from cassandra.cluster import Session

from alcor.models import (Group,
                          Star)
from alcor.models.velocities import (Cloud,
                                     LepineCaseUVCloud,
                                     LepineCaseUWCloud,
                                     LepineCaseVWCloud)
from alcor.services.data_access import (insert,
                                        model_insert_statement)
from .utils import generate_clouds


def process_stars_group_velocities_clouds(*,
                                          stars: List[Star],
                                          group: Group,
                                          lepine_criterion: bool,
                                          session: Session
                                          ) -> None:
    if lepine_criterion:
        uv_clouds, uw_clouds, vw_clouds = generate_clouds(stars)
        clouds = {LepineCaseUVCloud: uv_clouds,
                  LepineCaseUWCloud: uw_clouds,
                  LepineCaseVWCloud: vw_clouds}
    else:
        clouds = {Cloud: [Cloud(group_id=group.id,
                                velocity_u=star.velocity_u,
                                velocity_v=star.velocity_v,
                                velocity_w=star.velocity_w)
                          for star in stars]}
    for model, instances in clouds.items():
        statement = model_insert_statement(model)
        insert(instances=instances,
               statement=statement,
               session=session)
