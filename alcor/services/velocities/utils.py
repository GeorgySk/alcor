from typing import (Tuple,
                    List)

from alcor.cassandra_models import CStar
from alcor.cassandra_models.velocities import (CLepineCaseUVCloud,
                                               CLepineCaseUWCloud,
                                               CLepineCaseVWCloud)


def generate_clouds(stars: List[CStar],
                    uv_cls,
                    uw_cls,
                    vw_cls
                    ) -> Tuple[List[CLepineCaseUVCloud],
                               List[CLepineCaseUWCloud],
                               List[CLepineCaseVWCloud]]:
    uv_clouds = []
    uw_clouds = []
    vw_clouds = []
    for star in stars:
        highest_coordinate = max(star.coordinate_x,
                                 star.coordinate_y,
                                 star.coordinate_z)
        if star.coordinate_x == highest_coordinate:
            vw_clouds.append(vw_cls(group_id=star.group_id,
                                    velocity_v=star.velocity_v,
                                    velocity_w=star.velocity_w))
        elif star.coordinate_y == highest_coordinate:
            uw_clouds.append(uw_cls(group_id=star.group_id,
                                    velocity_u=star.velocity_u,
                                    velocity_w=star.velocity_w))
        else:
            uv_clouds.append(uv_cls(group_id=star.group_id,
                                    velocity_u=star.velocity_u,
                                    velocity_v=star.velocity_v))
    return uv_clouds, uw_clouds, vw_clouds
