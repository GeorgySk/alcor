from typing import Union, Tuple, List

from alcor.models import Star
from alcor.models.velocities import (LepineCaseUVCloud,
                                     LepineCaseUWCloud,
                                     LepineCaseVWCloud)


def star_cloud(stars: List[Star]) -> Tuple[List[LepineCaseUVCloud],
                                           List[LepineCaseUWCloud],
                                           List[LepineCaseVWCloud]]:
    uv_cloud = []
    uw_cloud = []
    vw_cloud = []
    for star in stars:
        highest_coordinate = max(abs(star.coordinate_x),
                                 abs(star.coordinate_y),
                                 abs(star.coordinate_z))
        if abs(star.coordinate_x) == highest_coordinate:
            vw_cloud.append(LepineCaseVWCloud(group_id=star.group_id,
                                              velocity_v=star.velocity_v,
                                              velocity_w=star.velocity_w))
        elif abs(star.coordinate_y) == highest_coordinate:
            uw_cloud.append(LepineCaseUWCloud(group_id=star.group_id,
                                              velocity_u=star.velocity_u,
                                              velocity_w=star.velocity_w))
        else:
            uv_cloud.append(LepineCaseUVCloud(group_id=star.group_id,
                                              velocity_u=star.velocity_u,
                                              velocity_v=star.velocity_v))
    return uv_cloud, uw_cloud, vw_cloud
