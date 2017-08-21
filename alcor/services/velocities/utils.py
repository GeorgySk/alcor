from typing import Union

from alcor.models import Star
from alcor.models.velocities import (LepineCaseUVCloud,
                                     LepineCaseUWCloud,
                                     LepineCaseVWCloud)


def star_cloud(star: Star) -> Union[LepineCaseUVCloud,
                                    LepineCaseUWCloud,
                                    LepineCaseVWCloud]:
    highest_coordinate = max(abs(star.coordinate_x),
                             abs(star.coordinate_y),
                             abs(star.coordinate_z))
    if abs(star.coordinate_x) == highest_coordinate:
        yield LepineCaseVWCloud(group_id=star.group_id,
                                velocity_v=star.velocity_v,
                                velocity_w=star.velocity_w)
    elif abs(star.coordinate_y) == highest_coordinate:
        yield LepineCaseUWCloud(group_id=star.group_id,
                                velocity_u=star.velocity_u,
                                velocity_w=star.velocity_w)
    else:
        yield LepineCaseUVCloud(group_id=star.group_id,
                                velocity_u=star.velocity_u,
                                velocity_v=star.velocity_v)
