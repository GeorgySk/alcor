from typing import Union

from alcor.models import Star
from alcor.models.velocities import (LepineCaseUVCloud,
                                     LepineCaseUWCloud,
                                     LepineCaseVWCloud)


def star_cloud(star: Star) -> Union[LepineCaseUVCloud,
                                    LepineCaseUWCloud,
                                    LepineCaseVWCloud]:
    if abs(star.coordinate_x) == star.max_coordinates_modulus:
        return LepineCaseVWCloud(group_id=star.group_id,
                                 velocity_v=star.velocity_v,
                                 velocity_w=star.velocity_w)
    elif abs(star.coordinate_y) == star.max_coordinates_modulus:
        return LepineCaseUWCloud(group_id=star.group_id,
                                 velocity_u=star.velocity_u,
                                 velocity_w=star.velocity_w)
    else:
        return LepineCaseUVCloud(group_id=star.group_id,
                                 velocity_u=star.velocity_u,
                                 velocity_v=star.velocity_v)
