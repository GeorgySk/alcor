from typing import Union

from alcor.models import Star
from alcor.models.velocities import (LepineCaseUVCloud,
                                     LepineCaseUWCloud,
                                     LepineCaseVWCloud)


def star_cloud(star: Star) -> Union[LepineCaseUVCloud,
                                    LepineCaseUWCloud,
                                    LepineCaseVWCloud]:
    max_coordinates_modulus = star.max_coordinates_modulus

    if abs(star.x_coordinate) == max_coordinates_modulus:
        return LepineCaseVWCloud(group_id=star.group_id,
                                 v_velocity=star.v_velocity,
                                 w_velocity=star.w_velocity)
    elif abs(star.y_coordinate) == max_coordinates_modulus:
        return LepineCaseUWCloud(group_id=star.group_id,
                                 u_velocity=star.u_velocity,
                                 w_velocity=star.w_velocity)
    else:
        return LepineCaseUVCloud(group_id=star.group_id,
                                 u_velocity=star.u_velocity,
                                 v_velocity=star.v_velocity)
