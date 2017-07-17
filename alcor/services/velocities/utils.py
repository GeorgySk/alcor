from typing import (Tuple,
                    List)

from alcor.models import Star
from alcor.models.velocities import (LepineCaseUVCloud,
                                     LepineCaseUWCloud,
                                     LepineCaseVWCloud)


def generate_clouds(stars: List[Star],
                    ) -> Tuple[List[LepineCaseUVCloud],
                               List[LepineCaseUWCloud],
                               List[LepineCaseVWCloud]]:
    for star in stars:
        highest_coordinate = max(star.coordinate_x,
                                 star.coordinate_y,
                                 star.coordinate_z)
        if star.coordinate_x == highest_coordinate:
            yield LepineCaseVWCloud(group_id=star.group_id,
                                    velocity_v=star.velocity_v,
                                    velocity_w=star.velocity_w)
        elif star.coordinate_y == highest_coordinate:
            yield LepineCaseUWCloud(group_id=star.group_id,
                                    velocity_u=star.velocity_u,
                                    velocity_w=star.velocity_w)
        else:
            yield LepineCaseUVCloud(group_id=star.group_id,
                                    velocity_u=star.velocity_u,
                                    velocity_v=star.velocity_v)
