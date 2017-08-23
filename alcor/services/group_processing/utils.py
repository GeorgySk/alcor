from typing import List

from alcor.models import Star


def copy_velocities(src_stars: List[Star],
                    dst_stars: List[Star]) -> None:
    for src_star, dst_star in zip(src_stars, dst_stars):
        dst_star.u_velocity = src_star.u_velocity
        dst_star.v_velocity = src_star.v_velocity
        dst_star.w_velocity = src_star.w_velocity
