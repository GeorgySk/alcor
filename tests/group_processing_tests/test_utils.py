from typing import List

from alcor.models import Star
from alcor.services.group_processing.utils import copy_velocities


def test_copy_velocities(src_stars: List[Star],
                         dst_stars: List[Star]) -> None:
    copy_velocities(src_stars,
                    dst_stars)

    assert all(dst_star.u_velocity == src_star.u_velocity and
               dst_star.v_velocity == src_star.v_velocity and
               dst_star.w_velocity == src_star.w_velocity
               for src_star, dst_star in zip(src_stars, dst_stars))
