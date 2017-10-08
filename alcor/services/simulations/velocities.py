import copy
from functools import partial
from typing import List

import numpy as np
from math import (sin,
                  cos)

from alcor.models import Star


def modify(*,
           stars: List[Star],
           u_velocity_dispersion: float,
           v_velocity_dispersion: float,
           w_velocity_dispersion: float,
           u_peculiar_solar_velocity: float = -11.,
           v_peculiar_solar_velocity: float = -12.,
           w_peculiar_solar_velocity: float = -7.,
           solar_galactocentric_distance: float,
           oort_a_const: float,
           oort_b_const: float
           ) -> List[Star]:
    update_velocities = partial(
        update_star_velocities,
        oort_a_const=oort_a_const,
        oort_b_const=oort_b_const,
        solar_galactocentric_distance=solar_galactocentric_distance,
        u_peculiar_solar_velocity=u_peculiar_solar_velocity,
        v_peculiar_solar_velocity=v_peculiar_solar_velocity,
        w_peculiar_solar_velocity=w_peculiar_solar_velocity,
        u_velocity_dispersion=u_velocity_dispersion,
        v_velocity_dispersion=v_velocity_dispersion,
        w_velocity_dispersion=w_velocity_dispersion)

    yield from map(update_velocities, stars)


def update_star_velocities(star: Star,
                           *,
                           oort_a_const: float,
                           oort_b_const: float,
                           solar_galactocentric_distance: float,
                           u_peculiar_solar_velocity: float,
                           v_peculiar_solar_velocity: float,
                           w_peculiar_solar_velocity: float,
                           u_velocity_dispersion: float,
                           v_velocity_dispersion: float,
                           w_velocity_dispersion: float) -> Star:
    uop = uom(u_peculiar_solar_velocity,
              star.r_cylindric_coordinate,
              star.th_cylindric_coordinate,
              solar_galactocentric_distance,
              oort_a_const,
              oort_b_const)
    u_velocity = uop + u_velocity_dispersion * np.random.normal()
    vop = vom(v_peculiar_solar_velocity,
              star.r_cylindric_coordinate,
              star.th_cylindric_coordinate,
              solar_galactocentric_distance,
              oort_a_const,
              oort_b_const)
    v_velocity = (vop + v_velocity_dispersion * np.random.normal()
                  - u_velocity_dispersion ** 2 / 120.)
    w_velocity = (w_peculiar_solar_velocity
                  + w_velocity_dispersion * np.random.normal())

    return star.modify(u_velocity=u_velocity,
                       v_velocity=v_velocity,
                       w_velocity=w_velocity)


def uom(u_peculiar_solar_velocity: float,
        radius: float,
        theta: float,
        solar_galactocentric_distance: float,
        oort_a_const: float,
        oort_b_const: float) -> float:
    return (u_peculiar_solar_velocity
            + ((3. - (2. * radius) / solar_galactocentric_distance)
               * oort_a_const - oort_b_const) * radius * sin(theta))


def vom(v_peculiar_solar_velocity: float,
        radius: float,
        theta: float,
        solar_galactocentric_distance: float,
        oort_a_const: float,
        oort_b_const: float) -> float:
    return (v_peculiar_solar_velocity
            + ((3. - (2. * radius) / solar_galactocentric_distance)
               * oort_a_const - oort_b_const) * radius * cos(theta)
            - (oort_a_const - oort_b_const) * solar_galactocentric_distance)
