from copy import deepcopy
from functools import partial
from typing import List

import numpy as np
from math import (sin,
                  cos)

from alcor.models.star import (Star,
                               GalacticDiskEnum)


def generate_velocities(*,
                        stars: List[Star],
                        u_peculiar_solar_velocity: float = -11.,
                        v_peculiar_solar_velocity: float = -12.,
                        w_peculiar_solar_velocity: float = -7.,
                        thin_disk_u_velocity_dispersion: float = 32.4,
                        thin_disk_v_velocity_dispersion: float = 23.,
                        thin_disk_w_velocity_dispersion: float = 18.1,
                        thick_disk_u_velocity_dispersion: float = 50.,
                        thick_disk_v_velocity_dispersion: float = 56.,
                        thick_disk_w_velocity_dispersion: float = 34.,
                        solar_galactocentric_distance: float,
                        oort_a_const: float,
                        oort_b_const: float
                        ) -> List[Star]:
    thin_disk_stars = [star
                       for star in stars
                       if star.disk_belonging == GalacticDiskEnum.thin]
    thick_disk_stars = [star
                        for star in stars
                        if star.disk_belonging == GalacticDiskEnum.thick]

    update_velocities = partial(
        update_star_velocities,
        oort_a_const=oort_a_const,
        oort_b_const=oort_b_const,
        solar_galactocentric_distance=solar_galactocentric_distance,
        u_peculiar_solar_velocity=u_peculiar_solar_velocity,
        v_peculiar_solar_velocity=v_peculiar_solar_velocity,
        w_peculiar_solar_velocity=w_peculiar_solar_velocity)

    new_thin_disk_stars = [update_velocities(
            star,
            u_velocity_dispersion=thin_disk_u_velocity_dispersion,
            v_velocity_dispersion=thin_disk_v_velocity_dispersion,
            w_velocity_dispersion=thin_disk_w_velocity_dispersion)
        for star in thin_disk_stars]
    new_thick_disk_stars = [update_velocities(
            star,
            u_velocity_dispersion=thick_disk_u_velocity_dispersion,
            v_velocity_dispersion=thick_disk_v_velocity_dispersion,
            w_velocity_dispersion=thick_disk_w_velocity_dispersion)
        for star in thick_disk_stars]

    return new_thin_disk_stars + new_thick_disk_stars


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
    new_star = deepcopy(star)

    uop = uom(u_peculiar_solar_velocity,
              star.r_cylindric_coordinate,
              star.th_cylindric_coordinate,
              solar_galactocentric_distance,
              oort_a_const,
              oort_b_const)
    new_star.u_velocity = uop + u_velocity_dispersion * np.random.normal()
    vop = vom(v_peculiar_solar_velocity,
              star.r_cylindric_coordinate,
              star.th_cylindric_coordinate,
              solar_galactocentric_distance,
              oort_a_const,
              oort_b_const)
    new_star.v_velocity = (vop + v_velocity_dispersion * np.random.normal()
                           - u_velocity_dispersion ** 2 / 120.)
    new_star.w_velocity = (w_peculiar_solar_velocity
                           + w_velocity_dispersion * np.random.normal())

    return new_star


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
