from typing import List

import numpy as np

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
                        ) -> None:
    for star in stars:
        if star.disk_belonging == GalacticDiskEnum.thin:
            u_velocity_dispersion = thin_disk_u_velocity_dispersion
            v_velocity_dispersion = thin_disk_v_velocity_dispersion
            w_velocity_dispersion = thin_disk_w_velocity_dispersion
        else:
            u_velocity_dispersion = thick_disk_u_velocity_dispersion
            v_velocity_dispersion = thick_disk_v_velocity_dispersion
            w_velocity_dispersion = thick_disk_w_velocity_dispersion
        # TODO: find out what is going on here
        uop = uom(u_peculiar_solar_velocity,
                  star.r_cylindric_coordinate,
                  star.th_cylindric_coordinate,
                  solar_galactocentric_distance,
                  oort_a_const,
                  oort_b_const)
        star.u_velocity = uop + u_velocity_dispersion * np.random.normal()

        vop = vom(v_peculiar_solar_velocity,
                  star.r_cylindric_coordinate,
                  star.th_cylindric_coordinate,
                  solar_galactocentric_distance,
                  oort_a_const,
                  oort_b_const)
        star.v_velocity = (vop + v_velocity_dispersion * np.random.normal()
                           - u_velocity_dispersion ** 2 / 120.)

        star.w_velocity = (w_peculiar_solar_velocity
                           + w_velocity_dispersion * np.random.normal())


def uom(u_peculiar_solar_velocity: float,
        radius: float,
        theta: float,
        solar_galactocentric_distance: float,
        oort_a_const: float,
        oort_b_const: float) -> float:
    return (u_peculiar_solar_velocity
            + ((3. - (2. * radius) / solar_galactocentric_distance)
               * oort_a_const - oort_b_const) * radius * np.sin(theta))


def vom(v_peculiar_solar_velocity: float,
        radius: float,
        theta: float,
        solar_galactocentric_distance: float,
        oort_a_const: float,
        oort_b_const: float) -> float:
    return (v_peculiar_solar_velocity
            + ((3. - (2. * radius) / solar_galactocentric_distance)
               * oort_a_const - oort_b_const) * radius * np.cos(theta)
            - (oort_a_const - oort_b_const) * solar_galactocentric_distance)
