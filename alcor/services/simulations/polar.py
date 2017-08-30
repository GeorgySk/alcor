from math import (pi,
                  sqrt,
                  exp,
                  cos,
                  sin)
from typing import List

import numpy as np

from alcor.models.star import Star


# TODO: this belongs to generating only sphere stars
def assign_polar_coordinates(stars: List[Star],
                             sector_radius: float,
                             sector_covering_angle: float,
                             max_sector_radius: float,
                             min_sector_radius: float,
                             scale_length: float,
                             solar_galactocentric_distance: float) -> None:
    squared_sector_radius = sector_radius ** 2
    squared_max_sector_radius = max_sector_radius ** 2
    squared_min_sector_radius = min_sector_radius ** 2
    squared_radii_difference = (squared_max_sector_radius
                                - squared_min_sector_radius)
    sector_diameter = max_sector_radius - min_sector_radius

    for star in stars:
        assign_coordinates(
            star=star,
            sector_covering_angle=sector_covering_angle,
            min_sector_radius=min_sector_radius,
            sector_diameter=sector_diameter,
            scale_length=scale_length,
            squared_min_sector_radius=squared_min_sector_radius,
            squared_radii_difference=squared_radii_difference,
            solar_galactocentric_distance=solar_galactocentric_distance,
            squared_sector_radius=squared_sector_radius)


def get_random_valid_radius(min_sector_radius: float,
                            sector_diameter: float,
                            scale_length: float) -> float:
    while True:
        random_valid_radius = min_sector_radius + (sector_diameter
                                                   * np.random.rand()[0])
        zzy = 0.16 * np.random.rand()
        zzr = exp(-random_valid_radius / scale_length)
        if zzy <= zzr:
            return random_valid_radius


def assign_coordinates(star: Star,
                       sector_covering_angle: float,
                       min_sector_radius: float,
                       sector_diameter: float,
                       scale_length: float,
                       squared_min_sector_radius: float,
                       squared_radii_difference: float,
                       solar_galactocentric_distance: float,
                       squared_sector_radius: float
                       ) -> None:
    while True:
        star.th_cylindric_coordinate = th_cylindric_coordinate(
            sector_covering_angle)

        random_valid_radius = get_random_valid_radius(
            min_sector_radius=min_sector_radius,
            sector_diameter=sector_diameter,
            scale_length=scale_length)

        zz = (random_valid_radius - min_sector_radius) / sector_diameter
        xx = squared_min_sector_radius + squared_radii_difference * zz
        star.r_cylindric_coordinate = sqrt(xx)

        xc = (star.r_cylindric_coordinate
              * cos(star.th_cylindric_coordinate))
        yc = (star.r_cylindric_coordinate
              * sin(star.th_cylindric_coordinate))
        dist = (xc - solar_galactocentric_distance) ** 2 + yc ** 2

        if 1.5e-6 <= dist <= squared_sector_radius:
            # TODO: find out the meaning
            star.mysterious_x = xc
            star.mysterious_y = yc

            return


def th_cylindric_coordinate(sector_covering_angle: float) -> float:
    res = (sector_covering_angle * np.random.rand()[0]
           - sector_covering_angle / 2.)
    if res < 0.:
        res += 2. * pi
    return res
