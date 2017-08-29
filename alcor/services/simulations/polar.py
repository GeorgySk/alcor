from typing import List

import numpy as np

from alcor.models.star import Star


# TODO: this belongs to generating only sphere stars
def generate_polar_coordinates(stars: List[Star],
                               sector_radius: float,
                               sector_covering_angle: float,
                               max_sector_raius: float,
                               min_sector_radius: float,
                               scale_length: float,
                               solar_galactocentric_distance: float) -> None:
    squared_sector_radius = sector_radius ** 2
    squared_max_sector_radius = max_sector_raius ** 2
    squared_min_sector_radius = min_sector_radius ** 2
    squared_radii_difference = (squared_max_sector_radius
                                - squared_min_sector_radius)
    sector_diameter = max_sector_raius - min_sector_radius

    for star in stars:
        while True:
            star.th_cylindric_coordinate = (sector_covering_angle
                                            * np.random.rand()
                                            - sector_covering_angle / 2.)
            if star.th_cylindric_coordinate < 0.:
                star.th_cylindric_coordinate += 2. * np.pi

            while True:
                random_valid_radius = min_sector_radius + (sector_diameter
                                                           * np.random.rand())
                zzy = 0.16 * np.random.rand()
                zzr = np.exp(-random_valid_radius / scale_length)
                if zzy <= zzr:
                    break

            zz = (random_valid_radius - min_sector_radius) / sector_diameter
            xx = squared_min_sector_radius + squared_radii_difference * zz
            star.r_cylindric_coordinate = np.sqrt(xx)

            xc = (star.r_cylindric_coordinate
                  * np.cos(star.th_cylindric_coordinate))
            yc = (star.r_cylindric_coordinate
                  * np.sin(star.th_cylindric_coordinate))
            dist = (xc - solar_galactocentric_distance) ** 2 + yc ** 2

            if 0.0000015 <= dist <= squared_sector_radius:
                break

        # TODO: find out the meaning
        star.mysterious_x = xc
        star.mysterious_y = yc
