from math import (sqrt,
                  pi,
                  sin,
                  cos,
                  asin,
                  acos,
                  atan)
from typing import List

from alcor.models.star import Star


PC_PER_KPC = 1e3


# More info on conversions at: https://goo.gl/e7uiiZ
# ngp - North Galactic Pole
# theta - see BK angle at the link above
# kappa - km/sec in one au/year:
# transverse_velocity = KAPPA * proper_motion * distance
def calculate_coordinates(stars: List[Star],
                          solar_galactocentric_distance: float,
                          ngp_declination: float = 0.478,
                          theta: float = 2.147,
                          ngp_right_ascension: float = 3.35,
                          kappa: float = 4.74
                          ) -> None:
    for star in stars:
        distance_plane_projection = (opposite_triangle_side(
            solar_galactocentric_distance,
            star.r_cylindric_coordinate,
            star.th_cylindric_coordinate))
        star.distance = sqrt(distance_plane_projection ** 2
                             + star.z_coordinate ** 2)

        # TODO: implement function
        star.galactic_longitude = acos((solar_galactocentric_distance ** 2
                                        + distance_plane_projection ** 2
                                        - star.r_cylindric_coordinate ** 2)
                                       / (2. * distance_plane_projection
                                          * solar_galactocentric_distance))
        # Unfolding from 0-180 to 0-360
        if star.th_cylindric_coordinate > pi:
            star.galactic_longitude = 2. * pi - star.galactic_longitude

        # TODO: or use arctan2
        star.galactic_latitude = atan(abs(star.z_coordinate
                                          / distance_plane_projection))
        if star.z_coordinate < 0.:
            star.galactic_latitude = -star.galactic_latitude

        sin_longitude = sin(star.galactic_longitude)
        cos_longitude = cos(star.galactic_longitude)
        sin_latitude = sin(star.galactic_latitude)
        cos_latitude = cos(star.galactic_latitude)

        velocity_by_prop_motion = 1. / (kappa * star.distance * PC_PER_KPC)

        # TODO: find out if we need to divide by `cos_latitude`
        star.proper_motion_component_l = velocity_by_prop_motion * (
            -star.u_velocity * sin_longitude
            + star.v_velocity * cos_longitude) / cos_latitude
        star.proper_motion_component_b = velocity_by_prop_motion * (
            -star.u_velocity * cos_longitude * sin_latitude
            - star.v_velocity * sin_latitude * sin_longitude
            + star.w_velocity * cos_latitude)
        star.radial_velocity = (
            (cos_latitude * cos_longitude * star.u_velocity)
            + (cos_latitude * sin_latitude * star.v_velocity)
            + (sin_latitude * star.w_velocity))

        # TODO: find out if we need to multiply l-component by cos(b)
        star.proper_motion = sqrt(star.proper_motion_component_l ** 2
                                  + star.proper_motion_component_b ** 2)

        star.declination = (asin(sin(ngp_declination) * sin_latitude
                                 + cos(ngp_declination) * cos_latitude
                                   * cos(theta - star.galactic_longitude)))

        # TODO: give better names
        # These variables are for conversion from galactic to equatorial
        # coordinates. More info at the link above
        xs = ((cos_latitude * sin(theta - star.galactic_longitude))
              / cos(star.declination))
        xc = ((cos(ngp_declination) * sin_latitude
               - sin(ngp_declination) * cos_latitude
                 * cos(theta - star.galactic_longitude))
              / cos(star.declination))

        # TODO: find out what is going on here
        if xs >= 0.:
            if xc >= 0.:
                star.right_ascension = asin(xs) + ngp_right_ascension
            else:
                star.right_ascension = acos(xc) + ngp_right_ascension
        else:
            if xc >= 0.:
                star.right_ascension = 2. * pi + asin(xs) + ngp_right_ascension
            else:
                star.right_ascension = pi - asin(xs) + ngp_right_ascension
                
        if star.right_ascension > 2. * pi:
            star.right_ascension -= 2. * pi


def opposite_triangle_side(adjacent: float,
                           other_adjacent: float,
                           enclosed_angle: float) -> float:
    return sqrt(adjacent ** 2 + other_adjacent ** 2
                - 2. * adjacent * other_adjacent * cos(enclosed_angle))
