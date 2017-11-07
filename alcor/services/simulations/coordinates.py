import copy
from math import (sqrt,
                  pi,
                  sin,
                  cos,
                  asin,
                  acos,
                  atan,
                  copysign)

from alcor.models import Star


PC_PER_KPC = 1e3


# More info on conversions at: https://goo.gl/e7uiiZ
# ngp - North Galactic Pole
# theta - see BK angle at the link above
# kappa - km/sec in one au/year:
# transverse_velocity = KAPPA * proper_motion * distance
def modify(star: Star,
           solar_galactocentric_distance: float,
           ngp_declination: float = 0.478,
           theta: float = 2.147,
           ngp_right_ascension: float = 3.35,
           kappa: float = 4.74
           ) -> Star:
    new_star = copy.deepcopy(star)

    r_cylindrical_coordinate = new_star.r_cylindric_coordinate
    th_cylindrical_coordinate = new_star.th_cylindric_coordinate
    z_coordinate = new_star.z_coordinate
    u_velocity = new_star.u_velocity
    v_velocity = new_star.v_velocity
    w_velocity = new_star.w_velocity

    distance_plane_projection = (opposite_triangle_side(
        solar_galactocentric_distance,
        r_cylindrical_coordinate,
        th_cylindrical_coordinate))
    distance = sqrt(distance_plane_projection ** 2 + z_coordinate ** 2)

    galactic_longitude = triangle_angle(
        adjacent=solar_galactocentric_distance,
        other_adjacent=distance_plane_projection,
        opposite=r_cylindrical_coordinate)
    # Unfolding from 0-180 to 0-360
    if th_cylindrical_coordinate > pi:
        galactic_longitude = 2. * pi - galactic_longitude

    # TODO: or use arctan2
    galactic_latitude = atan(abs(z_coordinate / distance_plane_projection))
    galactic_latitude = copysign(galactic_latitude, z_coordinate)

    sin_longitude = sin(galactic_longitude)
    cos_longitude = cos(galactic_longitude)
    sin_latitude = sin(galactic_latitude)
    cos_latitude = cos(galactic_latitude)

    velocity_by_prop_motion = 1. / (kappa * distance * PC_PER_KPC)

    # TODO: find out if we need to divide by `cos_latitude`
    proper_motion_in_longitude = (velocity_by_prop_motion / cos_latitude) * (
        -u_velocity * sin_longitude + v_velocity * cos_longitude)
    proper_motion_in_latitude = velocity_by_prop_motion * (
        -u_velocity * cos_longitude * sin_latitude
        - v_velocity * sin_latitude * sin_longitude
        + w_velocity * cos_latitude)
    radial_velocity = (
        u_velocity * cos_latitude * cos_longitude
        + v_velocity * cos_latitude * sin_latitude
        + w_velocity * sin_latitude)

    # TODO: find out if we need to multiply l-component by cos(b)
    proper_motion = sqrt(proper_motion_in_longitude ** 2
                         + proper_motion_in_latitude ** 2)

    declination = (asin(sin(ngp_declination) * sin_latitude
                        + cos(ngp_declination) * cos_latitude
                          * cos(theta - galactic_longitude)))

    # TODO: give better names
    # These variables are for conversion from galactic to equatorial
    # coordinates. More info at the link above
    xs = cos_latitude * sin(theta - galactic_longitude) / cos(declination)
    xc = ((cos(ngp_declination) * sin_latitude
           - sin(ngp_declination) * cos_latitude
             * cos(theta - galactic_longitude))
          / cos(declination))

    # TODO: find out what is going on here
    if xs >= 0.:
        if xc >= 0.:
            right_ascension = asin(xs) + ngp_right_ascension
        else:
            right_ascension = acos(xc) + ngp_right_ascension
    else:
        if xc >= 0.:
            right_ascension = 2. * pi + asin(xs) + ngp_right_ascension
        else:
            right_ascension = pi - asin(xs) + ngp_right_ascension

    if right_ascension > 2. * pi:
        right_ascension -= 2. * pi

    new_star.distance = distance
    new_star.galactic_longitude = galactic_longitude
    new_star.galactic_latitude = galactic_latitude
    new_star.proper_motion_component_l = proper_motion_in_longitude
    new_star.proper_motion_component_b = proper_motion_in_latitude
    new_star.radial_velocity = radial_velocity
    new_star.proper_motion = proper_motion
    new_star.declination = declination
    new_star.right_ascension = right_ascension

    return new_star


def opposite_triangle_side(adjacent: float,
                           other_adjacent: float,
                           enclosed_angle: float) -> float:
    return sqrt(adjacent ** 2 + other_adjacent ** 2
                - 2. * adjacent * other_adjacent * cos(enclosed_angle))


def triangle_angle(adjacent: float,
                   other_adjacent: float,
                   opposite: float) -> float:
    return acos((adjacent ** 2 + other_adjacent ** 2 - opposite ** 2)
                / (2. * adjacent * other_adjacent))
